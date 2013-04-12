// compilation:
// rdmd --build-only -I/path/to/BioD/ -O -release -inline read_bam.d

import bio.bam.reader, bio.bam.md.core;
import std.range, std.algorithm, std.parallelism, std.conv, std.stdio;

struct ReadStats {
  enum kMinimumMappingQuality = 1;

  ubyte mapping_quality;
  size_t n_mismatches;
  size_t n_insertions;
  size_t n_deletions;
  size_t length_insertions;
  size_t length_deletions;
  size_t qlen;

  this(BamRead read) {
    auto cigar = read.cigar;
    auto insertions = cigar.filter!(op => op.type == 'I');
    auto deletions  = cigar.filter!(op => op.type == 'D');
    auto soft_clips = cigar.filter!(op => op.type == 'S');

    qlen = read.sequence.length - reduce!"a+b"(0, soft_clips.map!(op => op.length));
    mapping_quality = read.mapping_quality;
    
    if (read.is_unmapped || mapping_quality < kMinimumMappingQuality)
      mapping_quality = 0;

    if (mapping_quality == 0)
      return;

    n_insertions = insertions.walkLength();
    n_deletions = deletions.walkLength();
    length_insertions = reduce!"a+b"(0, insertions.map!(op => op.length));
    length_deletions = reduce!"a+b"(0, deletions.map!(op => op.length));

    auto md_ops = mdOperations(to!string(read["MD"]));
    n_mismatches = md_ops.filter!(op => op.is_mismatch).walkLength();
  }
}

void print(ReadStats stats, size_t read_id, string description, string reference_name)
{
  auto mapped = stats.mapping_quality > 0;
  writeln(description,        '\t',     
          reference_name,     '\t', read_id,                  '\t',     
          mapped ? 1 : 0,     '\t', stats.mapping_quality,    '\t',     
          stats.n_insertions, '\t', stats.length_insertions,  '\t',     
          stats.n_deletions,  '\t', stats.length_deletions,   '\t',     
          stats.qlen,         '\t', stats.n_mismatches);
}

void main(string[] args) {
  auto pool = new TaskPool(16);
  scope(exit) pool.finish();

  if (args.length < 4) {
    stderr.writeln("Usage: ", args[0], " <input.bam> <description> <reference name>");
    return;
  }
  
  auto bam = new BamReader(args[1], pool);
  auto description = args[2];
  auto reference_name = args[3];

  writeln("sample",     '\t', 
          "ref",        '\t',     "rid",          '\t', 
          "mapped",     '\t',     "mapq",         '\t',
          "insertions", '\t',     "l_insertions", '\t', 
          "deletions",  '\t',     "l_deletions",  '\t',
          "rlen",       '\t',     "subst");
  
  auto id = 0U;
  foreach (read; bam.reads)
    ReadStats(read).print(++id, description, reference_name);
}
