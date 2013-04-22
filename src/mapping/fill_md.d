import bio.bam.reader, bio.bam.writer;
import bio.core.fasta, bio.core.utils.format : write;
import std.array, std.stdio, std.conv, std.algorithm, std.range, std.parallelism, std.string;

__gshared string[string] refs;

BamRead fillmd(BamRead read) {
  auto new_read = read.dup;

  auto reference_info = read.reader.reference_sequences[read.ref_id];
  auto ref_seq = refs[reference_info.name];
  auto read_seq = read.sequence;

  auto nm = 0;
  auto ref_pos = read.position;
  auto read_pos = 0;
  auto md_match = 0;
  char[65536] md_buffer;
  char* md = md_buffer.ptr;

  foreach (op; read.cigar) {
    if (op.is_match_or_mismatch) {
      foreach (j; 0 .. op.length) {
        if (ref_pos + j >= ref_seq.length)
          break;
        auto ref_base = ref_seq[ref_pos + j];
        auto read_base = read_seq[read_pos + j];
        if ((ref_base == read_base && ref_base != 'N' && read_base != 'N') 
            || read_base == '=')
        {
          ++md_match;
        } else {
          md.write(md_match);
          md.write(ref_base);
          md_match = 0;
          nm += 1;
        }
      }
      if (ref_pos + op.length > ref_seq.length)
        break;
      ref_pos += op.length;
      read_pos += op.length;
    } else if (op.type == 'D') {
      md.write(md_match);
      md.write('^');
      md.write(ref_seq[ref_pos .. min(ref_pos + op.length, $)]);
      md_match = 0;
      if (ref_pos + op.length > ref_seq.length)
        break;
      nm += op.length;
      ref_pos += op.length;
    } else if (op.type == 'I' || op.type == 'S') {
      read_pos += op.length;
      if (op.type == 'I')
        nm += op.length;
    } else if (op.type == 'N') {
      ref_pos += op.length;
    }
  }

  md.write(md_match);

  new_read["NM"] = nm;
  new_read["MD"] = cast(string)md_buffer[0 .. md - md_buffer.ptr];
  return new_read;
}

void printUsage() {
  stderr.writeln("usage: fillmd <ref.fa> <input.bam> <output.bam>");
}

void main(string[] args) {
  if (args.length < 4) {
    printUsage();
    return;
  }
  foreach (reference; fastaRecords(args[1]))
    refs[reference.header.split()[0]] = reference.sequence;
  auto taskpool = new TaskPool(16);
  scope(exit) taskpool.finish();
  auto bam = new BamReader(args[2], taskpool);
  auto writer = new BamWriter(args[3], -1, taskpool);
  scope(exit) writer.finish();
  writer.writeSamHeader(bam.header);
  writer.writeReferenceSequenceInfo(bam.reference_sequences);
  foreach (read; taskpool.map!fillmd(bam.reads, 1024)) {
    writer.writeRecord(read);
  }
}
