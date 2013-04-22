import bio.sam.reader;
import bio.bam.writer;
import std.stdio, std.parallelism, std.math, std.algorithm, std.conv;

void main(string[] args) {
  if (args.length < 3) {
    stderr.writeln("usage: ./filter_unique <input.sam> <output.bam>");
    return;
  }

  auto tp = new TaskPool(8);
  scope(exit) tp.finish();
  auto src = new SamReader(args[1]);
  auto dst = new BamWriter(args[2], 1, tp);
  scope(exit) dst.finish();

  dst.writeSamHeader(src.header);
  dst.writeReferenceSequenceInfo(src.reference_sequences);

  uint total, mapped, unmapped, unique, not_unique;
  string last_id;
 
  auto reads = src.reads;
  while (!reads.empty) {
    auto read = reads.front;
    reads.popFront();
    
    if (read.is_unmapped) {
      total += 1;
      unmapped += 1;
      continue;
    }

    auto count = 1;
    while (!reads.empty && reads.front.name == read.name) {
      reads.popFront();
      ++count;
    }

    auto bwasw_equal = false;
    auto xs = read["XS"];
    auto as = read["AS"];
    if (xs.is_integer && as.is_integer)
      bwasw_equal = (xs.to!int() - as.to!int()).abs < 1;

    if (count == 1 && !read["XA"].is_string && !bwasw_equal) {
      dst.writeRecord(read);
      unique += 1;
    } else {
      not_unique += 1;
    }

    mapped += 1;
    total += 1;
  }

  stderr.writeln("Reads found:\t\t", total);
  stderr.writeln("Reads mapped:\t\t", mapped);
  stderr.writeln("  Uniquely mapped:\t", unique);
  stderr.writeln("  Not uniquely mapped:\t", not_unique);
  stderr.writeln("Not mapped:\t", unmapped);
}
