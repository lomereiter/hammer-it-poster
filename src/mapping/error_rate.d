import bio.bam.reader, bio.bam.read;
import std.stdio, std.parallelism;

void main(string[] args) {
  auto pool = new TaskPool(4);
  scope(exit) pool.finish();
  auto bam = new BamReader(args[1], pool);

  uint[65536] mismatches_at, deletions_before, insertions_starting_at, total;

  foreach (read; bam.reads) {
    if (read.is_unmapped) continue;
    total[0 .. read.sequence_length] += 1;
    CigarOperation[4096] buf = void;
    int buf_pos = 0;
    foreach (op; read.extended_cigar) buf[buf_pos++] = op;
    auto cigar = buf[0 .. buf_pos];
    if (read.strand == '-') cigar.reverse;
    int pos = 0;
    foreach (op; cigar) {
      if      (op.type == 'X') ++mismatches_at[pos];
      else if (op.type == 'D') ++deletions_before[pos];
      else if (op.type == 'I') ++insertions_starting_at[pos];
      if (op.is_query_consuming) pos += op.length;
    }
  }

  writeln("offset\ttotal\tsubst\tdel\tins");
  foreach (size_t i, n; total)
    if (n > 0)
      writeln(i, "\t", n, "\t",
              mismatches_at[i], "\t",
              deletions_before[i], "\t",
              insertions_starting_at[i], "\t");
}
