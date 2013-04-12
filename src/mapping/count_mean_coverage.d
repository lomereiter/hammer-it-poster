import bio.bam.reader, bio.bam.pileup, bio.bam.baseinfo;
import std.stdio, std.parallelism, std.conv, std.range, std.algorithm;

void main(string[] args) {
  auto pool = new TaskPool(16);
  scope(exit) pool.finish();

  auto bam = new BamReader(args[1], pool);

  auto n_total_bases = 0U;
  auto n_total_length = 0U;
  auto n_total_nonzero_coverage = 0U;

  foreach (int i, ref_seq; bam.reference_sequences) {
    n_total_length += ref_seq.length;
    auto n_ref_bases = 0U;
    auto n_ref_nonzero_coverage = 0U;
    
    auto reads = bam.reference(i)[];
    foreach (column; makePileup(reads)) {
      n_total_bases += column.coverage;
      n_ref_bases += column.coverage;
      n_ref_nonzero_coverage += 1;
      n_total_nonzero_coverage += 1;
    }

    writeln(ref_seq.name, '\t',
            ref_seq.length, '\t',
            n_ref_nonzero_coverage, '\t',
            n_ref_bases.to!double() / ref_seq.length);
  }
}
