#!/usr/bin/env python

from collections import Counter
import gzip
import sys
import time

from bioino.fasta import FastaCollection
from tqdm.auto import tqdm

def parse_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            yield header, seq, plus, qual


if __name__ == "__main__":
    fasta_filename = sys.argv[1]
    reads_filename = sys.argv[2]
    
    guide_index = {
        seq.sequence.casefold(): seq.name
        for seq in FastaCollection.from_file(fasta_filename).sequences
    }
    lengths = sorted(set(len(s) for s in guide_index), reverse=True)
    variable_lengths = len(lengths) > 1

    counts = Counter()
    length_counts = Counter()
    n_total = n_matched = 0

    with gzip.open("matched.fastq.gz", "wt") as matched, gzip.open("unmatched.fastq.gz", "wt") as unmatched:
        for i, (header, seq, plus, qual) in enumerate(tqdm(
        parse_fastq(reads_filename), 
        desc=f"Streaming {reads_filename}",
    )):
            seq_lowercase = seq.casefold()
            hit = guide_index.get(seq_lowercase)
            if variable_lengths and hit is None:
                for l in lengths:  # get longest match
                    hit = guide_index.get(seq_lowercase[:l])
                    if hit is not None:
                        hit_len = l
                        break
            else:
                hit_len = len(seq) 
            if hit is not None:
                length_counts[hit_len] += 1
                counts[hit] += 1
                header = f"{header.split()[0]} {hit}"
                handle = matched
            else:
                handle = unmatched
            sequence_to_write = f"{header}\n{seq}\n{plus}\n{qual}"
            print(sequence_to_write, file=handle)

    n_total = i + 1
    n_matched = counts.total()
    percent_matched = 100. * n_matched / max(n_total, 1)
    with open("demux.log", "w") as log:
        print(
            f"""Total reads: {n_total}
            Matched: {n_matched} ({percent_matched:.1f}%)
            Unmatched: {n_total - n_matched}

            Lengths found:
            """,
            file=log,
        )
        for l, count in sorted(length_counts.items(), key=lambda x: int(x[0])):
            print(f"\t{l} nt: {count}", file=log)
