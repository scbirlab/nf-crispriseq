#!/usr/bin/env python

"""Left-join stdin with first file on commandline, write to stdout."""

import pandas as pd
import sys

right_input = sys.argv[1]

try:
    sep2 = sys.argv[2]
except IndexError:
    sep2 = "\t"

right = pd.read_csv(right_input, sep=sep2)

for i, left_chunk in enumerate(pd.read_csv(sys.stdin, sep="\t", chunksize=10000)):

    if i == 0:
        common_cols = set(left_chunk.columns).intersection(set(right.columns))

        if len(common_cols) == 0:
            raise AttributeError(f'No shared columns!\nLeft: {", ".join(left_chunk.columns)}\nRight: {", ".join(right.columns)}')

        print("Joining on shared columns {}".format(", ".join(sorted(common_cols))),
             file=sys.stderr)
    
    (pd.merge(left_chunk, right, how='left')
     .to_csv(sys.stdout, 
             mode='w' if i == 0 else 'a',
             header=i == 0,
             index=False, 
             sep='\t'))