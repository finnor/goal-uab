import numpy as np
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Deduplicates repeat flags in repeat_regions bed')
parser.add_argument('file', help='BGZipped Repeat Regions')
parser.add_argument('output', help='Output file')
args = parser.parse_args()
inFile = args.file
outFile = args.output

def dedupRepeat(repeatFlags):
  return ",".join(np.sort(np.unique(np.array(repeatFlags.split(",")))))

os.makedirs(os.path.dirname(outFile), exist_ok=True)
inCols = [
  "CHROM",
  "START",
  "END",
  "REPEAT_TYPE"
]
dtypes = {
  'CHROM': str, 'START': int, 'END': int, 'REPEAT_TYPE': str
}

rowsProcessed = 0
chunkSize = 100000
for rowChunk in pd.read_csv(inFile, compression='gzip', header=None, sep='\t', chunksize=chunkSize, names=inCols, dtype=dtypes):
  #Remove duplicate values from repeat type
  rowChunk["REPEAT_TYPE"] = rowChunk["REPEAT_TYPE"].apply(dedupRepeat)

  rowChunk.to_csv(outFile, header=False, sep='\t', index=False, mode='a', columns = [
    "CHROM",
    "START",
    "END",
    "REPEAT_TYPE"
  ])

  rowsProcessed += chunkSize
  print("Processed " + str(rowsProcessed) + " lines")