import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Filters dbSNP by 1000 Genomes SNPs and subsets columns to values by GOAL CNV calling')
parser.add_argument('file', help='dbSNP file')
parser.add_argument('output', help='Output file')
args = parser.parse_args()
inFile = args.file
outFile = args.output

os.makedirs(os.path.dirname(outFile), exist_ok=True)
inCols = [
  "TEMP1",
  "CHROM",
  "START",
  "END",
  "rsID",
  "TEMP2",
  "TEMP3",
  "TEMP4",
  "TEMP5",
  "TEMP6",
  "TEMP7",
  "VARIANT_TYPE",
  "CALL_TYPE",
  "TEMP8",
  "TEMP9",
  "TEMP10",
  "TEMP11",
  "TEMP12",
  "TEMP13",
  "TEMP14",
  "TEMP15",
  "TEMP16",
  "TEMP17",
  "TEMP18",
  "TEMP19",
  "TEMP20",
]
useCols = [
  "CHROM",
  "START",
  "END",
  "rsID",
  "VARIANT_TYPE",
  "CALL_TYPE"
]

rowsProcessed = 0
chunkSize = 100000
for rowChunk in pd.read_csv(inFile, compression='gzip', header=None, sep='\t', comment='#', chunksize=chunkSize, names=inCols, usecols=useCols):
  # Filter to 1000 genomes SNPs
  rowChunk = rowChunk[(rowChunk['VARIANT_TYPE']=="single") & (rowChunk['CALL_TYPE'].str.contains("by-1000genomes"))]
  rowChunk.to_csv(outFile, header=False, sep='\t', index=False, mode='a', columns = [
    "CHROM",
    "START",
    "END",
    "rsID"
  ])

  rowsProcessed += chunkSize
  print("Processed " + str(rowsProcessed) + " lines")