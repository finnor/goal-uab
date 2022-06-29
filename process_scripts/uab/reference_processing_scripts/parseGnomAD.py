import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Subsets gnomAD columns to values by GOAL pipeline')
parser.add_argument('file', help='gnomAD bcf')
parser.add_argument('output', help='Output file')
parser.add_argument('-c', '--chr', action='store_true', help='Add chr to chromosome names')
args = parser.parse_args()
inFile = args.file
outFile = args.output
addChr = args.chr

def getVariantInfoDict(variantInfo):
  keys = ["nhomalt", "AF", "AF_popmax"]
  variantInfo = list(map(lambda infoPairs: infoPairs if "=" in infoPairs else infoPairs + '=', variantInfo))
  variantInfo = dict(map(lambda infoPairs: tuple(infoPairs.split("=")), variantInfo))
  return { key: variantInfo.get(key, None) for key in keys }

os.makedirs(os.path.dirname(outFile), exist_ok=True)
header = True
inCols = [
  "#CHROM",
  "POS",
  "ID",
  "REF",
  "ALT",
  "QUAL",
  "FILTER",
  "INFO"
]
useCols = [
  "#CHROM",
  "POS",
  "REF",
  "ALT",
  "INFO",
]
dtypes = {
  '#CHROM': str, 'POS': int, 'REF': str, 'ALT': str, 'INFO': str
}

rowsProcessed = 0
chunkSize = 100000
for rowChunk in pd.read_csv(inFile, compression='gzip', header=0, sep='\t', comment='#', chunksize=chunkSize, names=inCols, usecols=useCols, dtype=dtypes):
  #Get homozygous count, allele frequency, and max population allele frequency from INFO
  info = rowChunk["INFO"].str.split(";").apply(getVariantInfoDict)
  rowChunk = pd.concat([rowChunk.drop(['INFO'], axis=1), info.apply(pd.Series)], axis=1)

  #Process AFs
  rowChunk["AF"] = rowChunk["AF"].apply(pd.to_numeric)
  rowChunk = rowChunk.fillna({"AF": 0.00e+00, "AF_popmax": "0.00e+00"})
  rowChunk["AF"] = rowChunk['AF'].map(lambda val: "{:.2e}".format(val))
  rowChunk = rowChunk.rename(columns={
    "nhomalt": "GNOMAD_HOM",
    "AF": "GNOMAD_AF",
    "AF_popmax": "AF_POPMAX",
  })

  #Add hg19 coordinates
  rowChunk = rowChunk.assign(GNOMAD_HG19_VARIANT = rowChunk["#CHROM"].astype(str) + "-" + rowChunk["POS"].astype(str) + "-" + rowChunk["REF"] + "-" + rowChunk["ALT"])

  #Convert chromosomes to chr syntax, if flagged
  if(addChr):
    rowChunk["#CHROM"] = rowChunk["#CHROM"].apply(lambda chr: "chr" + str(chr))

  rowChunk.to_csv(outFile, header=header, sep='\t', index=False, mode='a', columns = [
    "#CHROM",
    "POS",
    "REF",
    "ALT",
    "GNOMAD_HOM",
    "GNOMAD_AF",
    "AF_POPMAX",
    "GNOMAD_HG19_VARIANT"
  ])
  header = False

  rowsProcessed += chunkSize
  print("Processed " + str(rowsProcessed) + " lines")