import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Subsets cancerhotspots columns to values by GOAL pipeline')
parser.add_argument('file', help='CancerHotspots maf')
parser.add_argument('output', help='Output file')
args = parser.parse_args()
inFile = args.file
outFile = args.output

os.makedirs(os.path.dirname(outFile), exist_ok=True)
header = True
inCols = [
  "Chromosome",
  "Start_Position",
  "End_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Hugo_Symbol",
  "HGVSp",
  "allele_freq",
  "oncotree_organtype",
  "oncotree_parent",
  "oncotree_detailed"
]
rowsProcessed = 0
chunkSize = 10000
for rowChunk in pd.read_csv(inFile, compression='gzip', header=0, sep='\t', quotechar='"', comment='#', chunksize=chunkSize, usecols=inCols):
  #Convert chromosomes to chr syntax
  rowChunk["Chromosome"] = rowChunk["Chromosome"].apply(lambda chr: "chr" + str(chr))

  rowChunk = rowChunk.rename(columns={
    "Chromosome": "#CHROM",
    "Start_Position": "FROM",
    "End_Position": "TO",
    "Reference_Allele": "OncoKB_REF",
    "Tumor_Seq_Allele2": "OncoKB_ALT",
    "Hugo_Symbol": "Gene",
    "HGVSp": "OncoKB_ProteinChange",
    "allele_freq": "OncoKB_AF",
    "oncotree_organtype": "OncoTree_Tissue",
    "oncotree_parent": "OncoTree_MainType",
    "oncotree_detailed": "OncoTree_Code"
  })
  rowChunk = rowChunk.assign(OncoKBHotspot = 'OncoKBHotspot')
  rowChunk.to_csv(outFile, header=header, sep='\t', index=False, mode='a', columns = [
    "#CHROM",
    "FROM",
    "TO",
    "OncoKB_REF",
    "OncoKB_ALT",
    "Gene",
    "OncoKB_ProteinChange",
    "OncoKB_AF",
    "OncoTree_Tissue",
    "OncoTree_MainType",
    "OncoTree_Code",
    "OncoKBHotspot"
  ])
  header = False

  rowsProcessed += chunkSize
  print("Processed " + str(rowsProcessed) + " lines")