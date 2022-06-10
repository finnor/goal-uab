import argparse
import csv

parser = argparse.ArgumentParser(description='Mutates the Illumina sample sheet to work for BCL convert')
parser.add_argument('original_sample_sheet', help='The original sample sheet')
parser.add_argument('override_cycles', help='String to insert for override cycles setting')
args = parser.parse_args()
originalSampleSheet = args.original_sample_sheet
overrideCycles = args.override_cycles

output = []
with open(originalSampleSheet) as sheetFile:
  csv_reader = csv.reader(sheetFile, delimiter=',')
  for row in csv_reader:
    if (row[0]=="[Settings]"):
      output.append(row)
      output.append(["OverrideCycles", overrideCycles])
    elif (row[0]=="Adapter"):
      row[0] = "AdapterRead1"
      output.append(row)
    else:
      output.append(row)

with open('NewSampleSheet.csv', 'w') as outfile:
  csv_writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
  for row in output:
    csv_writer.writerow(row)