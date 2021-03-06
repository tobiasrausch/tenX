#! /usr/bin/env python

from __future__ import print_function
import csv
import argparse
import numpy

# Parse command line
parser = argparse.ArgumentParser(description='Genotype comparison')
parser.add_argument('-t', '--table', metavar='svgeno.txt', required=True, dest='tenxTable', help='genotype table (required)')
parser.add_argument('-m', '--minsize', type=int, default=50, metavar='50', required=False, dest='minsize', help='min. SV size (optional)')
args = parser.parse_args()

geno = numpy.zeros((4, 4), dtype=numpy.int32)
if args.tenxTable:
    with open(args.tenxTable) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            if ('size' not in row.keys()) or (int(row['size']) > args.minsize):
                if row['calledgenotype'] == 'None':
                    row['calledgenotype'] = -1
                if row['genotype'] == 'None':
                    row['genotype'] = -1
                geno[int(row['genotype']) + 1, int(row['calledgenotype']) + 1] += 1

# Compute genotype concordances
for i in range(1, 4):
    if sum(geno[1:, i]) != 0:
        discordance = (float(sum(geno[1:, i])) - float(geno[i, i]))/float(sum(geno[1:, i]))
        discordance = numpy.round(discordance, decimals=4)
        print("Genotype " + str(i-1) + " discordance: ", discordance)
# FDR
denom = sum(geno[2:, 1]) + sum(geno[2:, 2]) + sum(geno[2:, 3])
fdr = 'NA'
if denom:
    fdr = float(sum(geno[2:, 1]))/float(denom)
print("FDR: ", fdr, " (#n=", denom, ")", sep="")
# FNR
denom = sum(geno[1, 2:]) + sum(geno[2, 2:]) + sum(geno[3, 2:])
fnr = 'NA'
if denom:
    fnr = float(sum(geno[1, 2:]))/float(denom)
print("FNR: ", fnr, " (#n=", denom, ")", sep="")

# Print genotype confusion matrix
print("GT\tNA", end="")
for i in range(3):
    print("\t" + str(i), end="")
print()
for i in range(4):
    for j in range(4):
        if j != 0:
            print("\t", end="")
        else:
            if i:
                print(str(i-1) + "\t", end="")
            else:
                print("NA\t", end="")
        print(geno[i, j], end="")
    print()


