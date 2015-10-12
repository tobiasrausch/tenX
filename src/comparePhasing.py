#! /usr/bin/env python

from __future__ import print_function
import csv
import argparse
import numpy

# Parse command line
parser = argparse.ArgumentParser(description='Phasing comparison')
parser.add_argument('-t', '--table', metavar='svgeno.txt', required=True, dest='tenxTable', help='genotype table (required)')
args = parser.parse_args()

geno = numpy.zeros((4, 4), dtype=numpy.int32)
if args.tenxTable:
    with open(args.tenxTable) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        psID = ""
        correctPhase = 0
        total = 0
        haplotypeTotal = 0
        haplotypeCooccurrence = 0
        for row in reader:
            if (int(row['genotype']) == 1) and (int(row['calledgenotype']) == 1):
                hGiven = None
                if '|' in row['haplotype']:
                    hGiven = [int(i) for i in row['haplotype'].split('|')]
                hCalled = [int(i) for i in row['calledhaplotype'].split('|')]
                if psID != row['phasedblockid']:
                    psID = row['phasedblockid']
                    oldHGiven = hGiven
                    oldHCalled = hCalled
                else:
                    if hGiven is not None:
                        total += 1
                        if ((hGiven[0] == hCalled[0]) and (oldHGiven[0] == oldHCalled[0])) or ((hGiven[0] == hCalled[1]) and (oldHGiven[0] == oldHCalled[1])):
                            correctPhase += 1
                    haplotypeTotal += 1
                    if hCalled[0] == oldHCalled[0]:
                        haplotypeCooccurrence +=1
        if total:
            print("Phasing accuracy: ", float(correctPhase)/float(total), " ( #n =", total, ")")
        print("Haplotype co-occurrence: ", float(haplotypeCooccurrence)/float(haplotypeTotal), " ( #n =", haplotypeTotal, ")")




