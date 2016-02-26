#! /usr/bin/env python

from __future__ import print_function
import csv
import argparse
import collections

# Parse command line
parser = argparse.ArgumentParser(description='Phasing comparison')
parser.add_argument('-t', '--table', metavar='svgeno.txt', required=True, dest='tenxTable', help='genotype table (required)')
args = parser.parse_args()

if args.tenxTable:
    with open(args.tenxTable) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        psID = ""
        total = collections.Counter()
        correctPhase = collections.Counter()
        haplotypeTotal = 0
        haplotypeCooccurrence = 0
        for row in reader:
            if ('genotype' not in row.keys()) or ((int(row['genotype']) == 1) and (int(row['calledgenotype']) == 1)):
                hGiven = None
                if '|' in row['haplotype']:
                    hGiven = [int(i) for i in row['haplotype'].split('|')]
                hCalled = [int(i) for i in row['calledhaplotype'].split('|')]
                if psID != row['phasedblockid']:
                    psID = row['phasedblockid']
                    if '/' in row['haplotype']:
                        hGiven = [int(i) for i in row['haplotype'].split('/')]
                else:
                    if hGiven is not None:
                        total[row['chr']] += 1
                        if ((hGiven[0] == hCalled[0]) and (oldHGiven[0] == oldHCalled[0])) or ((hGiven[0] == hCalled[1]) and (oldHGiven[0] == oldHCalled[1])):
                            correctPhase[row['chr']] += 1
                    haplotypeTotal += 1
                    if hCalled[0] == oldHCalled[0]:
                        haplotypeCooccurrence += 1
                oldHGiven = hGiven
                oldHCalled = hCalled
        for c in sorted(total.keys()):
            if total[c]:
                print("Phasing accuracy, ", c, ": ", float(correctPhase[c])/float(total[c]), " ( #n =", total[c], ")")
        genomictotal = sum(total.values())
        genomicCorrectPhase = sum(correctPhase.values())
        if genomictotal:
            print("Phasing accuracy: ", float(genomicCorrectPhase)/float(genomictotal), " ( #n =", genomictotal, ")")
        print("Haplotype co-occurrence: ", float(haplotypeCooccurrence)/float(haplotypeTotal), " ( #n =", haplotypeTotal, ")")




