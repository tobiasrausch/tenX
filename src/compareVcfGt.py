#! /usr/bin/env python

from __future__ import print_function
import banyan
import vcf
import argparse
import csv


# Parse command line
parser = argparse.ArgumentParser(description='VCF comparison')
parser.add_argument('-v', '--vcf', metavar='variants1.vcf', required=True, dest='vFile', help='input vcf file (required)')
parser.add_argument('-w', '--wvcf', metavar='variants2.vcf', required=True, dest='wFile', help='input vcf file (required)')
parser.add_argument('-s', '--sample', metavar='NA12878', required=True, dest='sample', help='sample name (required)')
parser.add_argument('-t', '--type', metavar='DEL', required=True, dest='svtype', help='SV type (required)')
args = parser.parse_args()

# Parse samples
sampleSet = []
if args.vFile:
    vcf1 = vcf.Reader(open(args.vFile), 'r', compressed=True) if args.vFile.endswith('.gz') else vcf.Reader(open(args.vFile), 'r', compressed=False)
    v1 = set(vcf1.samples)
    if args.wFile:
        vcf2 = vcf.Reader(open(args.wFile), 'r', compressed=True) if args.wFile.endswith('.gz') else vcf.Reader(open(args.wFile), 'r', compressed=False)
        sampleSet = v1.intersection(set(vcf2.samples))

if args.sample in sampleSet:
    sv = dict()
    if args.vFile:
        vcf_reader = vcf.Reader(open(args.vFile), 'r', compressed=True) if args.vFile.endswith('.gz') else vcf.Reader(open(args.vFile), 'r', compressed=False)
        for record in vcf_reader:
            svt = record.INFO['SVTYPE']
            if args.svtype != svt:
                continue
            svEnd = int(record.INFO['END'])
            size = svEnd - record.POS
            if svt == "INS":
                size = record.INFO['INSLEN']
            call = record.genotype(args.sample)
            if not sv.has_key(record.CHROM):
                sv[record.CHROM] = banyan.SortedDict(key_type=(int, int), alg=banyan.RED_BLACK_TREE, updator=banyan.OverlappingIntervalsUpdator)
            sv[record.CHROM][(record.POS, (record.POS+size))] = (record.ID, call['GT'], call.gt_type)
        if args.wFile:
            print("chr", "start", "end", "id", "size", "haplotype", "genotype", "calledhaplotype", "calledgenotype", sep="\t")
            vcf_reader = vcf.Reader(open(args.wFile), 'r', compressed=True) if args.wFile.endswith('.gz') else vcf.Reader(open(args.wFile), 'r', compressed=False)
            for record in vcf_reader:
                if record.CHROM in sv.keys():
                    try:
                        svt = record.INFO['SVTYPE']
                    except KeyError:
                        if len(record.ALT)>1:
                            continue
                        else:
                            if len(record.REF) < len(record.ALT[0]):
                                svt = "INS"
                            elif len(record.REF) > len(record.ALT[0]):
                                svt = "DEL"
                    if args.svtype != svt:
                        continue
                    try:
                        svEnd = int(record.INFO['END'])
                    except KeyError:
                        svEnd = record.POS + abs(len(record.REF) - len(record.ALT[0]))
                    call = record.genotype(args.sample)
                    for cStart, cEnd in sv[record.CHROM].overlap((record.POS, svEnd)):
                        cSvID, cHap, cGT = sv[record.CHROM][(cStart, cEnd)]
                        if (abs(record.POS - cStart)<25) and (abs(svEnd - cEnd)<25):
                            print(record.CHROM, cStart, cEnd, cSvID, (cEnd - cStart), call['GT'], call.gt_type, cHap, cGT, sep="\t")
                            break
