#!/usr/bin/env python
#
# FASTA Merge (c) 2011 Andrew Paprocki / MIT License
#
# The purpose of this script is to take a FASTA file with mtDNA HVR1+HVR2
# information and merge any CR SNPs found within a 23andMe, Inc. mtDNA
# export back into the FASTA file for more accurate analysis.
#
# The script takes as input two filenames. The first filename is a valid
# FASTA dump of a full mtDNA sequence. This script does not compare the
# input against the Cambridge reference, so any insertions/deletions must
# be removed prior to processing the file. The second filename is a valid
# 23andMe, Inc. export containing mtDNA ('MT' chromosome) SNP information.
# The file is tab delimited specifying rsid, chromosome, position, and
# genotype.
#
# Only CR SNPs (00575-16000) will be merged back into the original FASTA
# data and the resulting FASTA dump will be printed to standard output.
#
# Usage: $ dna-fasta-merge.py 123456-FASTA.fasta genome_Someone.txt
#
from sys import argv
import string

script, fasta, snpfile = argv

# Process FTDNA mtDNA FASTA dump

def chunk_fasta(l):
    return [c for c in l if c in "ACGNT" and l[0] != '<']

def read_fasta(x):
    data = open(x)
    header = data.readline().strip('\r\n')
    vals = list()
    map(vals.extend, [chunk_fasta(l) for l in data])
    return header, vals

header, mtdna = read_fasta(fasta)

# Process 23andMe, Inc. export, grabbing and converting 'MT' chromosome data

def read_snpfile(x):
    data = open(x)
    vals = list()
    map(vals.append, [l.strip('\r\n').split('\t') for l in data if l[0] != '#'])
    return dict([(int(l[2]), l[3]) for l in vals if l[1] == 'MT' and l[3] in "ACGT"])

snps = read_snpfile(snpfile)

# Loop over 23andMe SNPs and steal any in the coding region (575-16000)

for k, v in [(k, snps[k]) for k in sorted(snps.iterkeys()) if k >= 575 and k <= 16000]:
    # 23andMe Yoruba to Cambridge conversion
    # http://www.snpedia.com/index.php/MtDNA_Position_Conversions
    if k <= 309: k = k - 1
    elif k >= 3109 and k <= 16183: k = k - 1 - 1
    elif (k >= 312 and k <= 3108) or (k >= 16185 and k <= 16571): k = k - 2 - 1
    else: continue
    #if mtdna[k] != v: print k, mtdna[k-1:k+2], '=>', v
    mtdna[k] = v

def print_fasta():
    print header
    pos = 0
    while pos < len(mtdna):
        print string.joinfields(mtdna[pos:pos+80], '')
        pos = pos + 80

print_fasta()
