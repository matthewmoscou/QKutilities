#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Predicts the allele of a KASP marker based on a de novo assembled transcriptome, genome, or other sequence resource (such as GBS markers).
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Reads FASTA sequence
	2. Reads a KASP marker database
    3. For every marker, determine if present in the sequence data set and return the allele
Note: Script could substantially improve based on optimization of search process
"""


## modules
import optparse
from optparse import OptionParser

## global variables
DNA = ['A', 'C', 'G', 'T']
cmp_DNA = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'M':'K', 'R':'Y', 'W':'S', 'S':'W', 'Y':'R', 'K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V', 'X':'X', 'N':'N'}
IUPAC_convert = {'AG':'R', 'GA':'R', 'CT':'Y', 'TC':'Y', 'CA':'M', 'AC':'M', 'TG':'K', 'GT':'K', 'TA':'W', 'AT':'W', 'CG':'S', 'GC':'S'}
IUPAC_SNP = {'R':['A', 'G'], 'Y':['C', 'T'], 'M':['A', 'C'], 'K':['G', 'T'], 'W':['A', 'T'], 'S':['C', 'G']}

## functions
# reverse complement
# input : DNA sequence
# output : reverse complement of said DNA sequence
def reverse_complement(orig_DNA):
    rev_comp_DNA = ''
    for index in range(len(orig_DNA)):
        rev_comp_DNA += cmp_DNA[orig_DNA[len(orig_DNA)-index-1]]
    return rev_comp_DNA


## OptionParser
# import arguments and options
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta", default="", help="FASTA file")
parser.add_option("-k", "--kasp", action="store", type="string", dest="kasp", default="", help="KASP database (see README for structure)")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default="", help="Output for primer matches")
(options, args) = parser.parse_args()


## import fasta

fasta_file = open(options.fasta, 'r')

ID_sequence = {}

ID = ''
line = fasta_file.readline()

while line:
	if len(line) > 0:
		if line[0] == '>':
			ID = line.split()[0][1:]
			ID_sequence[ID] = ''
		else:
			for sline in line.split():
				ID_sequence[ID] += sline.upper()

	line = fasta_file.readline()

fasta_file.close()

## import KASP file
KASP_marker_database = open(options.kasp, 'r')

marker_sequence = {}
truth = False

for line in KASP_marker_database.readlines():
    line = line.replace('\n', '')
    sline = line.split('\t')

    if truth:
        if len(sline) >= 14:
            if len(sline[12]) == 41:
                marker_sequence[sline[0]] = [sline[12][21:40], sline[14], sline[12][40], sline[13][40], [], [], []]
    else:
        truth = True

KASP_marker_database.close()


for marker in marker_sequence.keys():
    marker_rc = reverse_complement(marker_sequence[marker][0])

    for contig in ID_sequence.keys():
        if marker_sequence[marker][0] in ID_sequence[contig]:
            if marker_sequence[marker][1] in ID_sequence[contig]:
                reverse_primer = True
            elif reverse_complement(marker_sequence[marker][1]) in ID_sequence[contig]:
                reverse_primer = True
            else:
                reverse_primer = False

            if len(ID_sequence[contig]) > (ID_sequence[contig].index(marker_sequence[marker][0]) + len(marker_sequence[marker][0])):
                marker_sequence[marker][4].append(contig)
                marker_sequence[marker][5].append(ID_sequence[contig][ID_sequence[contig].index(marker_sequence[marker][0]) + len(marker_sequence[marker][0])])
                marker_sequence[marker][6].append(reverse_primer)
        elif marker_rc in ID_sequence[contig]:
            if marker_sequence[marker][1] in ID_sequence[contig]:
                reverse_primer = True
            elif reverse_complement(marker_sequence[marker][1]) in ID_sequence[contig]:
                reverse_primer = True
            else:
                reverse_primer = False

            if ID_sequence[contig].index(reverse_complement(marker_sequence[marker][0])) - 1 >= 0:
                marker_sequence[marker][4].append(contig)
                marker_sequence[marker][5].append(cmp_DNA[ID_sequence[contig][ID_sequence[contig].index(reverse_complement(marker_sequence[marker][0])) - 1]])
                marker_sequence[marker][6].append(reverse_primer)

output_file = open(options.output, 'w')

for marker in marker_sequence.keys():
    if len(marker_sequence[marker][4]) > 0:
        if len(set(marker_sequence[marker][5])) == 1:
            output_file.write(marker + '\t' + marker_sequence[marker][2] + '\t' + marker_sequence[marker][3] + '\t' + marker_sequence[marker][4][0] + '\t' + marker_sequence[marker][5][0] + '\t' + str(marker_sequence[marker][6][0]) + '\n')
        else:
            print('Multiple alleles', marker)

output_file.close()
