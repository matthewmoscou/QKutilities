#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Determines the lengths of individual contigs/chromosomes from a FASTA file
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Reads FASTA and counts the number of nucleotide/amino acids in each sequence
	2. Exports sequence lengths
"""


## modules
import optparse
from optparse import OptionParser
import sets
import string


## OptionParser
# import arguments and options
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta", default="", help="FASTA file")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default="", help="Output for sequence lengths")
(options, args) = parser.parse_args()


## import fasta

fasta_file = open(options.fasta, 'r')

ID_seqlen = {}
ID_order = []

ID = ''
line = fasta_file.readline()

while line:
	if len(line) > 0:
		if line[0] == '>':
			ID = string.split(line)[0][1:]
			ID_seqlen[ID] = 0
			ID_order.append(ID)
		else:
			for sline in string.split(line):
				ID_seqlen[ID] += len(sline)

	line = fasta_file.readline()

fasta_file.close()

if len(ID_order) != len(sets.Set(ID_order)):
	print '\tError - FASTA file has redundant identifiers'


# export sequence lengths
length_out = open(options.output, 'w')

for ID in ID_order:
	length_out.write(ID + '\t' + str(ID_seqlen[ID]) + '\n')

length_out.close()
