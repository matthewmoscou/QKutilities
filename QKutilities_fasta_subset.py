#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Extract a subset of sequences from a FASTA file
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Reads a list of sequence identifiers
	2. Reads the FASTA file and exports the subset of sequences

Improvements from existing script set (7 March 2017):

Future improvements to include:

"""

## modules
import optparse
from optparse import OptionParser
import string


## OptionParser
# import arguments and options
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta", default='', help="FASTA file")
parser.add_option("-s", "--subset", action="store", type="string", dest="subset", default='', help="List of sequence identifiers")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default=-1, help="FASTA file for the subset of sequences")
(options, args) = parser.parse_args()


## import list of subset
contigs = []

subset_file = open(options.subset, 'r')

for line in subset_file.readlines():
	sline = string.split(line)

	contigs.append(sline[0])

subset_file.close()


## read FASTA and export subset
fasta_file = open(options.fasta, 'r')
output_file = open(options.output, 'w')

line = fasta_file.readline()

while line:
	if len(line) > 0:
		if line[0] == '>':
			sline = string.split(line)

			if sline[0][1:] in contigs:
				file_index = contigs.index(sline[0][1:])
				output_file.write(sline[0] + '\n')
				contigs.remove(sline[0][1:])
			else:
				file_index = -1

		elif file_index >= 0:
			output_file.write(line)

	line = fasta_file.readline()

fasta_file.close()
output_file.close()

## error message if not all sequences were found
if len(contigs) > 0:
	print 'Error - The following sequences were not found in the FASTA file:'

	for contig in contigs:
		print '\t' + contig
