#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Provides a detailed analysis of a protein sequence.
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Scans a protein sequence for the abundance of individual amino acids
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
parser.add_option("-w", "--window", action="store", type="int", dest="window", default=10, help="Window size in amino acids")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default="", help="Output for window-based coverage analysis")
(options, args) = parser.parse_args()

aminoacids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

## import fasta
if len(options.fasta) > 0:
	fasta_file = open(options.fasta, 'r')
	
	ID_sequence = {}
	
	ID = ''
	line = fasta_file.readline()
	
	while line:
		sline = string.split(line)

		if len(line) > 0:
			if line[0] == '>':
				ID = string.split(line)[0][1:]
				ID_sequence[ID] = ''
			else:
				ID_sequence[ID] += sline[0]
	
		line = fasta_file.readline()
	
	fasta_file.close()


## read and analyze coverage
output = open(options.output, 'w')

output.write('gene' + '\t' + 'index' + '\t' + 'aminoacid' + '\t' + 'frequency' + '\n')

for ID in ID_sequence.keys():
	for index in range(len(ID_sequence[ID]) - options.window):
		for aminoacid in aminoacids:
			output.write(ID + '\t' + str(index) + '\t' + aminoacid + '\t' + str(ID_sequence[ID][index:(index + options.window)].count(aminoacid)) + '\n')

output.close()

