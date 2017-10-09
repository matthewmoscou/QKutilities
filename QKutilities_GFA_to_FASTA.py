#! /usr/bin/python

"""
%prog gfa fasta 
Reads GFA files and exports a FASTA file

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
"""

# modules
import optparse
from optparse import OptionParser

import string


# import arguments and options
usage = "usage: %prog gfa fasta"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()


# import sequence (GFA), export FASTA
gfa_file = open(args[0], 'r')
fasta_file = open(args[1], 'w')
	
for line in gfa_file.readlines():
	if len(line) > 0:
		if line[0] == 'S':
			sline = string.split(line, '\t')
			fasta_file.write('>' + sline[1] + '\n')
			fasta_file.write(sline[2] + '\n')

gfa_file.close()
fasta_file.close()
