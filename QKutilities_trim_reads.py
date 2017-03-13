#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
Trims the 5' and 3' ends of reads, also permits size selection (lower and upper bounds).
This performs the following:
	1. Reads FASTA and trims/size selects reads 

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
parser = OptionParser()
parser.add_option("-d", "--distance", action="store", type="int", dest="distance", default=0, help="Sequence to trim from start and end of reads")
parser.add_option("-l", "--lower", action="store", type="int", dest="lower", default=-1, help="Lower bound for size selection on reads")
parser.add_option("-u", "--upper", action="store", type="int", dest="upper", default=-1, help="Upper bound for size selection on reads")
parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", default="", help="Input FASTA file")
parser.add_option("-o", "--output", action="store", type="string", dest="outputfile", default="", help="Output FASTA file")
(options, args) = parser.parse_args()


# read transcripts
# after import, export trimmed and size selected sequence
reads = open(options.inputfile, 'r')
reads_trimmed = open(options.outputfile, 'w')

line = reads.readline()
sequence = ''

while line:
	if len(line) > 0:
		sline = string.split(line)

		if line[0] == '>':
			if len(sequence) > 0:
				trimmed_sequence = sequence[int(options.distance):len(sequence) - int(options.distance)]
				export = False

				# logic control for size selection
				if options.upper >= 0:
					if options.lower >= 0:
						if len(sequence) < options.upper:
							if len(sequence) > options.lower:
								export = True
							else:
								export = False
						else:
							export = False
					elif len(sequence) < options.upper:
						export = True
				elif options.lower >= 0:
					if len(sequence) > options.lower:
						export = True
				else:
					export = True

				if export:
					reads_trimmed.write('>' + ID + '\n')
					reads_trimmed.write(trimmed_sequence + '\n')

			ID = sline[0][1:]
			sequence = ''
		else:
			sequence += sline[0]

	line = reads.readline()

# export final sequence, if it meets the requirements
# logic control for size selection
trimmed_sequence = sequence[int(options.distance):len(sequence) - int(options.distance)]
export = False

if options.upper >= 0:
	if options.lower >= 0:
		if len(sequence) < options.upper:
			if len(sequence) > options.lower:
				export = True
			else:
				export = False
		else:
			export = False
	elif len(sequence) < options.upper:
		export = True
elif options.lower >= 0:
	if len(sequence) > options.lower:
		export = True
else:
	export = True

if export:
	reads_trimmed.write('>' + ID + '\n')
	reads_trimmed.write(trimmed_sequence + '\n')

reads.close()
reads_trimmed.close()
