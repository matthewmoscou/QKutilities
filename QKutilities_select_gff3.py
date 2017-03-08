#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Extracts a region within a pair of FASTA and GFF3 files
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
Extracts a region from sequence and annotations based on user-defined start and end positions
This performs the following:
	1. Reads FASTA and selects only the region of interest
	2. Reads GFF3 and exports all features within the region of interest

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
parser.add_option("-g", "--gff", action="store", type="string", dest="gff", default='', help="GFF3 file")
parser.add_option("-c", "--contig", action="store", type="string", dest="contig", default='', help="Chromosome or contig identifier")
parser.add_option("-s", "--start", action="store", type="int", dest="start", default=-1, help="Start position within chromosome/contig")
parser.add_option("-e", "--end", action="store", type="int", dest="end", default=-1, help="End position within chromosome/contig")
parser.add_option("-p", "--prefix", action="store", type="string", dest="prefix", default=-1, help="Prefix for FASTA and GFF3 files with region")
(options, args) = parser.parse_args()


## import fasta
# for memory efficiency, individual lines are evaluated
# the following logic tests must be evaluated at each line
#    before start
#    after end
#    midpoint
#    at start, not end
#    at end, not start
#    start and end

fasta_file = open(options.fasta, 'r')

ID = ''
line = fasta_file.readline()

while line:
	if len(line) > 0:
		if line[0] == '>':
			ID = string.split(line)[0][1:]
			position = 0
		else:
			if ID == options.contig:
				if position > options.start:
					if (position + len(string.split(line)[0])) < options.end:
						# midpoint
						sequence[1] += string.split(line)[0]
					elif position < options.end: 
						# at end, not start
						sequence[1] += string.split(line)[0][:(options.end - position)]
				elif (position + len(string.split(line)[0])) > options.start:
					if (position + len(string.split(line)[0])) < options.end:
						# at start, not end
						sequence = [ID, string.split(line)[0][(options.start - position):]]
					else:
						# start and end
						sequence = [ID, string.split(line)[0][(options.start - position):(options.end - position)]]

				position += len(string.split(line)[0])

	line = fasta_file.readline()

fasta_file.close()


# export fasta
fasta_out = open(options.prefix + '.fa', 'w')

fasta_out.write('>' + sequence[0] + '_' + str(options.start) + '_' + str(options.end) + '\n')
fasta_out.write(sequence[1] + '\n')

fasta_out.close()


## import/export gff3
# two ways to do this
# simple way, export everything within the interval
# complicated way, evaluate each gene model and comment if truncated within a region
#    perform analysis line by line, but only export a full gene model once it is complete

gff_file = open(options.gff, 'r')
gff_out = open(options.prefix + '.gff3', 'w')

line = gff_file.readline()

while line:
	if len(line) > 0:
		if line[0] == '#':
			gff_out.write(line)
		else:
			line = string.replace(line, '\n', '')
			sline = string.split(line, '\t')

			if sline[0] == options.contig:
				if int(sline[3]) >= options.start:
					if int(sline[4]) <= options.end:
						gff_out.write(sline[0] + '\t' + sline[1] + '\t' + sline[2] + '\t' + str(int(sline[3]) - options.start) + '\t' + str(int(sline[4]) - options.start) + '\t' + sline[5] + '\t' + sline[6] + '\t' + sline[7] + '\t' + sline[8] + '\n')

#				else:
#					print 'truncated gene model - start - ', sline[0], sline[3], sline[4], sline[8]
#			elif int(sline[4]) > options.end:
#				print 'truncated gene model - end - ', sline[0], sline[3], sline[4], sline[8]

	
	line = gff_file.readline()

gff_file.close()
gff_out.close()

