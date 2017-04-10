#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] 
Assesses read coverage over large stretches of sequence space
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
	1. Reads genome coverage file and uses a window-based analysis for assessment
	2. (Optional) Use a masked FASTA file to skip masked regions in analysis
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
parser.add_option("-c", "--coverage", action="store", type="string", dest="coverage", default="", help="Genome coverage file from bedtools")
parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta", default="", help="FASTA file")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default="", help="Output for window-based coverage analysis")
parser.add_option("-w", "--window", action="store", type="int", dest="window", default=0, help="Window size in base pairs")
(options, args) = parser.parse_args()


## import fasta
if len(options.fasta) > 0:
	fasta_file = open(options.fasta, 'r')
	output_file = open(options.fasta + '.w' + str(options.window) + '.coverage', 'w')
	
	masked = {}
	
	ID = ''
	line = fasta_file.readline()
	
	while line:
		sline = string.split(line)

		if len(line) > 0:
			if line[0] == '>':
				ID = string.split(line)[0][1:]
				masked[ID] = []
				position = 0
			else:
				for index in range(len(sline[0])):
					while len(masked[ID]) <= ((position + index) / options.window):
						masked[ID].append(0)
				
					if sline[0][index] == 'N':
						masked[ID][(position + index) / options.window] += 1

				position += len(sline[0])
	
		line = fasta_file.readline()
	
	output_file.write('contig' + '\t' + 'window' + '\t' + 'masking' + '\n')

	for ID in masked.keys():
		for position in range(len(masked[ID])):
			output_file.write(ID + '\t' + str(position) + '\t' + str(masked[ID][position]) + '\n')

	fasta_file.close()
	output_file.close()


## read and analyze coverage

coverage_file = open(options.coverage, 'r')

coverage = {}
evaluted_sites = {}

line = coverage_file.readline()

while line:
	sline = string.split(line)

	if sline[0] not in coverage.keys():
		coverage[sline[0]] = []
		evaluted_sites[sline[0]] = []

	while len(coverage[sline[0]]) <= (int(sline[1]) / options.window):
		coverage[sline[0]].append(0)
		evaluted_sites[sline[0]].append(0)

	coverage[sline[0]][int(sline[1]) / options.window] += int(sline[2])
	evaluted_sites[sline[0]][int(sline[1]) / options.window] += 1

	line = coverage_file.readline()

coverage_file.close()

output = open(options.output, 'w')

output.write('contig' + '\t' + 'window' + '\t' + 'coverage' + '\n')

for ID in coverage.keys():
	for index in range(len(coverage[ID])):
		output.write(ID + '\t' + str(index) + '\t' + str(float(coverage[ID][index]) / float(evaluted_sites[ID][index])) + '\n')

output.close()
