#!/usr/bin/python
from __future__ import print_function, division
from collections import defaultdict
from optparse import OptionParser
import re,sys

usage = '%prog -[i]'
p = OptionParser(usage = usage)
p.add_option('-i','--input',dest='input',help='given an input table file')
(options,args) = p.parse_args()

if options.input:
	infile = open(options.input)
	outfile = open('percent.%s' % options.input, 'w')
	sampledict = defaultdict(float)
	for i in infile:
		line = re.split(r'\t', i.rstrip())
		if re.match(r'#', line[0]) and not re.match(r'#OTU|#taxa', line[0]):
			pass
		elif re.match(r'#OTU|#taxa', line[0]):
			sample = {}
			m = 1
			while m <= len(line) -2:
				sample[m] = line[m]
				sampledict[line[m]]
				m += 1
		else:
			m = 1
			while m <= len(line) -2:
				sampledict[sample[m]] += float(line[m])
				m += 1
	infile.close()

	infile = open(options.input)
	for i in infile:
		line = re.split(r'\t', i.rstrip())
		if re.match(r'#', line[0]) and not re.match(r'#OTU|#taxa', line[0]):
			pass
		elif re.match(r'#OTU|#taxa', line[0]):
			temp = '\t'
			m = 1
			while m <= len(line) -2:
				temp = temp + line[m] + '\t'
				m += 1
			print ('%s' % temp.rstrip(), file = outfile)
		else:
			m = 1
			total = 0
			temp = line[0] + '\t'
			while m <= len(line) -2:
				total += float(line[m])
				temp = temp + '%s' % format(float(line[m])/sampledict[sample[m]], 'e') + '\t'
				m += 1
			if total == 0:
				pass
			else:
				print ('%s' % temp.rstrip(), file = outfile)

	infile.close()
	outfile.close()

else:
	print ('error:must given a input file')
	sys.exit()
