#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import sys
import os
from Bio import SeqIO

if len(sys.argv) == 3:
	pass
else:
	print "============ genbank2fna ============"
	print "USAGE: genbank2faa.py <INPUT: genbank file > <OUTPUT: non-coding DNA fasta file>"
	print "EG: genbank2faa.py leptospira.gbk leptospira.faa"
	sys.exit()

genbankHandle = open(sys.argv[1])
genbankParser = SeqIO.parse(genbankHandle,'genbank')
genbankCDSText = ""
cds_count = 0

for record in genbankParser:
	for feature in record.features:
		if not feature:
			continue
		if feature.type == 'source':
			continue
		if feature.type == 'CDS':
			continue
		if feature.type == 'gene':
			continue
		locus_tag = 'NO_LOCUS_TAG'
		if 'locus_tag' in feature.qualifiers.keys():
			locus_tag = feature.qualifiers['locus_tag'][0]
		cds_count += 1
		product = 'MISSING_PRODUCT_NAME'
		if 'product' in feature.qualifiers.keys():
			product = feature.qualifiers['product'][0]
		genbankCDSText += ('>{0}|{7}|{1}|{2}:{3}:{4}| {5} {8}\n'
				                   '{6}\n').format(str(cds_count),locus_tag,feature.location.start,
					       feature.location.end,feature.location.strand,product,
					       str(feature.location.extract(record.seq)),record.id,feature.type)
fastaHandleOut = open(sys.argv[2],'w')
fastaHandleOut.write(genbankCDSText)
fastaHandleOut.close()
			