#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")

import sys
import os
from Bio import SeqIO

if len(sys.argv) == 3:
	pass
else:
	print "============ genbank2faa ============"
	print "USAGE: genbank2faa.py <INPUT: genbank file > <OUTPUT: protein fasta file>"
	print "EG: genbank2faa.py leptospira.gbk leptospira.faa"
	sys.exit()

genbankHandle = open(sys.argv[1])
genbankParser = SeqIO.parse(genbankHandle,'genbank')
genbankCDSText = ""
cds_count = 0

for record in genbankParser:
	for feature in record.features:
		if feature.type == 'gene':
			locus_tag = 'NO_LOCUS_TAG'
			if 'locus_tag' in feature.qualifiers.keys():
				locus_tag = feature.qualifiers['locus_tag'][0]
		if feature.type == 'CDS':
			cds_count += 1
			product = 'MISSING_PRODUCT_NAME'
			if 'product' in feature.qualifiers.keys():
				product = feature.qualifiers['product'][0]
			genbankCDSText += ('>{0}|{7}|{1}|{2}:{3}:{4}| {5}\n'
				                   '{6}\n').format(str(cds_count),locus_tag,feature.location.start,
					       feature.location.end,feature.location.strand,product,
					       str(feature.location.extract(record.seq).translate()).replace('*',''),record.id)
fastaHandleOut = open(sys.argv[2],'w')
fastaHandleOut.write(genbankCDSText)
fastaHandleOut.close()
			
