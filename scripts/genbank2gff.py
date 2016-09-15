from Bio import SeqIO
from Bio.SeqFeature import *
import sys
import os

if len(sys.argv) == 0:
	pass

if os.path.isfile(sys.argv[1]):
	pass

else:
	sys.exit()

gff_file_content = '##gff-version 3\n'

count = 0

for record in SeqIO.parse(open(sys.argv[1]),'genbank'):
	count += 1
	for feature in record.features:
		if not feature:
			continue
		if feature.type == 'source':
			continue
		seqid = record.name
		partial = False
		start = feature.location.start
		end = feature.location.end
		strand = {1:'+',-1:'-'}[feature.location.strand]
		program = 'genix'
		phase = '.'
		if feature.type == 'CDS':
			if 'codon_start' in feature.qualifiers:
				phase = int(feature.qualifiers['codon_start'][0]) - 1 
		if isinstance(feature.location.start, AfterPosition):
			partial = True
		if isinstance(feature.location.end, AfterPosition):
			partial = True
		gff_file_content += seqid + '\t'
		gff_file_content += program + '\t'
		gff_file_content += feature.type + '\t'
		if strand == '+':
			gff_file_content += str(start.numerator+1) + '\t'
			gff_file_content += str(end.numerator) + '\t'
		else:
			gff_file_content += str(start.numerator+1) + '\t'
			gff_file_content += str(end.numerator) + '\t'
		gff_file_content += '.' + '\t'
		gff_file_content += strand + '\t'
		gff_file_content += str(phase) + '\t'
		atributos = []
		if 'locus_tag' in feature.qualifiers:
			atributos.append('ID=%s'%feature.qualifiers['locus_tag'][0])
		if 'product' in feature.qualifiers:
			atributos.append('Name=%s'%feature.qualifiers['product'][0].split(',')[0])
		gff_file_content += ';'.join(atributos)
		gff_file_content += '\n'
sys.stdout.write(gff_file_content)
