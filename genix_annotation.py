#!/usr/bin/env python

print '''
############################################################
#////////// Genix Annotation //////////////////////////////#
# An automated pipeline for bacterial genome annotation    #
# Version: 0.4                                             #
# Reference: KREMER et al, 2016                            #
############################################################
'''

import warnings
warnings.filterwarnings("ignore")

import re
import argparse
import sys
import os
import gc
import shutil
import subprocess
import sqlite3
import getpass
import glob
import shlex
import string

try:

	from Bio.Seq import *
	from Bio.Alphabet import *
	from Bio.SeqFeature import *
	from Bio.SeqRecord import *
	from Bio import SeqIO
	from Bio import SearchIO
	from Bio.Blast import NCBIXML
	from Bio.Alphabet import generic_dna

except:

	print 'ERROR: GENIX requires the Biopython package, which was not'
	print '       found in your machine. Please, install biopython'
	print '       before proceed to the annotation.\n'

	sys.exit()

#  -----------------------------------------------------------------------------
#  Genix:
#  -----------------------------------------------------------------------------
#  NOTE:
#  The this is an almost draft version of the pipeline and was not initially
#  designed to be published as a stand-alone application. Therefore, there are
#  many aspects that must be updated in future version, especially in terms of
#  code structure and standatization. As soon as possible, all code will be re-
#  written to take advantage of object-oriented programming, concurrence,
#  Biopython features, SQLite / MySQL and other this.
#  Current Features:
#  - gene prediction with Prodigal
#  - rRNA prediction with RNAmmer
#  - tRNA prediction with tRNAscan-SE
#  - ncRNA / regulatory / other stuff prediction with INFERNAL + Rfam
#
#  TODOS:
#  -----------------------------------------------------------------------------

# pega o local do script

PERL_PATH = '/usr/bin/perl'
PYTHON_PATH = '/usr/bin/python'
RNAMMER_PATH = '/home/cdtec/Frederico/programas/RNAmmer/rnammer'
script_location = os.path.dirname(os.path.abspath(sys.argv[0]))

#snippet para o devnull

devnull = open(os.devnull)

def format_sequence(sequence):

	sequence = sequence.upper()
	temp_seq = open('temp_seq_file.txt','w')
	temp_seq.write('>1\n'+sequence)
	temp_seq.close()
	global script_location
	os.system(('{1} {0}/scripts/genbank_sequence_format.pl temp_seq_file.txt '
	           ' temp_out_file.txt').format(script_location,PERL_PATH))

	resultado = open('temp_out_file.txt').read()
	os.remove('temp_out_file.txt')
	return resultado

def process_fasta(sequence):

	processed = ''

	for i in range(len(sequence)):

		if sequence[i] in 'ATCG':

			processed += sequence[i]

		else:

			processed += 'N'

	return processed

def extend_cds(prodigal_gene,scaffold,limit):
	global start_codons
	global stop_codons

	new_prodigal_gene = list(prodigal_gene)

	if prodigal_gene[4] == '+':
		old_start = prodigal_gene[2]
		new_start = old_start
		stop  = prodigal_gene[3]
		stop_codon = scaffold[2][prodigal_gene[3]-1:prodigal_gene[3]+2]

		while new_start >= limit:
			new_start -= 3

			new_start_codon = scaffold[2][new_start-1:new_start+2]
			old_start_codon = scaffold[2][old_start-1:old_start+2]

			if new_start_codon in stop_codons:
				break

			new_prodigal_gene_sequence = scaffold[2][new_start-1:prodigal_gene[3]]

			if new_start_codon in start_codons and not re.match("(.+?)?N+(.+?)?",new_prodigal_gene_sequence):
				new_prodigal_gene[2] = int(new_start)
				yield new_prodigal_gene

			else:

				pass

	if prodigal_gene[4] == '-':

		old_start = prodigal_gene[3]
		new_start = old_start
		stop = prodigal_gene[2]

		while new_start <= limit:

			new_start += 3

			new_start_codon = str(reverse_complement(Seq(scaffold[2][new_start-3:new_start])))
			old_start_codon = str(reverse_complement(Seq(scaffold[2][old_start-3:old_start])))

			if new_start_codon in stop_codons:

				break

			new_prodigal_gene_sequence = str(reverse_complement(Seq(scaffold[2][stop-1:new_start])))

			if new_start_codon in start_codons and not re.match("[ACTG]+?N+?[ACTG]+?",
				                                                 new_prodigal_gene_sequence):

				new_prodigal_gene[3] = int(new_start)

				yield new_prodigal_gene

			else:

				pass

def reduce_cds(prodigal_gene,scaffold,limit=0):
	global database_connection
	global database_cursor
	limit=prodigal_gene[3] if prodigal_gene[4] == '+' else prodigal_gene[2]
	global start_codons
	global stop_codons

	new_prodigal_gene = list(prodigal_gene)

	if prodigal_gene[4] == '+':

		old_start = prodigal_gene[2]
		new_start = old_start
		stop  = prodigal_gene[3]
		stop_codon = scaffold[2][prodigal_gene[3]-1:prodigal_gene[3]+2]

		while new_start <= limit-3:

			new_start += 3

			new_start_codon = scaffold[2][new_start-1:new_start+2]
			old_start_codon = scaffold[2][old_start-1:old_start+2]
			new_prodigal_gene_sequence = scaffold[2][new_start-1:prodigal_gene[3]]

			if new_start_codon in start_codons and not re.match("[ACTG]+?N+?[ACTG]+?",
				                                                 new_prodigal_gene_sequence):
				new_prodigal_gene[2] = int(new_start)

				yield new_prodigal_gene

			else:

				pass

	if prodigal_gene[4] == '-':

		old_start = prodigal_gene[3]
		new_start = old_start
		stop = prodigal_gene[2]

		while new_start >= limit+3:

			new_start -= 3

			new_start_codon = str(reverse_complement(Seq(scaffold[2][new_start-3:new_start])))
			old_start_codon = str(reverse_complement(Seq(scaffold[2][old_start-3:old_start])))

			if new_start_codon in stop_codons:

				break

			new_prodigal_gene_sequence = str(reverse_complement(Seq(scaffold[2][stop-1:new_start])))

			if new_start_codon in start_codons and not re.match("[ACTG]+?N+?[ACTG]+?",
				                                                 new_prodigal_gene_sequence):

				new_prodigal_gene[3] = int(new_start)

				yield new_prodigal_gene

			else:

				pass

#PARSE ARGUMENTS

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", help='Contigs in FASTA format',
	                type=str,required=True)

parser.add_argument("-o","--output_dir",help='The output dir',
	                default=os.getcwd())

parser.add_argument("-t","--threads", help='Number of threads',
	                default=1)

parser.add_argument("-dbp","--protein_database",
	                help='BLAST protein database',
	                required=True)

parser.add_argument("-antie","--antifam_evalue",
	                help=('Minimal e-value to a ORF be considered spurious '
	                	  'in the AntiFam analysis'),
	                default=0.00001)
parser.add_argument("-blaste","--blast_evalue",
	                default=0.0001)

parser.add_argument("-blastrnae","--blast_rna_evalue",
	                default=0.0001)

parser.add_argument("-nfernale","--infernal_evalue",default=0.0001)

parser.add_argument("-sub","--submission_template",
	                help='Template file for NCBI submission')

parser.add_argument("-mg","--min_gap",
	                help='Minimum length of a gap to annotate',
	                default=1)

parser.add_argument("-gc","--genetic_code", help='genetic code table',
	                default='11',choices=['11','4'])

parser.add_argument("-sl","--scaffold_label", help='Tag name for the scaffolds',
	                default='SCAFFOLD')

parser.add_argument("-lt","--locus_tag", help='Locus_tag prefix',default='LC')

parser.add_argument("-nst","--institution", help='Name of the research center',
	                type=str,default='genix')

# "search_program" it will future allow other tools, including USEARCH / UBLAST, BLAT and WU-BLAST.
# USEARCH was added in a previous version, but i removed because the implementation was very "spaghetti"

parser.add_argument("-sp","--search_program",default="BLAST",choices=['BLAST'])

parser.add_argument("-source","--source_information",
	                help='String containing the source information for tbl2asn',
	                type=str,default='')

parser.add_argument("-gap","--gap_evidence",
	                help='Information used for scaffolding',
	                choices=['paired-ends', 'align-genus', 'align-xgenus',
	                         'align-trnscpt', 'within-clone', 'clone-contig',
	                         'map', 'strobe'],
	                default='paired-ends')

parser.add_argument("-ugl","--unknown_gap_length",
	                help=('Use if the length of the gaps '
	                	  'inside de scaffolds is unknow'))

arguments = parser.parse_args()


#start_codon tables and stop_codon_tables

if arguments.genetic_code == '11':

	start_codons = ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
	stop_codons = ['TAA', 'TAG', 'TGA']

if arguments.genetic_code == '4':

	start_codons = ['TTA','TTG','CTG','ATT','ATC','ATA','ATG','GTG']
	stop_codons = ['TAA','TAG']

#CHECK DIR

if os.path.isdir(arguments.output_dir):

	pass

else:

	try:

		os.mkdir(arguments.output_dir)

	except:

		sys.stderr.write('ERROR: could not create the "{0}" directory.'.format(
			              arguments.output_dir))
		sys.stderr.write('\n')
		sys.exit()

#CHECK FILES

#Checking sequence

if os.path.isfile(arguments.input):

	seq_handle = open(arguments.input)
	sequences = []

	for record in SeqIO.parse(seq_handle,'fasta'):

		sequences.append(record)

	if len(sequences) == 0:

		sys.stderr.write('ERROR: "{0}" has no sequence.'.format(
			              arguments.input))
		sys.stderr.write('\n')
		sys.exit()

else:

	sys.stderr.write('ERROR: "%s" was not found.\n'%arguments.input)
	sys.exit()

#  -----------------------------------------------------------------------------
#  Prepare the database
#  -----------------------------------------------------------------------------
#  -----------------------------------------------------------------------------

database_connection = sqlite3.connect(':memory:')
database_cursor = database_connection.cursor()

database_cursor.execute('''
	CREATE TABLE sequences (id integer PRIMARY KEY AUTOINCREMENT,
	header text, sequence text)
    ''')

database_cursor.execute('''
	CREATE TABLE gaps (id integer PRIMARY KEY AUTOINCREMENT,
	sequence_id integer, start integer, end integer)
	''')

database_cursor.execute('''
	CREATE TABLE prodigal (id integer PRIMARY KEY AUTOINCREMENT,
	sequence_id integer, start integer, end integer, strand text)
	''')

database_cursor.execute('''
	CREATE TABLE blast_hits (id integer PRIMARY KEY AUTOINCREMENT,
	prodigal_id integer, hit_name text, query_start integer,
	query_end integer, query_len integer, hit_start integer,
	hit_end integer, hit_len integer, evalue float,evidence string,
	coverage float,uniprot text)
	''')

database_cursor.execute('''
	CREATE TABLE antifam (id integer PRIMARY KEY AUTOINCREMENT,
	prodigal_id integer, correct text)
	''')

database_cursor.execute('''
	CREATE TABLE aragorn (id integer PRIMARY KEY AUTOINCREMENT,
	sequence_id integer, start integer, end integer, strand string,
	tag_peptide text)
	''')

database_cursor.execute('''
	CREATE TABLE trnascan (id integer PRIMARY KEY AUTOINCREMENT,
	sequence_id integer, start integer, end integer, strand string,
	anticodon string, aminoacid string, score float)
	''')

database_cursor.execute('''
	CREATE TABLE rnammer (id integer PRIMARY KEY AUTOINCREMENT,
	sequence_id integer, start integer, end integer, strand string,
	definition strings)
	''')

database_cursor.execute('''
	CREATE TABLE rfam (id integer PRIMARY KEY AUTOINCREMENT, rfam_id text,
	rfam_symbol text, rfam_product text, rfam_class text)
	''')

database_cursor.execute('''
	CREATE TABLE infernal (id integer PRIMARY KEY AUTOINCREMENT,
	sequence_id integer, start integer, end integer, strand string,
	definition string)
	''')

database_cursor.execute('''
	CREATE TABLE raw_features (id integer PRIMARY KEY AUTOINCREMENT,
	sequence_id integer, start integer, end integer, strand text,
	type text, data_id integer, spurious text)
	''')

database_cursor.execute('''
	CREATE TABLE final_features (id integer PRIMARY KEY AUTOINCREMENT,
	raw_feature_id integer)
	''')

database_cursor.execute('''
	CREATE TABLE putative_frameshifted (id integer PRIMARY KEY AUTOINCREMENT,
	raw_feature_id integer)
	''')

database_connection.commit()

#  -----------------------------------------------------------------------------
#  Index the Rfam database
#  -----------------------------------------------------------------------------
#  -----------------------------------------------------------------------------

sys.stdout.write('[Indexing]: Indexing The Rfam database.\n')

rfam_family_file = open('{0}/databases/rfam/family.txt'.format(
	                    script_location)).read().split('\n')

rfam_family_file.pop(-1)

rfam_dados = [entry.split('\t') for entry in rfam_family_file]

for dados in rfam_dados:
	database_cursor.execute("INSERT INTO rfam VALUES (NULL,?,?,?,?)",
                            (dados[0], dados[1], dados[3], dados[18]))

database_connection.commit()

del rfam_family_file
del rfam_dados

#  -----------------------------------------------------------------------------
#  Genix:
#  -----------------------------------------------------------------------------
#  -----------------------------------------------------------------------------

sys.stdout.write('[Indexing]: Indexing sequences.\n')

for sequence in sequences:

	sequence_raw = str(sequence.seq)
	for whitespace in string.whitespace:

		sequence_raw = sequence_raw.replace(whitespace,'')
		sequence_raw = sequence_raw.upper()
		processed_sequence = process_fasta(sequence_raw)

	database_cursor.execute("INSERT INTO sequences VALUES (NULL,?,?)",
		(sequence.name,processed_sequence))

del sequences
gc.collect()
database_connection.commit()

#MAP GAPS
sys.stdout.write('[Indexing]: Indexing gaps.\n')

for sequence_entry in database_cursor.execute('SELECT * FROM sequences'):

	sequence = sequence_entry[2]
	sequence_id = sequence_entry[0]
	pattern_iter = re.finditer('[ATCG]N+[ATCG]',sequence)

	for pattern in pattern_iter:

		pattern_start = pattern.start()+2
		pattern_end = pattern.end()-1
		gap_len = len(sequence[pattern_start:pattern_end])

		if gap_len >= arguments.min_gap:

			database_cursor.execute("INSERT INTO gaps VALUES (NULL,?,?,?)",
				                     (sequence_id,pattern_start,pattern_end))
database_connection.commit()

#  -----------------------------------------------------------------------------
#  RUN Prodigal
#  -----------------------------------------------------------------------------
#
#  -----------------------------------------------------------------------------

sys.stdout.write('[Prodigal]: Predicting CDSs.\n')

processed_fasta_string = ''
processed_fasta_handle = open('{0}/processed_sequences.fasta'.format(
	                          arguments.output_dir),'w')

for sequence_entry in database_cursor.execute('SELECT * FROM sequences'):

	processed_fasta_string += '>%s\n%s\n'%(str(sequence_entry[0]),
		                       sequence_entry[2])

processed_fasta_handle.write(processed_fasta_string)
processed_fasta_handle.close()

prodigal_return_code = subprocess.call(('{0}/bin/prodigal/./prodigal '
	                                    '-i {1}/processed_sequences.fasta '
	                                    '-o {1}/prodigal_genes.sco -f sco '
	                                    '-m -g {2}').format(script_location,
	                                    	arguments.output_dir,
	                                    	arguments.genetic_code),
                                        shell=True,
                                        stderr=devnull,
                                        stdout=devnull)

if prodigal_return_code != 0:

	sys.stderr.write('ERROR: Error while running Prodigal.\n')
	sys.exit()

prodigal_output = open('%s/prodigal_genes.sco'%arguments.output_dir)

for line in prodigal_output:

	if re.match('^#.+seqhdr=.+$',line):

		sequence_id = int(line.split('seqhdr="')[1].split('"')[0])

	if re.match('^>[0-9]+_[0-9]+_[0-9]+_[\+-]\n$',line):

		prodigal_gene = line.replace('\n','').split('_')
		prodigal_gene_start = int(prodigal_gene[1])
		prodigal_gene_end = int(prodigal_gene[2])
		prodigal_gene_strand = prodigal_gene[3]
		database_cursor.execute("INSERT INTO prodigal VALUES (NULL,?,?,?,?)",
			                     (sequence_id,prodigal_gene_start,
			                      prodigal_gene_end,prodigal_gene_strand))
database_connection.commit()

#os.remove('%s/prodigal_genes.sco'%arguments.output_dir)

#  -----------------------------------------------------------------------------
#  RUN tRNAscan-SE
#  -----------------------------------------------------------------------------
#
#  -----------------------------------------------------------------------------

sys.stdout.write('[tRNAscan-SE]: Predicting tRNAs.\n')

abs_output_dir = os.path.abspath(arguments.output_dir)

current_dir = os.getcwd()
os.chdir('%s/bin/tRNAscan-SE'%script_location)
os.environ["PATH"] = os.environ["PATH"]+':'+os.getcwd()

if "PERL5LIB" in os.environ.keys():
	os.environ["PERL5LIB"] = os.environ["PERL5LIB"]+':/usr/%s/bin'%getpass.getuser()
else:
	os.environ["PERL5LIB"] = '/usr/%s/bin'%getpass.getuser()

os.environ["MANPATH"] = '/usr/%s/man'%getpass.getuser()

if os.path.isfile('%s/trnascan.out'%abs_output_dir):
	os.remove('%s/trnascan.out'%abs_output_dir)

trnascan_return_code = subprocess.call(('./tRNAscan-SE -Q -B '
	                                   '-o {0}/trnascan.out '
	                                   '{0}/processed_sequences.fasta').format(
	                                   abs_output_dir,abs_output_dir),
	                                   shell=True,
	                                   stderr=devnull,
	                                   stdout=devnull)
if trnascan_return_code != 0:

	sys.stderr.write('ERROR: Error while running tRNAscan-SE.\n')
	sys.exit()

os.chdir(current_dir)

trnascan_output = open('%s/trnascan.out'%arguments.output_dir)

for line in trnascan_output:

	data = line.split("\t")

	if data[0] == 'Sequence':

		continue

	elif data[0] == 'Name    ':

		continue

	elif data[0] == '--------':

		continue

	else:

		sequence_id = int(data[0])

		if int(data[3]) > int(data[2]):

			trnascan_trna_start = int(data[2])
			trnascan_trna_end = int(data[3])
			trnascan_trna_strand = '+'

		else:

			trnascan_trna_start = int(data[3])
			trnascan_trna_end = int(data[2])
			trnascan_trna_strand = '-'

		trnascan_trna_aminoacid = data[4]
		trnascan_trna_anticodon = data[5]
		trnascan_trna_coven_sco = data[8]

		if trnascan_trna_aminoacid != 'Pseudo':

			database_cursor.execute(
				"""INSERT INTO trnascan VALUES (NULL,?,?,?,?,?,?,?)""",
				   (sequence_id,trnascan_trna_start,trnascan_trna_end,
				    trnascan_trna_strand,trnascan_trna_anticodon,
				    trnascan_trna_aminoacid,trnascan_trna_coven_sco))

database_connection.commit()
os.remove('%s/trnascan.out'%arguments.output_dir)

#  -----------------------------------------------------------------------------
#  RUN RNAmmer
#  -----------------------------------------------------------------------------
#  TODO: A better parser for the RNAmmer XML. (i hate XML)#
#  -----------------------------------------------------------------------------

if os.path.isfile(RNAMMER_PATH):

	sys.stdout.write('[RNAmmer]: Predicting rRNAs.\n')
	rnammer_return_code = subprocess.call(
		                  ('{2} {0} -S bac -m lsu,ssu,tsu -xml '
		                   '{1}/rnammer.xml < {1}/processed_sequences.fasta'
	                      ).format(RNAMMER_PATH,arguments.output_dir,PERL_PATH),
	                      shell=True,
	                      stderr=devnull,
	                      stdout=devnull)

	if rnammer_return_code != 0:

		sys.stderr.write('ERROR: Error while running RNAmmer.\n')
		sys.exit()

	file_patterns = ['*.neg.fsa.hmmsearchresult','*.tsu.xml','*.ssu.xml',
	                 '*.ssu.cf','*.lsu.xml','*.pos.fsa','*temp.*']

	for file_pattern in file_patterns:

		for file in glob.glob('{0}/*{1}'.format(os.getcwd(),file_pattern)):

			os.remove(file)

	rnammer_output = open('%s/rnammer.xml'%arguments.output_dir)

	for line in rnammer_output:

		if re.match('^.+?<entry>$',line):

			rnammer_entry = {}

		if re.match('^.+?<mol>.+</mol>$',line):

			rnammer_entry['mol_def'] = line.split('<mol>')[1].split(
				                                  '</mol>')[0]

		if re.match('^.+?<start>.+</start>$',line):

			rnammer_entry['start'] = int(line.split('<start>')[1].split(
				                                    '</start>')[0])

		if re.match('^.+?<stop>.+</stop>$',line):

			rnammer_entry['end'] = int(line.split('<stop>')[1].split(
				                                  '</stop>')[0])

		if re.match('^.+?<direction>.+</direction>$',line):

			rnammer_entry['direction'] = line.split('<direction>')[1].split(
				                                    '</direction>')[0]

		if re.match('^.+?<sequenceEntry>.+</sequenceEntry>$',line):

			rnammer_entry['sequenceEntry'] = int(line.split('<sequenceEntry>'
				                                 )[1].split('</sequenceEntry>'
				                                 )[0])

		if re.match('^.+?</entry>$',line):

			database_cursor.execute(
				"INSERT INTO rnammer VALUES (NULL,?,?,?,?,?)",
				(rnammer_entry['sequenceEntry'],rnammer_entry['start'],
				rnammer_entry['end'],rnammer_entry['direction'],
				rnammer_entry['mol_def']))

database_connection.commit()

#  -----------------------------------------------------------------------------
#  RUN Aragorn
#  -----------------------------------------------------------------------------
#  TODO: A better parser for the RNAmmer XML. (i hate XML)#
#  -----------------------------------------------------------------------------

sys.stdout.write('[ARAGORN]: Predicting tmRNAs.\n')

aragorn_return_code = subprocess.call(('{0}/bin/aragorn/aragorn -m -gcbact '
	                                   '-w -o {1}/aragorn_output.txt '
	                                   '{1}/processed_sequences.fasta').format(
	                                    script_location,arguments.output_dir),
	                                   shell=True,
	                                   stderr=devnull,
	                                   stdout=devnull)

if aragorn_return_code != 0:

	sys.stderr.write('ERROR: Error while running Aragorn.\n')
	sys.exit()

saida_aragorn = open('%s/aragorn_output.txt'%(arguments.output_dir)).read()

for line in saida_aragorn.split('\n'):

	if '>' in line:

		header = line.split(">")[1]

	elif 'genes found' in line:

		pass

	elif 'tmRNA' in line:

		tm_rna_data = line.split(" ")[-1].split("\t")
		tm_rna_location = tm_rna_data[0]
		tm_rna_peptide = tm_rna_data[2]

		if 'c[' in tm_rna_location:

			strand = '-'

		else:

			strand = '+'

		tm_rna_start = tm_rna_location.split('[')[1].split(',')[0]
		tm_rna_end = tm_rna_location.split(',')[1].split(']')[0]
		database_cursor.execute("INSERT INTO aragorn VALUES (NULL,?,?,?,?,?)",
			                    (header,tm_rna_start,tm_rna_end,strand,
								 tm_rna_peptide))

database_connection.commit()


#  -----------------------------------------------------------------------------
#  RUN INFERNAL
#  -----------------------------------------------------------------------------
#  These families of RNAs are already predicted by other tools, which are more
#  purpose-specific, so they are not included when found by INFERNAL.
#  -----------------------------------------------------------------------------

excluded_families = ['RF00177','RF02541','RF00005','RF00001','RF01959',
                     'RF01960','RF01959','RF01960','RF01118','RF00023']

sys.stdout.write('[INFERNAL]: Predicting ncRNAs.\n')
database_cursor.execute('SELECT * FROM sequences')
sequence_entries = database_cursor.fetchall()

SQL_inserts_blastn = []

for sequence_entry in sequence_entries:

	blast_hits_accs = []
	sequence = sequence_entry[2]
	sequence_id = sequence_entry[0]
	blast_fasta_handle = open('%s/temp.fasta'%arguments.output_dir,'w')
	blast_fasta_handle.write('>%s\n%s\n'%(sequence_id,sequence))
	blast_fasta_handle.close()

	subprocess.call(('/usr/bin/blastn -word_size 10 -query {0}/temp.fasta '
			            '-db {1}/databases/rfam/rfam -out {0}/temp.xml '
                                    '-num_alignments 100000 '
			            '-outfmt 5 -num_threads {2}').format(
			            	arguments.output_dir,script_location,
			            	arguments.threads),
			            shell=True,
			            stdout=devnull,
			            stderr=devnull)
	try:

		rfam_blast_parser = NCBIXML.parse(open('{0}/temp.xml'.format(
			                arguments.output_dir)))

	except:
		rfam_blast_parser = []
		continue

	for record in rfam_blast_parser:

		for alignment in record.alignments:

			adicionar_record = False

			for hsp in alignment.hsps:

				if hsp.expect <= arguments.blast_rna_evalue:

					adicionar_record = True

				break


			if adicionar_record:

				if re.match('.+RF[0-9]+;.+',alignment.title):

					for pattern_found in re.finditer('RF[0-9]+;',

						alignment.title):
						blast_hits_accs.append(
							alignment.title[
							pattern_found.start():pattern_found.end()-1
							])

	blast_hits_accs = set(blast_hits_accs)

	for accs in blast_hits_accs:

		subprocess.call(('{0}/bin/infernal/src/cmsearch '
			         '--tblout {1}/temp_infernal.txt --rfam --cpu {3} '
			         '{0}/databases/rfam/CMs/{2}.cm {1}/temp.fasta').format(
			            	script_location,arguments.output_dir,accs,arguments.threads),
		               	shell=True,
		               	stdout=devnull,
		               	stderr=devnull)

		infernal_out_handle = open('%s/temp_infernal.txt'%arguments.output_dir
			                      ).read().replace("'",'').split('\n')

		infernal_out_lines = [shlex.split(line) for line in infernal_out_handle\
		                      if not re.match('[(^#(.+$)(^$)]',line)]

		infernal_out_lines.pop(-1)

		for infernal_line in infernal_out_lines:

			if infernal_line[3] not in excluded_families:

				if float(infernal_line[15]) <= arguments.infernal_evalue:

					if infernal_line[9] == '+':

						database_cursor.execute("""INSERT INTO infernal VALUES
							                    (NULL,?,?,?,?,?)""",
							                    (sequence_entry[0],
							                    infernal_line[7],
							                    infernal_line[8],
							                    infernal_line[9],
							                    infernal_line[3]))

					else:

						database_cursor.execute("""INSERT INTO infernal VALUES
							                     (NULL,?,?,?,?,?)""",
							                     (sequence_entry[0],
							                      infernal_line[8],
							                      infernal_line[7],
							                      infernal_line[9],
							                      infernal_line[3]))

database_connection.commit()
#  -----------------------------------------------------------------------------
#  RUN BLAST vs. Uniprot and HMMER vs. Antifam
#  -----------------------------------------------------------------------------
#  These families of RNAs are already predicted by other tools, which are more
#  purpose-specific, so they are not included when found by INFERNAL.
#  -----------------------------------------------------------------------------

sys.stdout.write('[Uniprot-BLAST]: Annotating proteins.\n')
database_cursor.execute('SELECT * FROM sequences')
sequence_entries = database_cursor.fetchall()

SQL_inserts_antifam = []
SQL_inserts_blastp = []

for sequence_entry in sequence_entries:
	sequence = sequence_entry[2]
	sequence_id = sequence_entry[0]
	database_cursor.execute('SELECT * FROM prodigal \
		                     WHERE sequence_id = {0}'.format(sequence_id))

	prodigal_genes = database_cursor.fetchall()

	for prodigal_gene in prodigal_genes:

		gene_sequence = Seq(sequence[prodigal_gene[2]-1:prodigal_gene[3]],
			                IUPAC.ambiguous_dna)

		if prodigal_gene[4] == '-':

			gene_sequence = gene_sequence.reverse_complement()

		protein_sequence = str(gene_sequence.translate(
			                   table=11)).replace("*",'')
		blast_fasta_handle = open('%s/temp.fasta'%arguments.output_dir,'w')
		blast_fasta_handle.write('>{0}_{1}\n{2}\n'.format(sequence_id,
			                      prodigal_gene[0],protein_sequence))
		blast_fasta_handle.close()

		#RUN AntiFam

		prodigal_correct = 'True'
		subprocess.call(('{0}/bin/hmmer3/src/./hmmscan --cpu {1} '
			           '-o {2}/temp_hmmer.txt {0}/databases/antifam/AntiFam.hmm '
			           '{2}/temp.fasta').format(script_location,
			           arguments.threads,arguments.output_dir),
			           shell=True,
			           stdout=devnull,
			           stderr=devnull)

		antifam_output = SearchIO.parse('%s/temp_hmmer.txt'%arguments.output_dir,
			                            'hmmer3-text')
		for resultado in antifam_output:

			for hsp in resultado.hsps:

				if hsp.evalue <=  arguments.antifam_evalue:

					prodigal_correct = 'False'
					break

		#database_cursor.execute("INSERT INTO antifam VALUES (NULL,?,?)",
		#	                    (prodigal_gene[0],prodigal_correct))
		SQL_inserts_antifam.append((prodigal_gene[0],prodigal_correct))
		#database_connection.commit()

		#RUN database search

		hit_name_name = 'No Hit'
		hit_uniprot_id = ''
		blast_hit_q_start = 0
		blast_hit_q_end = 0
		blast_hit_q_len = 0
		blast_hit_s_start = 0
		blast_hit_s_end = 0
		blast_hit_s_len = 0
		blast_evalue = 1
		hit_evidence = 'No'

		if arguments.search_program == 'BLAST':

			subprocess.call('{0}/bin/ncbi-blast/blastp -word_size 7 -query {1}/temp.fasta '
				            '-db {2} -out {1}/temp.xml -outfmt 5 '
				            '-num_threads {3}'.format(
					            script_location,arguments.output_dir,
					            arguments.protein_database,arguments.threads),
				            shell=True)
			hit_name_name = 'No Hit'
			hit_uniprot_id = ''
			blast_hit_q_start = 0
			blast_hit_q_end = 0
			blast_hit_q_len = 0
			blast_hit_s_start = 0
			blast_hit_s_end = 0
			blast_hit_s_len = 0
			blast_evalue = 1
			has_hit = False

			try:
				parser_blast_out = NCBIXML.parse(open('{0}/temp.xml'.format(
                                                 arguments.output_dir)))
				for record in parser_blast_out:

					for alignment in record.alignments:

						alignment.title = alignment.title.replace("#",' ')
						blast_hit_name = alignment.title


						for hsp in alignment.hsps:

							if hsp.expect <= arguments.blast_evalue:

								has_hit = True
								blast_hit_q_start = hsp.query_start
								blast_hit_q_end = hsp.query_end
								blast_hit_q_len = record.query_length
								blast_hit_s_start = hsp.sbjct_start
								blast_hit_s_end = hsp.sbjct_end
								blast_hit_s_len = alignment.length
								blast_evalue = hsp.expect
								hit_name_data = blast_hit_name.split('|')
								hit_uniprot_id = ''
								hit_uniprot_id = hit_name_data[3]
								hit_name_name = hit_name_data[-1].split(' ',1)[1]
								hit_name_name = re.split('[A-Z]+=',hit_name_name)[0].replace('[','(').replace(']',')')
								hit_evidence = hit_name_data[-1].split(' ',1)[1].split("PE=")[1].split(' ')[0]

								#database_cursor.execute((
								#	"INSERT INTO blast_hits VALUES "
								#	"(NULL,?,?,?,?,?,?,?,?,?,?,?,?)"),
								#	(prodigal_gene[0],hit_name_name,
								#		blast_hit_q_start,blast_hit_q_end,
				                #        blast_hit_q_len,blast_hit_s_start,
				                #        blast_hit_s_end,blast_hit_s_len,
						        #        blast_evalue,hit_evidence,1,hit_uniprot_id))
						        SQL_inserts_blastp.append((prodigal_gene[0],hit_name_name,
										                   blast_hit_q_start,blast_hit_q_end,
				                                           blast_hit_q_len,blast_hit_s_start,
				                                           blast_hit_s_end,blast_hit_s_len,
						                                   blast_evalue,hit_evidence,1,hit_uniprot_id))
							break

			except:

				print "ERROR WHILE PARSING BLAST OUTPUT"
				continue

			if not has_hit:

						#database_cursor.execute(("INSERT INTO blast_hits VALUES "
						#	                    "(NULL,?,?,?,?,?,?,?,?,?,?,?,?)"),
						#						(prodigal_gene[0],hit_name_name,
						#						blast_hit_q_start,blast_hit_q_end,
				        #                         blast_hit_q_len,blast_hit_s_start,
				        #                         blast_hit_s_end,blast_hit_s_len,
						#                         blast_evalue,hit_evidence,0,hit_uniprot_id))
						SQL_inserts_blastp.append((prodigal_gene[0],hit_name_name,
												blast_hit_q_start,blast_hit_q_end,
				                                 blast_hit_q_len,blast_hit_s_start,
				                                 blast_hit_s_end,blast_hit_s_len,
						                         blast_evalue,hit_evidence,0,hit_uniprot_id))

        try:

        	os.remove('%s/temp.fasta'%arguments.output_dir)
        	os.remove('%s/temp.xml'%arguments.output_dir)
        	os.remove('%s/temp_hmmer.txt'%arguments.output_dir)

        except:

        	pass
database_cursor.executemany("INSERT INTO antifam VALUES (NULL,?,?)",
		                    SQL_inserts_antifam)

database_cursor.executemany(("INSERT INTO blast_hits VALUES "
	                         "(NULL,?,?,?,?,?,?,?,?,?,?,?,?)"),
                             SQL_inserts_blastp)

del SQL_inserts_antifam
del SQL_inserts_blastp
database_connection.commit()
gc.collect()

#  -----------------------------------------------------------------------------
#  Prepare table of "raw features"
#  -----------------------------------------------------------------------------
#  These families of RNAs are already predicted by other tools, which are more
#  purpose-specific, so they are not included when found by INFERNAL.
#  -----------------------------------------------------------------------------

sys.stdout.write('[Final processing]: Processing annotation.\n')

database_cursor.execute('SELECT * FROM gaps')
entries = database_cursor.fetchall()

for gap in entries:

	database_cursor.execute("INSERT INTO raw_features VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",(
		                    gap[1],gap[2],gap[3],'+','gap',gap[0],'True'))

database_cursor.execute('SELECT * FROM prodigal')
entries = database_cursor.fetchall()

for orf in entries:

	database_cursor.execute("SELECT * FROM antifam WHERE prodigal_id = %s"%(orf[0]))
	resultado_antifam = database_cursor.fetchone()
	database_cursor.execute("SELECT * FROM blast_hits WHERE prodigal_id = %s"%(orf[0]))
	resultado_blast = database_cursor.fetchall()

	if len(resultado_blast) > 0:

		if resultado_blast[0][2] != 'No Hit':

			database_cursor.execute("INSERT INTO raw_features VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",
					               (orf[1],orf[2],orf[3],orf[-1],'orf',orf[0],resultado_antifam[2]))
			database_connection.commit()

database_cursor.execute('SELECT * FROM trnascan')
entries = database_cursor.fetchall()

for trna in entries:

	database_cursor.execute("INSERT INTO raw_features VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",
		                    (trna[1],trna[2],trna[3],trna[4],'trna',trna[0],'True'))
        database_connection.commit()

database_cursor.execute('SELECT * FROM rnammer')
entries = database_cursor.fetchall()

for rrna in entries:

	database_cursor.execute("INSERT INTO raw_features VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",
		                    (rrna[1],rrna[2],rrna[3],rrna[4],'rrna',rrna[0],'True'))
        database_connection.commit()

database_cursor.execute('SELECT * FROM infernal')
entries = database_cursor.fetchall()

for ncrna in entries:

	database_cursor.execute("INSERT INTO raw_features VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",
		                    (ncrna[1],ncrna[2],ncrna[3],ncrna[4],'ncrna',ncrna[0],'True'))
	database_connection.commit()

database_cursor.execute('SELECT * FROM aragorn')
entries = database_cursor.fetchall()

for tmrna in entries:

	database_cursor.execute("INSERT INTO raw_features VALUES (NULL, ?, ?, ?, ?, ?, ?, ?)",
		                     (tmrna[1],tmrna[2],tmrna[3],tmrna[4],'tmrna',tmrna[0],'True'))
	database_connection.commit()

#  -----------------------------------------------------------------------------
#  Correct start codons based on BLAST results
#  -----------------------------------------------------------------------------
#  These families of RNAs are already predicted by other tools, which are more
#  purpose-specific, so they are not included when found by INFERNAL.
#  -----------------------------------------------------------------------------

sys.stdout.write('[Final processing]: Correcting start positions.\n')

database_cursor.execute('SELECT * FROM sequences')
sequencias = database_cursor.fetchall()

# Iterate all sequences

for sequencia in sequencias:
	database_cursor.execute('SELECT * FROM prodigal WHERE sequence_id = %s'%sequencia[0])
	prodigal_genes = database_cursor.fetchall()

	# iterate all prodigal genes predicted in the given sequence

	for prodigal_gene_index, prodigal_gene in enumerate(prodigal_genes):
		corrected = False
		database_cursor.execute('SELECT * FROM raw_features WHERE start <= "{0}"'.format(prodigal_gene[2]))
		database_cursor.execute('SELECT * FROM blast_hits WHERE prodigal_id = "{0}"'.format(prodigal_gene[0]))
		blast_hits = [(blast_hit[3],blast_hit[6],blast_hit[10]) for blast_hit in database_cursor.fetchall()]
		blast_hits_query = sorted(blast_hits,
		                          key=lambda hit:[
								             	blast_hit[0] for blast_hit in blast_hits
											 ].count(hit[0]))[::-1][0:5]
		blast_hits_hit = sorted(blast_hits,
		                          key=lambda hit:[
								                blast_hit[0] for blast_hit in blast_hits
										     ].count(hit[1]))[::-1][0:5]

		# correct the gene if it has a hit on antifam

		database_cursor.execute("SELECT * FROM raw_features WHERE data_id=%s AND type = 'orf'"%(prodigal_gene[0]))
		feature_data =  database_cursor.fetchone()
		if feature_data:
			if feature_data != 'False':
				for shorter_gene in reduce_cds(prodigal_gene,sequencia):
					if prodigal_gene[4] == '+':
						shorter_protein = Seq(sequencia[2][shorter_gene[2]-1:shorter_gene[3]]).translate()
					if prodigal_gene[4] == '-':
						shorter_protein = Seq(str(Seq(sequencia[2][shorter_gene[2]-1:shorter_gene[3]]).reverse_complement())).translate()
					temp_fasta = open('%s/temp.fasta'%arguments.output_dir,'w')
					temp_fasta.write(">{0} {1}\n{2}\n".format(prodigal_gene, shorter_gene,shorter_protein))
					temp_fasta.close()
					subprocess.call(('{0}/bin/hmmer3/src/./hmmscan --cpu {1} '
				   '-o {2}/temp.fasta {0}/databases/antifam/AntiFam.hmm '
				   '{2}/temp.fasta').format(script_location,
				   	arguments.threads,arguments.output_dir),
				   	shell=True,
				   	stdout=devnull,
				   	stderr=devnull)
					antifam_hit = False
					antifam_output = SearchIO.parse('%s/temp.fasta'%arguments.output_dir,'hmmer3-text')
					for resultado in antifam_output:
						for hsp in resultado.hsps:
							if hsp.evalue <= arguments.antifam_evalue:
								antifam_hit = True
					if not antifam_hit:
						break
				if not antifam_hit:
					corrected = True
					if prodigal_gene[4] == '+':
						database_cursor.execute("UPDATE raw_features SET start=%s,spurious='False' WHERE data_id=%s AND type = 'orf'"%(shorter_gene[2], prodigal_gene[0]))
					if prodigal_gene[4] == '-':
						database_cursor.execute("UPDATE raw_features SET end=%s,spurious='False' WHERE data_id=%s AND type = 'orf'"%(shorter_gene[3], prodigal_gene[0]))

					break

		# Reduce gene by shifting the 5' terminus

		if blast_hits_query[0][0] > 1:
			closest_gene = prodigal_gene
			closest_distance = blast_hits_query[0][0]
			for shorter_gene in reduce_cds(prodigal_gene,sequencia):
				protein_reduction = abs(float((prodigal_gene[3] - prodigal_gene[2] + 1) \
				                        - (shorter_gene[3] - shorter_gene[2] + 1)) / 3)
				distance = abs(blast_hits_query[0][0] - protein_reduction)
				if distance < closest_distance:
					closest_distance = distance
					closest_gene = shorter_gene
				else:
					break


		if prodigal_gene[4] == '+':
			database_cursor.execute('SELECT * FROM raw_features WHERE start < %s ORDER BY start'%(prodigal_gene[2]))
			downtream_features = database_cursor.fetchall()
			if downtream_features:
				extension_limit = downtream_features[-1][2]
			else:
				extension_limit = 1

		if prodigal_gene[4] == '-':
			database_cursor.execute('SELECT * FROM raw_features WHERE start > %s ORDER BY start'%(prodigal_gene[3]))
			downtream_features = database_cursor.fetchall()
			if downtream_features:
				extension_limit = downtream_features[0][2]
			else:
				extension_limit = len(sequencia[2])

		if not corrected and blast_hits_query[0][1] > 1:
			closest_gene = prodigal_gene
			closest_distance = blast_hits_query[0][0]
			initial_distance = closest_distance
			for longer_gene in extend_cds(prodigal_gene,sequencia,extension_limit):
				protein_reduction = abs(float((prodigal_gene[3] - prodigal_gene[2] + 1) \
				                             - (shorter_gene[3] - shorter_gene[2] + 1)) / 3)
				distance = abs(blast_hits_query[0][0] - protein_reduction)
				if distance < closest_distance:
					closest_distance = distance
					closest_gene = shorter_gene
				else:
					break

			if initial_distance > closest_distance:
				pass
		if not corrected and blast_hits_hit[0][1]:
			closest_gene = prodigal_gene
			closest_distance = blast_hits_hit[0][0]
			initial_distance = closest_distance
			antifam_hit = False
			for extended_gene in extend_cds(prodigal_gene,sequencia,limit=extension_limit):
				if prodigal_gene[4] == '+':
					extended_protein = Seq(sequencia[2][extended_gene[2]-1:extended_gene[3]]).translate()
				if prodigal_gene[4] == '-':
					extended_protein = Seq(str(Seq(sequencia[2][extended_gene[2]-1:extended_gene[3]]).reverse_complement())).translate()
				temp_fasta = open('%s/temp.fasta'%arguments.output_dir,'w')
				temp_fasta.write(">{0} {1}\n{2}\n".format(prodigal_gene, extended_gene,extended_protein))
				temp_fasta.close()
				subprocess.call(('{0}/bin/hmmer3/src/./hmmscan --cpu {1} '
					           '-o {2}/temp.fasta {0}/databases/antifam/AntiFam.hmm '
					           '{2}/temp.fasta').format(script_location,
					           arguments.threads,arguments.output_dir),
					           shell=True,
					           stdout=devnull,
					           stderr=devnull)
				antifam_output = SearchIO.parse('%s/temp.fasta'%arguments.output_dir,
					                            'hmmer3-text')
				for resultado in antifam_output:
					for hsp in resultado.hsps:
						if hsp.evalue <= arguments.antifam_evalue:
							antifam_hit = True
				if not antifam_hit:
					database_cursor.execute("UPDATE raw_features SET spurious='True' WHERE data_id=%s"%(prodigal_gene[0]))
					break



#  -----------------------------------------------------------------------------
#  Identify frameshifted genes
#  -----------------------------------------------------------------------------
#
#
#  -----------------------------------------------------------------------------

sys.stdout.write('[Final processing]: Identifying frameshifted genes (pseudogenes).\n')

database_cursor.execute('SELECT * FROM sequences')
sequencias = database_cursor.fetchall()

for sequencia in sequencias:

	database_cursor.execute('SELECT * FROM raw_features WHERE sequence_id = %s AND type = "orf"'%sequencia[0])
	prodigal_genes = database_cursor.fetchall()

	for prodigal_gene_index, prodigal_gene in enumerate(prodigal_genes):
		if prodigal_gene_index < (len(prodigal_genes)-1):
			database_cursor.execute('SELECT * FROM blast_hits WHERE prodigal_id={0}'.format(prodigal_gene[6]))
			current_gene_blast_hits = database_cursor.fetchall()
			current_gene_hits_list = []
			database_cursor.execute('SELECT * FROM blast_hits WHERE prodigal_id={0}'.format(int(prodigal_gene[6])+1))
			next_gene_blast_hits = database_cursor.fetchall()
			next_gene_hits_list = []
			for hit in current_gene_blast_hits:
				current_gene_hits_list.append((hit[12],hit[6],hit[7],hit[8]))
			for hit in next_gene_blast_hits:
				next_gene_hits_list.append((hit[12],hit[6],hit[7],hit[8]))
			for hit_c in current_gene_hits_list:
				for hit_n in next_gene_hits_list:
					if hit_c[0] == hit_n[0]:
						hit_c_qc = set(range(hit_c[1],hit_c[2]))
						hit_n_qc = set(range(hit_n[1],hit_n[2]))
						if float(len(hit_c_qc.intersection(hit_n_qc)))/hit_c[3] < 0.1:
							database_cursor.execute("INSERT INTO putative_frameshifted VALUES (NULL,%s)"%(prodigal_gene[0]))
							database_cursor.execute("INSERT INTO putative_frameshifted VALUES (NULL,%s)"%(prodigal_genes[prodigal_gene_index+1][0]))
							database_connection.commit()

#  -----------------------------------------------------------------------------
#  Remove overlapping CDSs
#  -----------------------------------------------------------------------------
#  These families of RNAs are already predicted by other tools, which are more
#  purpose-specific, so they are not included when found by INFERNAL.
#  -----------------------------------------------------------------------------

sys.stdout.write('[Final processing]: Removing CDSs that overlap ncRNAs.\n')

database_cursor.execute('SELECT * FROM sequences ORDER BY id')
sequences = database_cursor.fetchall()

for sequence in sequences:

	database_cursor.execute('SELECT * FROM raw_features WHERE type = "orf" AND sequence_id = %s ORDER BY start'%sequence[0])
	orf_features = database_cursor.fetchall()
        database_cursor.execute('SELECT * FROM raw_features WHERE type <> "orf" AND sequence_id = %s ORDER BY start'%sequence[0])
        non_orf_features = database_cursor.fetchall()

	for non_orf_feature in non_orf_features:

		database_cursor.execute('INSERT INTO final_features VALUES (NULL,%s)'%non_orf_feature[0])

	for orf_feature in orf_features:

		if len(non_orf_features) > 0:

			for non_orf_feature in non_orf_features:

				antisense_rna = False
				if non_orf_feature[5] == 'ncrna':

					database_cursor.execute('SELECT * FROM infernal where id = %s'%(non_orf_feature[6]))
					non_orf_infernal_data = database_cursor.fetchone()

					if non_orf_infernal_data:
						database_cursor.execute('SELECT rfam_class FROM rfam where rfam_id = "%s"'%non_orf_infernal_data[-1])
						rfam_data = database_cursor.fetchone()

						if rfam_data:

							if rfam_data[0].lower().find('antisense'):

								antisense_rna = True

				orf_feature_range = set(range(orf_feature[2],orf_feature[3]+1))
				non_orf_feature_range = set(range(non_orf_feature[2],non_orf_feature[3]+1))
				intersection = orf_feature_range.intersection(non_orf_feature_range)

				if len(intersection) >= 3 and antisense_rna == False:

					database_cursor.execute('SELECT * FROM blast_hits WHERE prodigal_id={0}'.format(orf_feature[6]))
					hits = database_cursor.fetchall()
					database_cursor.execute("SELECT * FROM prodigal WHERE ID=%s"%orf_feature[6])
					prodigal_gene = database_cursor.fetchone()
					
					for shorter_gene in reduce_cds(prodigal_gene,sequence):
						new_orf_feature_range = set(range(shorter_gene[2],shorter_gene[3]+1))
						intersection = new_orf_feature_range.intersection(non_orf_feature_range)
						if not intersection:
							database_cursor.execute('SELECT * FROM blast_hits WHERE prodigal_id={0}'.format(orf_feature[6]))
							hits = database_cursor.fetchall()
							break
					
					#else:
					#
					#	break

			else:

				if orf_feature[7] == 'True':

					database_cursor.execute('INSERT INTO final_features VALUES (NULL,%s)'%str(orf_feature[0]))
					database_connection.commit()

		else:

			if orf_feature[7] == 'True':

				database_cursor.execute('INSERT INTO final_features VALUES (NULL,%s)'%str(orf_feature[0]))
				database_connection.commit()

gc.collect()

#  -----------------------------------------------------------------------------
#  Prepare the output files
#  -----------------------------------------------------------------------------
#  These families of RNAs are already predicted by other tools, which are more
#  purpose-specific, so they are not included when found by INFERNAL.
#  -----------------------------------------------------------------------------

sys.stdout.write('[Final processing]: Preparing output files.\n')

database_cursor.execute('SELECT * FROM final_features')
all_final_features = database_cursor.fetchall()
final_features = [entry[1] for entry in all_final_features]
number_features = len(final_features)
database_cursor.execute('SELECT * FROM sequences')
sequences = database_cursor.fetchall()
number_features = len(final_features)
feature_number = 0

gc.collect()

file_genbank = ''
file_table = ''
file_fasta = ''
locus_count = 0

for index,sequence in enumerate(sequences):

	sequence_id = arguments.scaffold_label + '00' + '0'*(len(str(len(sequences)))-len(str(index+1)))+str(index+1)
	file_fasta += '>%s\n%s\n'%(sequence_id,sequence[2])
	database_cursor.execute('SELECT * FROM raw_features WHERE sequence_id=%s ORDER BY start'%sequence[0])
	raw_features = database_cursor.fetchall()
	features = [raw_feature for raw_feature in raw_features if raw_feature[0] in final_features]
	sys.stdout.write('\t - Annotation for sequence "%s" ... '%(sequence_id))
	sys.stdout.flush()
	file_table += '>Features '+sequence_id+'\n'
	file_genbank += 'LOCUS       %s\n'%(sequence_id)
	file_genbank += 'DEFINITION  "%s" annotated by Genix\n'%(arguments.input)
	file_genbank += 'FEATURES             Location/Qualifiers\n'
	file_genbank += '     source          1..%s\n'%(len(sequence[2]))

	for feature_index, feature in enumerate(features):

		feature_number += 1
		feature_start = feature[2]
		feature_end = feature[3]

		if feature[4] == '+':

			strand = 1
			feature_sequence = sequence[2][feature_start-1:feature_end]

		if feature[4] == '-':

			strand = -1
			feature_sequence = Seq(sequence[2][feature_start-1:feature_end],generic_dna)
			feature_sequence = str(feature_sequence.reverse_complement())

		if feature[5] == 'orf':

			locus_count += 1
			locus_tag = arguments.locus_tag + '_' +'0'*(len(str(len(final_features)))-len(str(locus_count)))+str(locus_count)
			codon_start = 1
			partial_5 = False
			partial_3 = False
			database_cursor.execute('SELECT * FROM blast_hits WHERE prodigal_id = %s'%feature[6])
			blast_result = database_cursor.fetchone()
			protein_name = blast_result[2]
			uniprot_hit = blast_result[12]

			#analisa caso haja apenas uma feature na scaffold

			if len(features) == 1:

				if feature_sequence[0:3] not in start_codons:

					if feature[4] == '+':

						partial_5 = True
						codon_start = ((feature[3])%3)+1
						feature_start=1

					if feature[4] == '-':
						partial_5 = True
						codon_start = (len(sequence[2])-feature[3])%3+1
						feature_end = len(sequence[2])

				if feature_sequence[-3:] not in stop_codons:

					if feature[4] == '+':

						partial_3 = True
						feature_end = len(sequence[2])

					if feature[4] == '-':

						partial_3 = True
						feature_start = 1


			#analisa a primeira feature da scaffold

			elif feature_index == 0 and len(features) > 1:

				if feature_sequence[0:3] not in start_codons:

					if feature[4] == '+':

						partial_5 = True
						codon_start = ((feature_end)%3)+1
						feature_start=1

				if feature_sequence[-3:] not in stop_codons:

					if feature[4] == '-':

						partial_3 = True
						feature_start = 1

			#se for a ultima feature da scaffold

			elif feature_index+1 == len(features) and len(features) > 1:

				if feature_sequence[0:3] not in start_codons:

					if feature[4] == '-':

						partial_5 = True
						codon_start = (((len(sequence[2]) - feature_start + 1))%3)+1
						feature_end = len(sequence[2])

				if feature_sequence[-3:] not in stop_codons:

					if feature[4] == '+':
						partial_3 = True
						feature_end = len(sequence[2])

			# genes quebrados dentro da scaffold

			elif feature_index > 0 and feature_index < len(features)-1:

				if feature_sequence[0:3] not in start_codons:

					partial_5 = True

				if feature_sequence[-3:] not in stop_codons:

					partial_3 = True

				if feature[4] == '+':

					tbl_start = feature_start
					tbl_end = feature_end

				if feature[4] == '-':

					tbl_start = feature_end
					tbl_end = feature_start

				if partial_5 or partial_3:

					if feature[4] == '+':

						file_genbank += '     gene            %s..%s\n'%(feature_start,feature_end)

					else:

						file_genbank += '     gene            complement(%s..%s)\n'%(feature_start,feature_end)

					file_genbank += '                     /locus_tag="%s"\n'%(locus_tag)
					file_table += "%s\t%s\tgene\n"%(tbl_start,tbl_end)
					file_table += "\t\t\tlocus_tag\t%s\n"%(locus_tag)

					if partial_5 and not partial_3:

						file_genbank += '                     /note="\%s\' identified by Prodigal and BLAST with missing 5\' region."\n'%(protein_name[:-1])
						file_table += '\t\t\tnote\t\%s\ identified by Prodigal and BLAST with missing 5\' region.\n'%(protein_name[:-1])

					if partial_3 and not partial_5:

						file_genbank += '                     /note="\%s\' identified by Prodigal and BLAST with missing 3\' region."\n'%(protein_name[:-1])
						file_table += '\t\t\tnote\t\%s\ identified by Prodigal and BLAST with missing 3\' region.\n'%(protein_name[:-1])

					if partial_5 and partial_3:

						file_table += '\t\t\tnote\t\%s\ identified by Prodigal and BLAST with missing 5\' and 3\' regions.\n'%(protein_name[:-1])
						file_genbank += '                     /note="\%s\' identified by Prodigal and BLAST with missing 5\' and 3\' regions."\n'%(protein_name[:-1])
					continue


				else:

					if feature[4] == '+':

						file_genbank += '     gene            %s..%s\n'%(feature_start,feature_end)

					else:

						file_genbank += '     gene            complement(%s..%s)\n'%(feature_start,feature_end)

					file_genbank += '                     /locus_tag="%s"\n'%(locus_tag)

					if feature[4] == '+':

						file_genbank += '     CDS             %s..%s\n'%(feature_start,feature_end)

					else:

						file_genbank += '     CDS             complement(%s..%s)\n'%(feature_start,feature_end)

					file_genbank += '                     /product="%s"\n'%(protein_name[:-1])
					file_genbank += '                     /protein_id="gnl|%s|%s"\n'%(arguments.institution,locus_tag)
					file_genbank += '                     /codon_start=%s\n'%(codon_start)
					file_genbank += '                     /note="predicted by Prodigal"\n'
					file_genbank += '                     /note="Similar to Uniprot:%s"\n'%(uniprot_hit)

					database_cursor.execute('SELECT * FROM putative_frameshifted WHERE raw_feature_id =%s'%feature[0])
					if database_cursor.fetchall():
						file_genbank += '                     /note="This feature shares one or more BLAST hits with other neighbor features. We strongly recommend manual revision."\n'


					if strand == 1:

						tbl_start = str(feature_start)
						tbl_end = str(feature_end)

					else:

						tbl_start = str(feature_end)
						tbl_end = str(feature_start)

					file_table += "%s\t%s\tgene\n"%(tbl_start,tbl_end)
					file_table += "\t\t\tlocus_tag\t%s\n"%(locus_tag)
					file_table += "%s\t%s\tCDS\n"%(tbl_start,tbl_end)
					file_table += "\t\t\tproduct\t%s\n"%(protein_name)
					file_table += '\t\t\tprotein_id\tgnl|%s|%s\n'%(arguments.institution,locus_tag)
					file_table += "\t\t\tcodon_start\t%s\n"%(codon_start)
					file_table += "\t\t\tnote\tPredicted by Prodigal\n"
					file_table += '\t\t\tnote\tSimilar to Uniprot:%s\n'%(uniprot_hit)
					database_cursor.execute('SELECT * FROM putative_frameshifted WHERE raw_feature_id =%s'%feature[0])
					if database_cursor.fetchall():
						file_table += "\t\t\tnote\tThis feature shares one or more BLAST hits with other neighbor features. We strongly recommend manual revision\n"


				continue

			#add data to the output file

			if strand == 1:

				if partial_5:

					corrected_start = '<'+str(feature_start)

				else:

					corrected_start = str(feature_start)

				if partial_3:

					corrected_end = '>'+str(feature_end)

				else:

					corrected_end = str(feature_end)

			if strand == -1:

				if partial_5:

					corrected_start = 'complement(<'+str(feature_start)

				else:

					corrected_start = 'complement('+str(feature_start)

				if partial_3:

					corrected_end = '>'+str(feature_end)+')'

				else:

					corrected_end = str(feature_end)+')'

			# Certainly there must be a better way to construct a locus_tag, but
			# at least this is working.

			locus_tag = arguments.locus_tag + '_' +'0'*(len(str(len(final_features)))-len(str(locus_count)))+str(locus_count)

			file_genbank += '     gene            %s..%s\n'%(corrected_start,corrected_end)
			file_genbank += '                     /locus_tag="%s"\n'%(locus_tag)
			file_genbank += '     CDS             %s..%s\n'%(corrected_start,corrected_end)
			file_genbank += '                     /product="%s"\n'%(protein_name[:-1])
			file_genbank += '                     /protein_id="gnl|%s|%s"\n'%(arguments.institution,locus_tag)
			file_genbank += '                     /codon_start=%s\n'%(codon_start)
			file_genbank += '                     /note="predicted by Prodigal"\n'
			file_genbank += '                     /note="Similar to Uniprot:%s"\n'%(uniprot_hit)
			database_cursor.execute('SELECT * FROM putative_frameshifted WHERE raw_feature_id =%s'%feature[0])
			if database_cursor.fetchall():
				file_genbank += '                     /note="This feature shares one or more BLAST hits with other neighbor features. We strongly recommend manual revision."\n'

			if strand == 1:

				if partial_5:

					tbl_start = '<'+str(feature_start)

				else:

					tbl_start = str(feature_start)

				if partial_3:

					tbl_end = '>'+str(feature_end)

				else:

					tbl_end = str(feature_end)

			if strand == -1:

				if partial_5:

					tbl_start = '<'+str(feature_end)

				else:

					tbl_start = str(feature_end)

				if partial_3:

					tbl_end = '>'+str(feature_start)

				else:

					tbl_end = str(feature_start)

			file_table += "%s\t%s\tgene\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tlocus_tag\t%s\n"%(locus_tag)
			file_table += "%s\t%s\tCDS\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tproduct\t%s\n"%(protein_name)
			file_table += '\t\t\tprotein_id\tgnl|%s|%s\n'%(arguments.institution,locus_tag)
			file_table += "\t\t\tcodon_start\t%s\n"%(codon_start)
			file_table += "\t\t\tnote\tPredicted by Prodigal\n"
			file_table += "\t\t\tnote\tSimilar to Uniprot:%s\n"%(uniprot_hit)
			database_cursor.execute('SELECT * FROM putative_frameshifted WHERE raw_feature_id =%s'%feature[0])
			if database_cursor.fetchall():
				file_table += "\t\t\tnote\tThis feature shares one or more BLAST hits with other neighbor features. We strongly recommend manual revision\n"

		elif feature[5] == 'trna':

			locus_count += 1
			locus_tag = arguments.locus_tag + '_' + '0'*(len(str(len(final_features)))-len(str(locus_count)))+str(locus_count)
			database_cursor.execute('SELECT * FROM trnascan WHERE id = %s'%feature[6])
			trnascan_result = database_cursor.fetchone()

			if feature[4] == '+':

				file_genbank += '     gene            %s..%s\n'%(feature_start,feature_end)

			else:

				file_genbank += '     gene            complement(%s..%s)\n'%(feature_start,feature_end)

			file_genbank += '                     /locus_tag="%s"\n'%(locus_tag)

			if feature[4] == '+':

				file_genbank += '     tRNA            %s..%s\n'%(feature_start,feature_end)

			else:

				file_genbank += '     tRNA            complement(%s..%s)\n'%(feature_start,feature_end)

			file_genbank += '                     /product="tRNA-%s-%s"\n'%(trnascan_result[5],trnascan_result[6])
			file_genbank += '                     /note="predicted by tRNAscan-SE"\n'

			if feature[4] == '+':

				tbl_start = feature_start
				tbl_end = feature_end

			if feature[4] == '-':

				tbl_start = feature_end
				tbl_end = feature_start

			file_table += "%s\t%s\tgene\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tlocus_tag\t%s\n"%(locus_tag)
			file_table += "%s\t%s\ttRNA\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tproduct\ttRNA-%s-%s\n"%(trnascan_result[5],trnascan_result[6])
			file_table += "\t\t\tnote\tPredicted by tRNAscan-SE\n"

		elif feature[5] == 'rrna':

			locus_count += 1
			locus_tag = arguments.locus_tag + '_' + '0'*(len(str(len(final_features)))-len(str(locus_count)))+str(locus_count)
			database_cursor.execute('SELECT * FROM rnammer WHERE id = %s'%feature[6])
			rnammer_result = database_cursor.fetchone()

			if feature[4] == '+':

				file_genbank += '     gene            %s..%s\n'%(feature_start,feature_end)

			else:

				file_genbank += '     gene            complement(%s..%s)\n'%(feature_start,feature_end)

			file_genbank += '                     /locus_tag="%s"\n'%(locus_tag)

			if feature[4] == '+':

				file_genbank += '     rRNA            %s..%s\n'%(feature_start,feature_end)

			else:

				file_genbank += '     rRNA            complement(%s..%s)\n'%(feature_start,feature_end)

			file_genbank += '                     /product="%s"\n'%(rnammer_result[5])
			file_genbank += '                     /note="predicted by RNAmmer"\n'

			if feature[4] == '+':

				tbl_start = feature_start
				tbl_end = feature_end

			if feature[4] == '-':

				tbl_start = feature_end
				tbl_end = feature_start

			file_table += "%s\t%s\tgene\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tlocus_tag\t%s\n"%(locus_tag)
			file_table += "%s\t%s\trRNA\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tproduct\t%s\n"%(rnammer_result[5])
			file_table += "\t\t\tnote\tPredicted by RNAmmer\n"

		elif feature[5] == 'tmrna':

			locus_count += 1
			locus_tag = arguments.locus_tag + '_' + '0'*(len(str(len(final_features)))-len(str(locus_count)))+str(locus_count)
			database_cursor.execute('SELECT * FROM aragorn WHERE id = %s'%feature[6])

			if feature[4] == '+':

				file_genbank += '     gene            %s..%s\n'%(feature_start,feature_end)

			else:

				file_genbank += '     gene            complement(%s..%s)\n'%(feature_start,feature_end)

			file_genbank += '                     /locus_tag="%s"\n'%(locus_tag)

			if feature[4] == '+':
				file_genbank += '     tmRNA            %s..%s\n'%(feature_start,feature_end)


			else:

				file_genbank += '     tmRNA            complement(%s..%s)\n'%(feature_start,feature_end)

			file_genbank += '                     /note="predicted by aragorn"\n'

			if feature[4] == '+':

				tbl_start = feature_start
				tbl_end = feature_end

			if feature[4] == '-':

				tbl_start = feature_end
				tbl_end = feature_start

			file_table += "%s\t%s\tgene\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tlocus_tag\t%s\n"%(locus_tag)
			file_table += "%s\t%s\ttmRNA\n"%(tbl_start,tbl_end)
			file_table += "\t\t\tnote\tPredicted by Aragorn\n"

		elif feature[5] == 'ncrna':

			database_cursor.execute('SELECT * FROM infernal WHERE id = %s'%feature[6])
			infernal_result = database_cursor.fetchone()
			database_cursor.execute('SELECT * FROM rfam WHERE rfam_id = "%s"'%infernal_result[-1])
			rfam_data = database_cursor.fetchone()

			if 'Gene' in rfam_data[-1]:

				if rfam_data[4].split(';')[1] == '':

					continue

				else:

					locus_count += 1
					locus_tag = arguments.locus_tag + '_' + '0'*(len(str(len(final_features)))-len(str(locus_count)))+str(locus_count)

					if feature[4] == '+':

						file_genbank += '     gene            %s..%s\n'%(feature_start,feature_end)

					else:

						file_genbank += '     gene            complement(%s..%s)\n'%(feature_start,feature_end)

					file_genbank += '                     /locus_tag="%s"\n'%(locus_tag)

					if feature[4] == '+':

						file_genbank += '     ncRNA           %s..%s\n'%(feature_start,feature_end)

					else:

						file_genbank += '     ncRNA           complement(%s..%s)\n'%(feature_start,feature_end)

					file_genbank += '                     /product="%s"\n'%(rfam_data[3])
					file_genbank += '                     /ncRNA_class="%s"\n'%(rfam_data[4].split(';')[1])
					file_genbank += '                     /note="predicted by INFERNAL using RFam database."\n'

					if feature[4] == '+':

						tbl_start = feature_start
						tbl_end = feature_end

					if feature[4] == '-':

						tbl_start = feature_end
						tbl_end = feature_start

					file_table += "%s\t%s\tgene\n"%(tbl_start,tbl_end)
					file_table += "\t\t\tlocus_tag\t%s\n"%(locus_tag)
					file_table += "%s\t%s\tncRNA\n"%(tbl_start,tbl_end)
					file_table += "\t\t\tproduct\t%s\n"%(rfam_data[3])
					file_table += '\t\t\tncRNA_class\t%s\n'%(rfam_data[4].split(';')[1])
					file_table += "\t\t\tnote\tpredicted by INFERNAL using RFam database.\n"

			elif 'CRISPR RNA direct repeat element' in rfam_data[3]:

				if feature[4] == '+':
					file_genbank += '     repeat_region   %s..%s\n'%(feature_start,feature_end)
				else:

					file_genbank += '     repeat_region   complement(%s..%s)\n'%(feature_start,feature_end)

				file_genbank += '                     /note="CRISPR RNA direct repeat element"\n'
				file_genbank += '                     /note="%s (RFam:%s). Predicted by INFERNAL using RFam database."\n'%(rfam_data[3],rfam_data[1])

				if feature[4] == '+':

					tbl_start = feature_start
					tbl_end = feature_end

				if feature[4] == '-':

					tbl_start = feature_end
					tbl_end = feature_start

				file_table += "%s\t%s\trepeat_region\n"%(tbl_start,tbl_end)
				file_table += "\t\t\tnote\tCRISPR RNA direct repeat element\n"
				file_table += "\t\t\tnote\t%s (RFam:%s). Predicted by INFERNAL using RFam database.\n"%(rfam_data[3],rfam_data[1])

			elif 'riboswitch' in rfam_data[-1]:

				if feature[4] == '+':

					file_genbank += '     misc_feature    %s..%s\n'%(feature_start,feature_end)

				else:

					file_genbank += '     misc_feature    complement(%s..%s)\n'%(feature_start,feature_end)

				file_genbank += '                     /note="riboswitch"\n'
				file_genbank += '                     /note="%s (RFam:%s). Predicted by INFERNAL using RFam database."\n'%(rfam_data[3],rfam_data[1])

				if feature[4] == '+':

					tbl_start = feature_start
					tbl_end = feature_end

				if feature[4] == '-':

					tbl_start = feature_end
					tbl_end = feature_start

				file_table += "%s\t%s\tmisc_feature\n"%(tbl_start,tbl_end)
				file_table += "\t\t\tnote\triboswitch\n"
				file_table += "\t\t\tnote\t%s (RFam:%s). Predicted by INFERNAL using RFam database.\n"%(rfam_data[3],rfam_data[1])

			elif 'Cis-reg' in rfam_data[-1]:

				if feature[4] == '+':

					file_genbank += '     regulatory      %s..%s\n'%(feature_start,feature_end)

				else:

					file_genbank += '     regulatory      complement(%s..%s)\n'%(feature_start,feature_end)

				file_genbank += '                     /regulatory_class="other"\n'
				file_genbank += '                     /note="%s (RFam:%s). Predicted by INFERNAL using RFam database."\n'%(rfam_data[3],rfam_data[1])

				if feature[4] == '+':

					tbl_start = feature_start
					tbl_end = feature_end

				if feature[4] == '-':

					tbl_start = feature_end
					tbl_end = feature_start

				file_table += "%s\t%s\tregulatory\n"%(tbl_start,tbl_end)
				file_table += "\t\t\tnote\t%s (RFam:%s). Predicted by INFERNAL using RFam database.\n"%(rfam_data[3],rfam_data[1])

	file_genbank += 'ORIGIN'
	file_genbank += format_sequence(sequence[2])
	file_genbank += '\n//\n'
	sys.stdout.write('done!\n')

	#save the genbank output

	final_genbank_output = open('%s/annotation.gb'%arguments.output_dir,'w')
	final_genbank_output.write(file_genbank)
	final_genbank_output.close()

	#save the feature table

	final_table_output = open('%s/annotation.tbl'%arguments.output_dir,'w')
	final_table_output.write(file_table)
	final_table_output.close()

	#save the processed FASTA

	final_fasta_output = open('%s/annotation.fasta'%arguments.output_dir,'w')
	final_fasta_output.write(file_fasta)
	final_fasta_output.close()

#  -----------------------------------------------------------------------------
#  Genix:
#  -----------------------------------------------------------------------------
#  Current Features:
#  TODOS:
#  -----------------------------------------------------------------------------

os.system("{2} {0}/scripts/genbank2faa.py {1}/annotation.gb\
           {1}/annotation.faa".format(script_location,abs_output_dir,
           	                          PYTHON_PATH))

os.system("{2} {0}/scripts/genbank2ffn.py {1}/annotation.gb\
           {1}/annotation.ffn".format(script_location,abs_output_dir,
           	                          PYTHON_PATH))

os.system("{2} {0}/scripts/genbank2fna.py {1}/annotation.gb\
           {1}/annotation.fna".format(script_location,abs_output_dir,
           	                          PYTHON_PATH))
