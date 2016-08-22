import sys
import csv
import re
import commands
from Bio import SeqIO

def get_clusters_cdhit(prefix):
	cluster_list = []
	clusterFileHandle = open('{0}_clustered.clstr'.format(
		                      prefix))
	cluster_current = []
	for linha in clusterFileHandle:
		if re.match('^>Cluster.+$',linha):
			if cluster_current:
				cluster_list.append(cluster_current)
				cluster_current = []
		elif re.match('^[0-9]+.+$',linha):
			uniprot_db = linha.split('>')[1].split('|')[0]
			uniprot_id = linha.split('>')[1].split('|')[1]
			if ' at ' in linha:
				cluster_id = float(linha.split(' at ')[1].replace('%',''))*0.01
			else:
				cluster_id = 1
			cluster_current.append([uniprot_db, uniprot_id, cluster_id])
	if cluster_current:
		cluster_list.append(cluster_current)
	return cluster_list

def refine_clusters(cluster_list,prefix):
	final_list = []
	csv_parser = csv.reader(open('{0}_raw.tab'.format(prefix)),delimiter='\t')
	uniprot_dict = {}
	for linha_i, linha_data in enumerate(csv_parser):
		if linha_i:
			pass
		else:
			continue
		uniprot_dict[linha_data[0]] = (linha_data[2], linha_data[1], linha_data[3], int(linha_data[4].split(' ')[0]))
		#print uniprot_dict
	for cluster in cluster_list:
		best_annotation = None
		for member in reversed(sorted(cluster,key=lambda member:member[2])):
			if not best_annotation:
				best_annotation = member
			else:
				if uniprot_dict[member[1]][1] == 'reviewed' and uniprot_dict[best_annotation[1]][1] == 'unreviewed':
					best_annotation = member
				if uniprot_dict[member[1]][3] > uniprot_dict[best_annotation[1]][3]:
					best_annotation = member
		final_list.append(best_annotation[1])
	return final_list
		
prefix = sys.argv[1]
lista = get_clusters_cdhit(prefix)

refined_list = refine_clusters(lista,prefix)

fasta_raw_handle = open('{0}_raw.fasta'.format(prefix))
fasta_raw_parser = SeqIO.parse(fasta_raw_handle,'fasta')
sys.stdout = open('{0}_clustered.corrected'.format(prefix),'w')

for entry in fasta_raw_parser:
	if entry.id.split('|')[1] in refined_list:
		sys.stdout.write(entry.format('fasta')+'\n')
		sys.stdout.flush()

