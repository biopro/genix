#!/usr/bin/env python

import json
import argparse
import sys

try:
	import requests
except:
	print 'GENIX Annotation Web API requires the "requests" package'
	sys.exit()

############################################################
#////////// Genix Annotation Web API //////////////////////#
############################################################
# Unfinished version 
#
'''
This script is just an example and features that are used here
may be replaced or updated in future releases of the Genix
API. Therefore, it's important to always check if your script
is in accordance with the most recent version of the API.
'''

# SUBMIT A GENOME

submission_data = {
    # the email of your account
    'user_email' : 'fred.s.kremer@gmail.com', #<------ replace by your user name
    #you password
    'user_password' : '', # <------ add your password right here
    #an string or your FASTA filte
    #do not use parsed objects (eg: "SeqIO.parse" objects)
    'genome_FASTA' : open('/home/cdtec/Frederico/genix/genix_annotation/'
                          'versions/0.1/webAPI/leptospira.fasta').read(),
    #the locus_tag prefix that will be assigned to your gene
    'locus_tag_prefix':'GNX',
    #the sequencing_center that will be assigned to your protein_ids
    'sequencing_center':'GNX',
    #the genetic_code table that will be used (11 or 4)
    'genetic_code_table':'11',
    #information for the database construction
    'database_construction':{
        #tax id: in this case it is necessary to provide the tax_id
        'tax_id':'171',
        #the identity threshold that will be used in the clustering
        'clustering':'95'
    },
    # metadata about your organism
    'organism' : {
        'organism' : 'Leptospira interrogans',
        'serovar' : '',
        'serogroup' : '',
        'serotype' : '',
        'isolate' : '',
        'isolation_source':'',
        'strain':'',
        'host':''
    },
    # the evalue threshold for each searching program
    'evalue_thresholds':{
        'protein_blast':0.0005,  # float number, usually between 0 and 1
        'rfam_blast':0.0005,     # float number, usually between 0 and 1
        'antifam_hmmer':0.0005,  # float number, usually between 0 and 1
        'rfam_infernal':0.0005}
    }

return_data = requests.post("http://labbioinfo.ufpel.edu.br/cgi-bin/genix_api.py",
                            data={'mode':'submission','data': json.dumps(
                                  submission_data)})

#print return_data.text
job_id =  json.loads(return_data.text)['job_id']

# CHECK STATUS

get_status_data = {
    # the email of your account
    'user_email' : 'fred.s.kremer@gmail.com',
    #you password
    'user_password' : 'ZX',
    'job_id':job_id
}

return_data = requests.post("http://labbioinfo.ufpel.edu.br/cgi-bin/genix_api.py",
                            data={'mode':'get_status','data': json.dumps(
                                  get_status_data)})
#print return_data.text

# GET RESULTS

get_results_data = {
    # the email of your account
    'user_email' : 'fred.s.kremer@gmail.com',
    #you password
    'user_password' : 'ZX',
    'job_id':'2NG9KIY0S60OFQP7DNSR'
}

result_data = requests.post("http://labbioinfo.ufpel.edu.br/cgi-bin/genix_api.py",
                            data={'mode':'get_results','data': json.dumps(
                                  get_results_data)})

# IT'S is just an preliminary implementation and currently on the
# Genbank annotation file can be retrieved. 

genbank_file = json.loads(result_data.text)['file_genbank']
