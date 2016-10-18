#!/usr/bin/env bash

##### DATABASE PREPARATION

## :: create a directory for the database:

if [ ! -d database ]; then
	mkdir database
fi

## :: run the "get_database.sh" script
#
#  this script takes 3 arguments: [1] the tax_id of the organism 
#  of interest; [2] the path to the final database; [3] the identity 
#  threshold used to clusterize the database and reduce de redundance.
#

bash ../scripts/get_database.sh 561 database/e.coli 0.95

##### RUN THE ANNOTATION PIPELINE

## :: create a directory for the database:

if [ ! -d output ]; then
	mkdir output
fi

## :: run the "genix_annotation.py" script

python ../genix_annotation.py -i e.coli_K12.fasta -dbp database/e.coli \
       -o output/ --threads 4


