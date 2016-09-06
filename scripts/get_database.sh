#!/usr/bin/env bash


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

GENIX_DIR=$(readlink -f $(dirname $0)/../)

echo '--------------------------------------------------------------------------'
echo 'Genix: Get Database'
echo '--------------------------------------------------------------------------'

if [ ! -n $1 ] && [ ! -n $2 ] && [ ! -n $3 ]; then
  echo 'get_database.sh requires 3 arguments:'
  echo ' - TaxID'
  echo ' - Database Prefix'
  echo ' - Identity clustering threshold'
  echo ''
  echo 'USAGE: $ bash get_database.sh <tax id> <database prefix> <identity threshold>'
  exit
fi

echo 'Arguments:'
echo 'TAX_ID:' $1
echo 'PREFIX:' $2
echo 'IDENTI:' $3

echo 'Downloading Database'
wget_url='http://www.uniprot.org/uniprot/?query=taxonomy:'$1'+NOT+name:ORF+NOT+name:contig+NOT+name:scaffold+NOT+name:"genome%20shotgun"+NOT+name:"DNA%20for"+NOT+name:complete+NOT+name:partial&sort=score'
wget_url='http://www.uniprot.org/uniprot/?query=taxonomy:'$1'+NOT+name:ORF+NOT+name:contig+NOT+name:scaffold+NOT+name:"genome%20shotgun"+NOT+name:"DNA%20for"+NOT+name:complete+NOT+name:partial&sort=score'
wgetoutput=$2\_raw

wget -O $wgetoutput'.fasta' $wget_url'&format=fasta' 2> $2_get_database.log
wget -O $wgetoutput'.tab' $wget_url'&format=tab&columns=id,reviewed,protein names,genes,annotation score' 2>> $2_get_database.log

echo 'Running CD-HIT'
cdhit -i $wgetoutput.fasta -c $3 -o $2\_clustered -T 8 -M 4000 #>> $2_get_database.log
echo 'Refining clustering'
python $GENIX_DIR/scripts/get_database_optimizer.py $2 cdhit #>> $2_get_database.log
mv $2_clustered.corrected $2\_clustered
$GENIX_DIR/bin/ncbi-blast/./makeblastdb -in $2\_clustered -dbtype prot -out $2 >> $2_get_database.log

rm $2\_raw.fasta
