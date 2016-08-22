# Genix:  An automated pipeline for bacterial genome annotation

## About

#### The annotation process

Genix is an online automated pipeline for bacterial genome annotation. The program takes a `FASTA file` containing a set of sequences (eg: complete chromosomes, contigs or scaffolds). The pipeline  uses a combination of several bioinformatics tools, including `Prodigal` (Hyatt et al. 2010), `BLAST`, `tRNAscan-SE` (Lowe & Eddy 1997), `RNAmmer` (Lagesen et al. 2007), `Aragorn` (Laslett 2004), `HMMER` (Eddy 2011), `INFERNAL` (Nawrocki et al. 2009), `RFam` (Griffiths-Jones et al. 2003) and `Antifam` (Eberhardt et al. 2012).

#### the *get_database.sh* script

<<<<<<< HEAD
<<<<<<< HEAD
The user may also create an optimized protein database for a given taxonomic group of interest by using the `get_database.sh` by providing it's tax_id. First, a dataset of proteins associated to the tax_id Is downloaded from `Uniprot` and used to build a raw dataset, which may contain several redundances. `CD-HIT` (Li & Godzik 2006) is used to build a non-redundant dataset, which is paserd by an python script to identify the best annotation from each cluster based on the protein-score and the source from where the proteins were obtained (Uniprot-Swissprot or Uniprot-trEMBL). Finally, the refined dataset is used to generate the final protein `BLAST` database (Altschul et al. 1990; Camacho et al. 2009).
=======
The user may also create an optimized protein database for a given taxonomic group of interest by using the `get_database.sh` by providing it's tax_id. First, a dataset of proteins associated to the tax_id Is downloaded from `Uniprot` and used to build a raw dataset, which may contain several redundances. `CD-HIT` (Li & Godzik 2006) is used to build a non-redundant dataset, which is paserd by an python script to identify the best annotation from each cluster based on the protein-score and the source from where the proteins were obtained (Uniprot-Swissprot or Uniprot-trEMBL). Finally, the refined dataset is used to generate the final protein `BLAST` database (Altschul et al. 1990; Camacho et al. 2009). 
>>>>>>> 959a9042d0fbe2e79480a8ca74146db34ad45b6d
=======
The user may also create an optimized protein database for a given taxonomic group of interest by using the `get_database.sh` by providing it's tax_id. First, a dataset of proteins associated to the tax_id Is downloaded from `Uniprot` and used to build a raw dataset, which may contain several redundances. `CD-HIT` (Li & Godzik 2006) is used to build a non-redundant dataset, which is paserd by an python script to identify the best annotation from each cluster based on the protein-score and the source from where the proteins were obtained (Uniprot-Swissprot or Uniprot-trEMBL). Finally, the refined dataset is used to generate the final protein `BLAST` database (Altschul et al. 1990; Camacho et al. 2009). 
>>>>>>> refs/remotes/origin/master

## Installation

### clonning the repository

`git clone https://github.com/fredericokremer/genix/`

`cd genix`

### Editing the source code

Edit the variables `PERL_PATH` (path to the Perl executable), `PYTHON_PATH` (path to
the Python executable) and `RNAMMER\_PATH` (path to the RNAmmer directory) are
correctly specified. Check if the binaries Prodigal (“bin/prodigal/”), Aragorn (“bin/Aragorn”), hmmer3 (“bin/hmmer3”), INFERNAL (“bin/infernal/src”), BLAST (“bin/ncbi-blast”) and tRNAscan-SE (“/bin/tRNAscan-SE”) are set as executable files. If not, you can use the command `chmod +x DIR/*` in each of the directories to provide executable permission. In some cases it may necessary to recombine some of these binaries, especially if you running then in a non-64 bits platform.

## Running Genix

### preparing the database

### Running the annotation

`python genix_annotation.py`

## Webserver

Genix is also available as an webserver through the URL [http://labbioinfo.ufpel.edu.br/genix](http://labbioinfo.ufpel.edu.br/genix).
<<<<<<< HEAD
<<<<<<< HEAD

##  Contact and Feedback

In case of any problem, or if you have any idea to improve the pipeline, please, contact us!
Frederico Schmitt Kremer,
fred.s.kremer@gmail.com / fredericok.cdtec@ufpel.edu.br
Luciano da Silva Pinto,
dmpluc@ufpel.edu.br / ls_pinto@hotmail.com
=======
>>>>>>> 959a9042d0fbe2e79480a8ca74146db34ad45b6d
=======
>>>>>>> refs/remotes/origin/master
