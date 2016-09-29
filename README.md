# Genix:  An automated pipeline for bacterial genome annotation

## About

#### The annotation process

Genix is an online automated pipeline for bacterial genome annotation. The program takes a `FASTA file` containing a set of sequences (eg: complete chromosomes, contigs or scaffolds). The pipeline  uses a combination of several bioinformatics tools, including `Prodigal` (Hyatt et al. 2010), `BLAST`, `tRNAscan-SE` (Lowe & Eddy 1997), `RNAmmer` (Lagesen et al. 2007), `Aragorn` (Laslett 2004), `HMMER` (Eddy 2011), `INFERNAL` (Nawrocki et al. 2009), `RFam` (Griffiths-Jones et al. 2003) and `Antifam` (Eberhardt et al. 2012).

#### the *get_database.sh* script

The user may also create an optimized protein database for a given taxonomic group of interest by using the `get_database.sh` by providing it's tax_id. First, a dataset of proteins associated to the tax_id Is downloaded from `Uniprot` and used to build a raw dataset, which may contain several redundances. `CD-HIT` (Li & Godzik 2006) is used to build a non-redundant dataset, which is paserd by an python script to identify the best annotation from each cluster based on the protein-score and the source from where the proteins were obtained (Uniprot-Swissprot or Uniprot-trEMBL). Finally, the refined dataset is used to generate the final protein `BLAST` database (Altschul et al. 1990; Camacho et al. 2009).

The user may also create an optimized protein database for a given taxonomic group of interest by using the `get_database.sh` by providing it's tax_id. First, a dataset of proteins associated to the tax_id Is downloaded from `Uniprot` and used to build a raw dataset, which may contain several redundances. `CD-HIT` (Li & Godzik 2006) is used to build a non-redundant dataset, which is paserd by an python script to identify the best annotation from each cluster based on the protein-score and the source from where the proteins were obtained (Uniprot-Swissprot or Uniprot-trEMBL). Finally, the refined dataset is used to generate the final protein `BLAST` database (Altschul et al. 1990; Camacho et al. 2009).

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

##  Contact and Feedback

In case of any problem, or if you have any idea to improve the pipeline, please, contact us!

Frederico Schmitt Kremer,

fred.s.kremer@gmail.com / fredericok.cdtec@ufpel.edu.br

Luciano da Silva Pinto,

dmpluc@ufpel.edu.br / ls_pinto@hotmail.com

## References

 Altschul, S.F. et al., 1990. Basic local alignment search tool. Journal of molecular biology, 215(3), pp.403–10.

 Camacho, C. et al., 2009. BLAST+: architecture and applications. BMC bioinformatics, 10(1), p.421.

 Eberhardt, R.Y. et al., 2012. AntiFam: a tool to help identify spurious ORFs in protein annotation. Database : the journal of biological databases and curation, 2012, p.bas003.

 Eddy, S.R., 2011. Accelerated Profile HMM Searches. W. R. Pearson, ed. PLoS computational biology, 7(10), p.e1002195.

 Griffiths-Jones, S. et al., 2003. Rfam: an RNA family database. Nucleic acids research, 31(1), pp.439–41.

 Hyatt, D. et al., 2010. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics, 11, p.119.

 Lagesen, K. et al., 2007. RNAmmer: consistent and rapid annotation of ribosomal RNA genes. Nucleic acids research, 35(9), pp.3100–8.

 Laslett, D., 2004. ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. Nucleic Acids Research, 32(1), pp.11–16.

 Li, W. & Godzik, A., 2006. Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics (Oxford, England), 22(13), pp.1658–9.

 Lowe, T.M. & Eddy, S.R., 1997. tRNAscan-SE: a program for improved detection of transfer RNA genes in genomic sequence. Nucleic acids research, 25(5), pp.955–64.

 Nawrocki, E.P., Kolbe, D.L. & Eddy, S.R., 2009. Infernal 1.0: inference of RNA alignments. Bioinformatics (Oxford, England), 25(10), pp.1335–7.