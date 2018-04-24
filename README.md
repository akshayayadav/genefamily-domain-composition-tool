## Tool for analyzing Pfam-domain compositions of gene families
This tool assigns Pfam-domains to families of protein sequences. The domains are reported in order of their starting positions in the sequence. Domain composition for each family is summarised in terms of percentage of sequences containing the domains found in the corresponding family and average jaccard score calculated from all the pair-wise jaccard scores from the family. Jaccard score between any two family sequences is = no-of-unique-common-domains-between-seq1-&-seq2/((total-no-of-unique-domains-in-seq1)+(total-no-of-unique-domains-in-seq2)) 

#### Requirments
* UNIX based OS
* PYTHON 2.7
* PYTHON modules _re_, _sys_, _os_, _subprocess_, _argparse_, _operator_
* _Pfam-A_ database files _Pfam-A.hmm_ and _Pfam-A.hmm.dat_ files saved in the _Pfam\/database_ directory. The _Pfam-A.hmm_ file must be pressed using _hmmpress_ program from the _HMMER_ package to create the database index file in the same directory
* The _scripts_ and the _Pfam_ directories must be present in the same working directory
* Add the Pfam modules to your PERL5LIB using the following command:
```
bash% export PERL5LIB=/path/to/pfam_Dir:$PERL5LIB
```

#### Tool Options
```
usage: gene-family-domain-arch-analysis.py [-h] --fasta_dir
                                           FAMILY_FASTA_FILE_DIR --output_dir
                                           OUTPUT_DIRECTORY --name
                                           DATASET_NAME

Tool for calculating domain composition and domain compostion scores for given
set of gene families

Arguments:
  -h, --help            show this help message and exit
  --fasta_dir FAMILY_FASTA_FILE_DIR
                        Location of directory containing family fasta files
  --output_dir OUTPUT_DIRECTORY
                        Location of the output directory
  --name DATASET_NAME   Name of the dataset used to create output files

```

#### Output files
* pfamscan\_out directory contains raw pfamscan output files for each family.
* domain\_order\_results directory contains domain order files for each family. Each domain order file contains Pfam-domains for each sequence in the family in order of their starting positions in the sequences. A _\*\*NULL\*\*_ domain is reported for sequences where no Pfam-domain is detected.
* \*.family\_domain\_compositions file contains summarised domain compositions for each family. Format\: \<family\_id\> \<family-size\> \<domain\-1\>\-\<\% of sequences in the family containing the domain\> ... 
* \*.family\_domain\_jaccard\_scores file contains domain composition Jaccard scores for all the families. Format \<family\_id\> \<family-size\> \<Jaccard score\>


