## Tool for analyzing Pfam-domain compositions of gene families
This tool assign Pfam domains to protein sequences in fasta files of gene families. The domains are reported in order of their starting positions in the sequence. Domain compositions of each given family is summarised in terms of percentage of sequences containing the domains found in the corresponding family and an average jaccard score calculated from all the pair-wise domain sequences in the family.

#### Requirments
* UNIX based OS
* PYTHON 2.7
* PYTHON modules _re_, _sys_, _os_, _subprocess_, _argparse_, _operator_
* _Pfam-A_ database files _Pfam-A.hmm_ and _Pfam-A.hmm.dat_ files saved in the _Pfam\/database_ directory. The _Pfam-A.hmm_ file must be pressed using _hmmpress_ program from the _HMMER_ package to create the database index file in the same directory.
* The _scripts_ and the _Pfam_ must be present in the same working directory
