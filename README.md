## Tool for analyzing Pfam-domain compositions of gene families
This tool assigns Pfam-domains to families of protein sequences. The domains are reported in order of their starting positions in the sequence. Domain compositions of each given family is summarised in terms of percentage of sequences containing the domains found in the corresponding family and average jaccard score calculated from all the pair-wise jaccard scores calculated using domain sequences from the family.

#### Requirments
* UNIX based OS
* PYTHON 2.7
* PYTHON modules _re_, _sys_, _os_, _subprocess_, _argparse_, _operator_
* _Pfam-A_ database files _Pfam-A.hmm_ and _Pfam-A.hmm.dat_ files saved in the _Pfam\/database_ directory. The _Pfam-A.hmm_ file must be pressed using _hmmpress_ program from the _HMMER_ package to create the database index file in the same directory
* The _scripts_ and the _Pfam_ directories must be present in the same working directory
* Add the Pfam modules to your PERL5LIB using the following command _bash\% export PERL5LIB\=\/path\/to\/pfam\_Dir\:$PERL5LIB_

#### Output files
* pfamscan\_out directory contains raw pfamscan output files for each family.
* domain\_order\_results directory contains domain order files for each family. Each domain order file contains Pfam-domains for each sequence in the family in order of their starting positions in the sequences. A _\*\*NULL\*\*_ domain is reported for sequences where no Pfam-domain is detected.
* \*.family\_domain\_compositions file contains summarised domain compositions for each family. Format\: \<family\_id\> \<domain\-1\>\-\<\% of sequences in the family containing the domain\> ... 
* \*.family\_domain\_jaccard\_scores file contains domain composition Jaccard scores for all the families. Format \<family\_id\> \<Jaccard score\>
