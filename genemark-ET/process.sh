#!/bin/bash

#export PATH=/home/caix/miniconda3/envs/py2.7/bin/:${PATH}
#export PERL5LIB=/home/caix/miniconda3/envs/py2.7/lib/5.26.2/x86_64-linux-thread-multi:$PERL5LIB

#/40t_1/caix/software/biosoft/geneAnnotation/gm_et_linux_64/gmes_petap/star_to_gff.pl  --star  SJ.out.tab --gff SJ.gff --label STAR

perl /40t_1/caix/software/biosoft/geneAnnotation/gm_et_linux_64/gmes_petap/gmes_petap.pl   --sequence hn53.scf_masked.fasta  --ET  SJ.gff  --cores 30
