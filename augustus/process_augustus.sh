##-

## --species is trained species --workingdir is tmp folder
mkdir -p /40t_1/caix/annotation_A01/augustus/autoAugPred
~/miniconda3/envs/py3.5/bin/autoAugPred.pl  --genome=hn53.genome.scf.fasta.masked --species=BrapaV1.5  --workingdir=/40t_1/caix/annotation_A01/augustus/autoAugPred &

##- parallel run
ls aug* > augustus.file
perl run_augu.pl augustus.file > commands.list
perl run_augu.pl augustus.file > commands.list

##- cat results
cat aug*out > all.augustus.gff3

##  join results
perl  ~/miniconda3/envs/py3.5/bin/join_aug_pred.pl   <all.augustus.gff3>  augustus.gff3 


