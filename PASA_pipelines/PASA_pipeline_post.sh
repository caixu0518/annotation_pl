#!/bin/bash
##- by Xu Cai
## must run in 23 servier

SPECIES=$1
#maskedGenome=$2
#transcripts=$3
homolog_protein=$2

##--------------------------------------------------------------------------------------------------------------------------------
BLAST_F4="qseqid  sseqid  qlen  slen  evalue  pident  ppos  length  mismatch  gapopen  qstart  qend  sstart  send  bitscore"
blastp=/40t_1/caix/software/biosoft/ncbi-blast-2.2.25+/bin/blastp
#GMAP=/data/caix/miniconda3/envs/PASA/bin
#BLAT=/40t_1/caix/software/gene_prediction/PASA/blat/bin
#FASTA=/40t_1/caix/software/gene_prediction/PASA/fasta-35.4.12/bin
#PASA=/40t_1/caix/software/gene_prediction/PASA/PASApipeline-2.0.2

BIN=/40t_1/caix/annotation_pl/PASA_process
##--------------------------------------------------------------------------------------------------------------------------------

#export  PATH=${GMAP}:${BLAT}:${FASTA}:${PASA}:$PASA/seqclean/seqclean/bin/:${PATH}
#export  PASAHOME=$PASA
#export  PERL5LIB=/40t_1/caix/software/tools/PerlLib/:/40t_1/caix/software/gene_prediction/PASA/PASApipeline-2.0.2/PerlLib/:$PERL5LIB
#export  PERL5LIB=/data/caix/perl5/lib/perl5/x86_64-linux-thread-multi/:$PERL5LIB

##--------------------------------------------------------------------------------------------------------------------------------
#echo "Step 1 ...cp pasa.alignAssembly.Template.txt..."
#cp  $PASA/pasa_conf/pasa.alignAssembly.Template.txt  alignAssembly.config

#echo "Step 2 ...vim alignAssembly.config..."
#perl  -p -i -e  "s/MYSQLDB=.*/MYSQLDB=pasa_$SPECIES/"  alignAssembly.config

#echo "Step 3 ...Launch_PASA_pipeline.pl..."
#perl $PASA/scripts/Launch_PASA_pipeline.pl  -g  ${maskedGenome}   -t  ${transcripts}  -c  alignAssembly.config  -C  -R  --ALIGNERS  blat,gmap  --CPU  20  --stringent_alignment_overlap 30.0

#echo "Step 4 ...pasa_asmbls_to_training_set.dbi..."
#perl $PASA/scripts/pasa_asmbls_to_training_set.dbi  --pasa_transcripts_fasta  pasa_$SPECIES.assemblies.fasta  --pasa_transcripts_gff3  pasa_${SPECIES}.pasa_assemblies.gff3

#echo "Step 5 ...pasa_to_gff.pl..."
#perl  $BIN/pasa_to_gff.pl  pasa_${SPECIES}.assemblies.fasta.transdecoder.genome.gff3  -o  pasa_${SPECIES}.assemblies.fasta.transdecoder.genome.gff

#echo "Step 6 ...gene_stats.pl ..."
#perl $BIN/gene_stats.pl pasa_${SPECIES}.assemblies.fasta.transdecoder.genome.gff > stat.out

#echo "Step 7 ...pasa_asmbls_to_training_set.extract_reference_orfs.pl..."
#$PASA/scripts/pasa_asmbls_to_training_set.extract_reference_orfs.pl  pasa_$SPECIES.assemblies.fasta.transdecoder.genome.gff3  100  >  pasa.orfs.gff

echo "Step 8 ...pasa_filter.pl..."
perl  $BIN/pasa_filter.pl  pasa.orfs.gff  -o  pasa.filter.gff  -cds  2  -p  pasa_$SPECIES.assemblies.fasta.transdecoder.pep  -c  pasa_${SPECIES}.assemblies.fasta.transdecoder.cds

echo "Step 9 ...run blastp..."
makeblastdb -in  ${homolog_protein}  -parse_seqids -hash_index  -out  ${homolog_protein}  -dbtype  prot
${blastp}  -evalue  1e-10   -query  pasa.filter.pep  -out  pasa.filter.blast  -db  ${homolog_protein}   -max_target_seqs  10  -num_threads 30   -outfmt  "7  $BLAST_F4"

echo "Step 10 ...pasa_filter1.pl..."
perl  $BIN/pasa_filter1.pl  pasa.filter.gff  -o  pasa.train.gff  -b  pasa.filter.blast  -p  pasa_${SPECIES}.assemblies.fasta.transdecoder.pep  -c  pasa_${SPECIES}.assemblies.fasta.transdecoder.cds  -min  100  -max  1000

echo "All done ..."
