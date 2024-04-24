#!/usr/bin/bash

genome="hn53.scf_masked.fasta"
protein="genewise.gff3.form"
transcript="pasa.train.gff"
gene_predictions="augustus.form.gff3  GeneMark-ET.gff3.form"
threads=50
final_evm="final.evm.gff3"

##----------------------------------------------------------------------------------------------------
BIN=/40t_1/caix/annotation_pl/PASA_process
EVM=/40t_1/caix/software/gene_prediction/EVidenceModeler
##----------------------------------------------------------------------------------------------------
export PERL5LIB=/40t_1/caix/software/gene_prediction/EVidenceModeler/PerlLib/:${PERL5LIB}

echo "Step 1 : merge each evidence ..."
#cat ${protein}          > protein_alignments.gff
#cat ${transcript}       > transcript_alignments.gff
#cat ${gene_predictions} > gene_predictions.gff

echo "process weights.txt ..."
#cp $BIN/weights.txt ./

echo "process partition_EVM_inputs.pl ..."
#$EVM/EvmUtils/partition_EVM_inputs.pl  --genome  ${genome}  --gene_predictions  gene_predictions.gff  --protein_alignments  protein_alignments.gff  --transcript_alignments  transcript_alignments.gff  --segmentSize  10000000  --overlapSize  10000  --partition_listing  partitions_list.out

echo "process write_EVM_commands.pl ..."
$EVM/EvmUtils/write_EVM_commands.pl    --genome  ${genome}  --gene_predictions  gene_predictions.gff  --protein_alignments  protein_alignments.gff  --transcript_alignments  transcript_alignments.gff  --weights  `pwd`/weights.txt  --output_file_name  evm.out  --partitions  partitions_list.out  >  commands.lis

echo "process execute_EVM_commands.pl ..."
#$EVM/EvmUtils/execute_EVM_commands.pl  commands.lis
ParaFly -c commands.lis -CPU ${threads}

echo "process recombine_EVM_partial_outputs.pl ..."
$EVM/EvmUtils/recombine_EVM_partial_outputs.pl  --partitions  partitions_list.out  --output_file_name  evm.out

echo "process convert_EVM_outputs_to_GFF3.pl ..."
$EVM/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions  partitions_list.out  --output evm.out  --genome  ${genome}

echo "merge *evm.out.gff3 ..."
cat  */evm.out.gff3  >  evm.out

echo "process evm_to_gff.pl ..."
$BIN/evm_to_gff.pl  evm.out  -o  ${final_evm}

echo " all done!!!"
