#!/usr/bin/perl -w
use strict;

my $genome = $ARGV[0]; # genome fasta file

my $path = `pwd`;
chomp($path);

mkdir "$path/RepeatModeler" unless -e "$path/RepeatModeler";
chdir "$path/RepeatModeler";
system("ln -s $path/$genome .");
system("/disk02/caix/src/RepeatModeler-open-1.0.10/BuildDatabase -name $genome  -engine ncbi $genome &> BuildDatabase.log");
system("/disk02/caix/src/RepeatModeler-open-1.0.10/RepeatModeler -database $genome  -engine ncbi -pa 25 &> RepeatModeler.log");

mkdir "$path/RepeatMasker" unless -e "$path/RepeatMasker";
chdir "$path/RepeatMasker";
system("ln -s $path/$genome .");
system("cp $path/RepeatModeler/RM_\*/\*.classified .");
mkdir "Repeat_result" unless -e "Repeat_result";
system("/disk02/caix/src/RepeatMasker/RepeatMasker -e ncbi -pa 25 -lib consensi.fa.classified -dir Repeat_result -gff $genome &> RepeatMasker.log");
