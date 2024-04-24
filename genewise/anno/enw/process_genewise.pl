#!/usr/bin/perl -w
use strict;

my $protein = $ARGV[0]; ##-  protein 
my $ref = $ARGV[1]; ##- reference 
my $pathbin = "/40t_1/caix/software/tools/geta";
my $cpu = 8;


my $cmdString1 = "perl  $pathbin/homolog_genewise --cpu $cpu  --coverage_ratio 0.4 --evalue 1e-9  $protein   $ref   &> homolog_genewise.log";
`$cmdString1`;
print "Step1 done.\n";

my $cmdString2 = "perl $pathbin/homolog_genewiseGFF2GFF3  --min_score 15 --gene_prefix genewise --filterMiddleStopCodon   --input_genewise_start_info genewise.start_info.txt --output_start_and_stop_hints_of_augustus genewise.start_stop_hints.gff  --genome $ref  genewise.gff > genewise.gff3 2> genewise.gene_id_with_stop_codon.txt";
`$cmdString2`;
print "Step2 done.\n";

