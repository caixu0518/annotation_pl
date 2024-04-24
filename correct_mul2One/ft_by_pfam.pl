#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- gff3 file
my $in1 = $ARGV[1]; ##- genome fasta
my $cpus = 30;

my $PfamValidateABinitio = '/40t_1/caix/software/tools/geta/PfamValidateABinitio';
my $cmd;

##- ft gene by Pfam
    $cmd = "$PfamValidateABinitio --cpu $cpus $in0  $in1  2> PfamValidateABinitio.log";
    system("$cmd") == 0 or die "failed to execute: $cmd\n";
    `rm -rf hmmscan.tmp`;
 
