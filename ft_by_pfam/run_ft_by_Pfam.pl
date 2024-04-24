#!/usr/bin/perl -w
##- Xu Cai
use strict;

my $in0 = $ARGV[0]; ##- gff3 file
my $in1 = $ARGV[1]; ##- Non_SynGenes.list
my $in2 = $ARGV[2]; ##- genome fasta
my $gene_prefix = $ARGV[3]; ##- gene_prefix
   $gene_prefix = 'gene', if(not defined $gene_prefix);

my $PfamValidateABinitio = '/40t_1/caix/software/tools/geta/PfamValidateABinitio';
my $GFF3Clear = '/40t_1/caix/software/tools/geta/GFF3Clear';
my $cpus = 30;

    &main();
sub main {

    ##- index gff3
    my %gene2info = ();
    my $cmd;
    &read_gff3($in0, \%gene2info);

    ##- remove Syn genes info 
    my $SynGenesGff3 = "SynGenes.gff3";
    my $Non_SynGeneGff3 = "NonSynGenes.gff3";
       &devide_Syn_Non($in1, \%gene2info, $SynGenesGff3, $Non_SynGeneGff3);    

    ##- ft gene by Pfam
    $cmd = "$PfamValidateABinitio --cpu $cpus $Non_SynGeneGff3 $in2 2> PfamValidateABinitio.log";
    system("$cmd") == 0 or die "failed to execute: $cmd\n";
    `rm -rf hmmscan.tmp`;
    
    ##- merge Syn genes and filtered results
    `cat $SynGenesGff3 out.filter_pass.gff3 > genome.merge.gff3`;

    ##- gff3clean
    $cmd = "perl $GFF3Clear --gene_prefix $gene_prefix --genome $in2 --no_attr_add genome.merge.gff3 >  genome.gff3";
    system("$cmd") == 0 or die "failed to execute: $cmd\n"; 

}


sub  devide_Syn_Non {

    my ($nonSynFile, $gene, $SynGenesGff3, $Non_SynGeneGff3) = @_;
    
    my %NonSynGenes = ();
    open IN1, $nonSynFile;
    while(<IN1>){
      $NonSynGenes{$1} = "Y", if(/^(\S+)/);   ##- it depends__  modified non_syntenic genes id  
    }
    close IN1; 

    ##- output
    open OUT0, ">$SynGenesGff3";
    open OUT1, ">$Non_SynGeneGff3";
    for my $key(sort keys %{$gene}){
        if(exists $NonSynGenes{$key}){
           print OUT1 $gene ->{$key}, "\n";
        }
        else{
           print OUT0 $gene ->{$key}, "\n";
        } 
    }
    close OUT0;
    close OUT1;

}


sub read_gff3 {

    my ($gff3, $gene) = @_;
   
    open IN0, $gff3;
    while(<IN0>){
      next, if(/^\s+$/ || /^#/);
      chomp;
      my @temp = split(/\t/, $_);
      if($temp[2] =~ /gene/ && $temp[8] =~ /^ID=(\S+?);/){     
         $gene ->{$1} = $_;
      }
      else{
         $gene ->{$1} .= "\n".$_;
      }
    }
    close IN0;

}
