#!/usr/bin/perl -w
##- Xu Cai
use strict;

my $in0 = $ARGV[0]; ##- genome gff3 file
my $in1 = $ARGV[1]; ##- genewise gff3 file
my $in2 = $ARGV[2]; ##- replaced.genes.list
#my $genome = $ARGV[3]; ##- genome fasta
#my $gene_prefix = $ARGV[4]; ##- gene prefix
#my $gff3clean = '/40t_1/caix/software/tools/geta/GFF3Clear';

    &main();
sub main {

    ##- read gff3
    my (%evmgenesIndex, %genewiseIndex);
        &read_gff3($in0, \%evmgenesIndex);
        &read_gff3($in1, \%genewiseIndex);
    
    ##- read replaced file
    my %replaced = ();
       &replaced($in2, \%replaced); 

    ##- output 
    my $out = $in0.".replaced";
       &output(\%evmgenesIndex, \%genewiseIndex, \%replaced, $out);    

    #my $cmdString = "perl $gff3clean --gene_prefix $gene_prefix  --genome $genome  $out > $out.clean";
    #system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

}

sub output {

    my ($evmgenesIndex, $genewiseIndex, $replaced, $out) = @_;

    my %gff = ();
    for my $key(sort keys %{$evmgenesIndex}){
        if(not exists $replaced ->{$key}){
           $gff{$key} = $evmgenesIndex ->{$key};
        }
        else{
           $gff{$key} = $genewiseIndex ->{$replaced ->{$key}};
        }
    }

    open OUT0, ">$out";
    for my $element(sort keys %gff){
        print OUT0 $gff{$element};
    }
    close OUT0;
   
}

sub replaced {
 
    my ($replacedFile, $replaced) = @_;
   
    open IN0, $replacedFile;
    while(<IN0>){
      chomp;
      my ($genewise, $evmgene) = (split(/\t/, $_))[0,1];
      if(not exists $replaced ->{$evmgene}){
         $replaced ->{$evmgene} = $genewise;
      }
      else{
         print "Waning: same gene id.\n";
      }
    }
    close IN0;

}

sub read_gff3 {
 
    my ($gff3File, $gff3index) = @_;

    open IN0, $gff3File;
    my $geneid;
    while(<IN0>){
       next if m/^#/;
       next if m/^\s*$/;
       if (/\tgene\t.*ID=([^;\s]+)/) {
           $geneid = $1;
           $gff3index ->{$geneid} = $_;

       }
       else{
           $gff3index ->{$geneid} .= $_;
       }
    }
    close IN0;

}

