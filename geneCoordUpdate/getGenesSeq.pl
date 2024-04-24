#!/usr/bin/perl -w
use strict;
require("/home/caix/bin/all-sub.pl");

my $gff = $ARGV[0]; # gff3 file
my $genome = $ARGV[1]; # genome fasta
my $filtergenes = $ARGV[2]; ##--Optional fake genes list. like LTR related genes 

   &main();
sub main {

  if(defined $filtergenes){
     my %indexLTR_genes = ();
        &read_LTR_genes($filtergenes, \%indexLTR_genes);
        &getNofakeGenesSeq($gff, $genome, \%indexLTR_genes);
  }  
  else{
     &getGenesSeq($gff, $genome);
  }

}

##--------------------------------------------------------------------------------------------------
sub getNofakeGenesSeq {
  my ($gff3F, $fasta, $indexLTR_genes) = @_;
  
  my $out = "Genes.info";
  open IN0, $gff3F;
  open OUT0, ">$out";

  while(<IN0>){
    if(not /^\s+$/){
       if(/^[^#]/){
        chomp;
        my @temp = split("\t", $_);
        if($temp[2] eq "gene" && $temp[-1] =~ /^ID=(\S+?);/){
          print OUT0 join("\t", $temp[0], $temp[3], $temp[4], $temp[6], $1),"\n", if(not exists $indexLTR_genes ->{$1});
        }
       }
    }   
  }
  close IN0;
  close OUT0;
 
  my $sort = $out.".sort"; 
  system("sort -k1,1 -k2,2n $out > $sort");
  my $geneseq = "$gff3F.genes.fa";
  &getSeqs($fasta, $sort, $geneseq);
  system("rm $out $sort");
}

sub getGenesSeq {
  my ($gff3F, $fasta) = @_;

  my $out = "Genes.info";
  open IN0, $gff3F;
  open OUT0, ">$out";

  while(<IN0>){
    if(not /^\s+$/){
       if(/^[^#]/){
        chomp;
        my @temp = split("\t", $_);
        if($temp[2] eq "gene" && $temp[-1] =~ /^ID=(\S+?);/){
          print OUT0 join("\t", $temp[0], $temp[3], $temp[4], $temp[6], $1),"\n";
        }
      }
    }
  }
  close IN0;
  close OUT0;

  my $sort = $out.".sort";
  system("sort -k1,1 -k2,2n $out > $sort");
  my $geneseq = "$gff3F.genes.fa";
  &getSeqs($fasta, $sort, $geneseq);
  system("rm $out $sort");
}


sub read_LTR_genes {
  my ($LTR_genes, $indexLTR_genes) = @_;
  
  open IN10, $LTR_genes;
  while(<IN10>){
    $indexLTR_genes ->{$1} = "Y",  if(/^(\S+)/ && not exists $indexLTR_genes ->{$1});
  }
  close IN10;

}
