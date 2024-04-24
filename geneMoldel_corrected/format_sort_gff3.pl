#!/usr/bin/perl -w
use strict;
##- Xu Cai

my $in0 = $ARGV[0]; ##- gff3 file
my $out = "$in0.form";

my %gff3 = ();
my %gene2cdsNum = ();
my $temp = $in0.".tmp";
my $id;
open IN0, $in0;
open Temp, ">$temp";
while(<IN0>){
  next, if(/^\s+$/ || /^#/);
  next, if(/\tintron\t/ || /\tstart_codon\t/ || /\tstop_codon\t/ || /three_prime_UTR/ || /five_prime_UTR/);
  chomp;
  s/\t(GeneWise)\t/\tEVM\t/, if(m/\tGeneWise\t/);
  my @temp = split(/\t/, $_);
  if(m/\tgene\t.*ID=([^;\s]+)/){
     $id = $1;
     $gff3{$id} = $_;
     $gene2cdsNum{$id} = 0;
     print Temp join("\t", $temp[0], $temp[3], $temp[4], $temp[6], $id), "\n";
  }
  else{
     $gff3{$id} .= "\n".$_;
     $gene2cdsNum{$id} += 1, if(/\tCDS\t/);
  }
}
close Temp;
close IN0;

my %geneid = ();
`sort -k1,1 -k2,2n $temp > $temp.sort`;
open IN1, "$temp.sort";
open OUT1, ">$out";
while(<IN1>){
   chomp;
   my ($strand, $geneid) = (split(/\t/, $_))[3, 4];
   ##- remove one to multi-genes
   if(not exists $geneid{$geneid}){
     $geneid{$geneid} = "Y"; 
   ##-----------------------------------------
     if(exists $gff3{$geneid}){
        my @geneinfo = split(/\n/, $gff3{$geneid});
    
        open OUT0, ">$in0.Temp.txt";
        for (my $n = 2; $n <= $#geneinfo; $n++){
             print OUT0 $geneinfo[$n], "\n";
        }
        close OUT0;
        `sort -k3,3r -k4,4n $in0.Temp.txt > $in0.Temp.txt.sort`;
        my @gene = split(/\t/, $geneinfo[0]);
        my @mrna = split(/\t/, $geneinfo[1]);     
           $gene[8] = "ID=$geneid;";
           $mrna[8] = "ID=$geneid.m1;Parent=$geneid;";
        print OUT1 join("\t", @gene), "\n", join("\t", @mrna), "\n";
      
        my ($count_cds, $count_exon) = (0, 0);
           ($count_cds, $count_exon) = ($gene2cdsNum{$geneid}, $gene2cdsNum{$geneid}), if($strand eq "-");
        open IN2, "$in0.Temp.txt.sort";
        while(<IN2>){
             my @element = split(/\t/, $_);
             if($strand eq "+"){
                if($element[2] eq "CDS"){
                   $count_cds += 1;
                   $element[8] = "ID=$geneid.m1.cds;Parent=$geneid.m1;";
                }
                if($element[2] eq "exon"){
                   $count_exon += 1;
                   $element[8] = "ID=$geneid.m1.exon$count_exon;Parent=$geneid.m1;";
                }
             }          
             if($strand eq "-"){
                if(exists $gene2cdsNum{$geneid}){
                   if($element[2] eq "CDS"){
                      $element[8] = "ID=$geneid.m1.cds;Parent=$geneid.m1;";
                      $count_cds = $count_cds -1;
                   }
                   if($element[2] eq "exon"){
                      $element[8] = "ID=$geneid.m1.exon$count_exon;Parent=$geneid.m1;";
                      $count_exon = $count_exon -1;
                   }
                }
                else{
                      die "Error: code 1.\n";
                }
             }
             print OUT1 join("\t", @element), "\n";
        }
        close IN2;
        `rm $in0.Temp.txt $in0.Temp.txt.sort`;
     }
     else{
       die "cannot find $geneid.\n";
     }
     ##--------------------------------------------------------------------
   }
   else{
      print "same geneid: $geneid\n";
   }
}
close IN1;
close OUT1;
`rm $temp $temp.sort`;
