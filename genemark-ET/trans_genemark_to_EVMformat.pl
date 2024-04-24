#!/usr/bin/perl -w
#Xu Cai
use strict;

my $in0 = $ARGV[0]; ##- genemark.gff3
my $out = $in0.".EVM.format";

my %gene = ();
my %gene2Len = ();
my $id;
open IN0, $in0;
while(<IN0>){
  next, if(/^\s+$/);
  if(/^[^#]/){
     chomp;
     my @temp = split(/\t/, $_);
     if($temp[2] eq "gene" && $temp[8] =~ /^ID=(\S+?);/){
        $id = $1;
        $gene{$id} = $_;
        $gene2Len{$id} = 0; 
     } 
     else{
        $gene{$id} .= "\n".$_;
        if($temp[2] eq "CDS"){
           $gene2Len{$id} += $temp[4] - $temp[3] + 1;
        }
     }
  }
}
close IN0;

open OUT0, ">$out";
for my $key(sort keys %gene){
    if($gene2Len{$key} >= 150){
       print OUT0 $gene{$key}, "\n";
    }
    else{
       print "filtered gene id: ", $key, "\n";
    }
}
close OUT0;
