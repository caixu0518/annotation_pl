#!/usr/bin/perl -w
use strict;
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; ##- gene coordinates file       final.evm.gff3.HQ.sort.coords
my $in1 = $ARGV[1]; ##- update gene matched file    gene.matched.list
my $in2 = $ARGV[2]; ##- scaffold fasta
my $in3 = $ARGV[3]; ##- chromosome fasta

   my ($GeneCoordstmp, $matchedCoordsTmp) = ("$in0.Tmp", $in1.".Tmp");
   `awk '{print \$2"\t"\$3"\t"\$4"\t"\$5"\t"\$1}' $in0 | sort -k1,1 -k2,2n > $GeneCoordstmp`;
   `awk '{print \$2"\t"\$3"\t"\$4"\t"\$5"\t"\$1}' $in1 | sort -k1,1 -k2,2n > $matchedCoordsTmp`;

   my ($scfGenesSeq, $chrGeneSeq) = ($GeneCoordstmp.'.fa', $matchedCoordsTmp.'.fa');
      &getSeqs($in2, $GeneCoordstmp, $scfGenesSeq);
      &getSeqs($in3, $matchedCoordsTmp, $chrGeneSeq);

   my %scfSeqindex;
   my %chrSeqindex;
      &readFasta($scfGenesSeq, \%scfSeqindex); 
      &readFasta($chrGeneSeq, \%chrSeqindex);

   ##- generate check info
   &check(\%scfSeqindex, \%chrSeqindex);

   ##- remove tmp file
   `rm $GeneCoordstmp $matchedCoordsTmp $scfGenesSeq $chrGeneSeq`;  
 

sub check {

    my ($scfSeqindex, $chrSeqindex) = @_;

    open OUT0, ">check.gene.info";
    for my $key(sort keys %{$chrSeqindex}){
        if(exists $chrSeqindex ->{$key} && exists $scfSeqindex ->{$key}){
           if($chrSeqindex ->{$key} eq $scfSeqindex ->{$key}){
              print OUT0 $key, "\t", "Check OK!\n";
              #print OUT0 ">$key\n", $scfSeqindex ->{$key}, "\n", $chrSeqindex ->{$key}, "\n";
           }
           else{
              print OUT0 $key, "\t", "Check Error!\n";
           }
        }
        else{
           die "cannot find index $key\n\n";
        }
    }
    close OUT0;

}
