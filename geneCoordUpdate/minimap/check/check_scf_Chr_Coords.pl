#!/usr/bin/perl -w
use strict;
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; ##- scaffold fasta
my $in1 = $ARGV[1]; ##- chromosome fasta
my $in2 = $ARGV[2]; ##- scf to chromosome coordinates

   my %scf2Seq = ();
   my %chr2seq = ();
      &readFasta($in0, \%scf2Seq);
      &readFasta($in1, \%chr2seq);

   my ($scafoldRegions, $chrRegions) = ("scf.regions.Lis", "chr.regions.Lis");
   my %index2info = ();   
      &read_paf($in2, $scafoldRegions, $chrRegions, \%index2info);

   ##- check cut length
   my ($scfSeq, $chrSeq) = ('scf.regions.Lis.fa', 'chr.regions.Lis.fa');
      &getSeqs($in0, $scafoldRegions, $scfSeq);
      &getSeqs($in1, $chrRegions, $chrSeq);

   my %scfSeqindex;
   my %chrSeqindex;
      &readFasta($scfSeq, \%scfSeqindex);
      &readFasta($chrSeq, \%chrSeqindex);
   
   ##- generate check info
      &check(\%scfSeqindex, \%chrSeqindex, \%index2info);

sub check {

    my ($scfSeqindex, $chrSeqindex, $index2info) = @_;
   
    for my $key(sort keys %{$scfSeqindex}){
        if(exists $chrSeqindex ->{$key} && exists $scfSeqindex ->{$key}){
           if($chrSeqindex ->{$key} eq $scfSeqindex ->{$key}){
              print join("\t", @{$index2info ->{$key}}), "\t", "Check OK!\n";
           }
           else{
              print join("\t", @{$index2info ->{$key}}), "\t", "Check Error!\n";
           }
        }
        else{
              die "cannot find index $key\n\n";
        }
    }

}


sub read_paf {

    my ($paf, $out1, $out2, $index2info) = @_;

    my $count = 0;
    open IN0, $paf;
    open OUT0, ">$out1";
    open OUT1, ">$out2";
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      $count += 1;
      my ($s1, $e1) = ($temp[7], $temp[8]);
         ($s1, $e1) = ($temp[8], $temp[7]), if($temp[7] > $temp[8]);
         ($temp[7], $temp[8]) = ($s1, $e1);
      $index2info ->{$count} = \@temp;
      print OUT0 join("\t", $temp[0], $temp[2]+1, $temp[3], "+", $count), "\n";      
      print OUT1 join("\t", $temp[5], $temp[7]+1, $temp[8], $temp[4], $count), "\n";
    }
    close IN0;
    close OUT0;
    close OUT1;

}

