#!/usr/bin/perl -w
use strict;
##- Xu Cai
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; ##- old_new_genes.list
my $in1 = $ARGV[1]; ##- evm gene pep file
my $in2 = $ARGV[2]; ##- update gene pep file


   my %scfSeqindex;
   my %chrSeqindex;
      &readFasta($in1, \%scfSeqindex);
      &readFasta($in2, \%chrSeqindex);

   ##- generate check info
   &check(\%scfSeqindex, \%chrSeqindex, $in0);
   

sub check {

    my ($scfSeqindex, $chrSeqindex, $matchedFile) = @_;

    open IN0, $matchedFile;
    while(<IN0>){
      chomp;
      my ($newid, $oldid) = (split(/\s+/, $_))[0,1];
      if(exists $chrSeqindex ->{$newid} && exists $scfSeqindex ->{$oldid}){
         if($chrSeqindex ->{$newid} eq $scfSeqindex ->{$oldid}){
            print "$newid-$oldid", "\t", "Check OK!\n";
         }
         else{
            print "$newid-$oldid", "\t", "Check Error!\n";
         }

      }
      else{
         die "cannot find gene id: $newid and $oldid\n\n";
      }

    }
    close IN0;

}










