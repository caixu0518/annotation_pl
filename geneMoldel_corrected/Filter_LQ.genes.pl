#!/usr/bin/perl -w
##- Xu Cai
use strict;
require("/home/caix/bin/all-sub.pl");

my $gff3 = $ARGV[0];   ##-  gff3
my $genome = $ARGV[1]; ##-  genome
my $out = $ARGV[2];    ##- HQ gff3

   &main();
sub main {
   
  my $pep = "first.pep.fa";
  my $cds = "first.cds.fa";
  `rm $pep $cds`, if(-e $cds && $pep);
  `sh /40t_1/caix/software/tools/bin/iTools.sh  $genome $gff3 first`;
  `gunzip first.cds.fa.gz first.pep.fa.gz`;

  my (%pepid2seq, %cdsid2seq) = @_;
      &readFasta($cds, \%cdsid2seq);
      &readFasta($pep, \%pepid2seq);

  my %fakegeneid = ();
  ##- detecte N in cds
  &cds_process(\%cdsid2seq, \%fakegeneid);
   
  ##- detect stop and short genes (< 50)
  &detect_fakeGenes(\%pepid2seq, \%fakegeneid);

  my $fakeList = "fakeGenes.list";
  &output(\%fakegeneid, $fakeList);

  ##- generate clean gff3
  my $cleangff3 = $out;
     &read_gff3($gff3, \%fakegeneid, $cleangff3);

}

sub read_gff3 {

   my ($gff3, $fakegenes, $cleangff3) = @_;

   my %gff3;
   my $id;
   open IN0, $gff3;
   while(<IN0>){
     next, if(/^\s+$/ || /^#/);
     chomp;
     my @temp = split(/\t/, $_);
     if($temp[2] eq "gene" && $temp[8] =~ /^ID=(\S+?);/){
        $id = $1;
        if(not exists $gff3{$id}){
           $gff3{$id} = $_;
        }
        else{
           print "$id\tgff3 file have multi-iso-from, please select representative genes first.\n";
        }
     }
     else{
        $gff3{$id} .= "\n".$_;
     }
   }
   close IN0;
   
   ##- format fake genes
   my %cleanFake = ();
   for my $key(sort keys %{$fakegenes}){
       my $fakegeneID;
       if($key =~ /(\S+?)\.m/){    ##- it depends_
          $fakegeneID = $1;
          if(not exists $cleanFake{$fakegeneID}){
             $cleanFake{$fakegeneID} = $fakegenes ->{$key};
          }
          else{
             print "gff3 file have multi-iso-from, please select representative genes first.\n";
          }
       }

   }

   open OUT1, ">$cleangff3";
   for my $key(sort keys %gff3){
       if(not exists $cleanFake{$key}){
          print OUT1 $gff3{$key}, "\n";
       }
   }
   close OUT1;

}


sub output {

   my ($fake, $out) = @_;

   open OUT0, ">$out";
   for my $key1(sort keys %{$fake}){
       print OUT0 $key1, "\n";
   }
   close OUT0;
 
}


sub detect_fakeGenes {

   my ($pep, $fake) = @_;

   for my $key(sort keys %{$pep}){
     my $seq = $pep ->{$key};
     my $len = length($seq);
     my $start=substr ($seq,0,1);
     my $end=substr ($seq,$len-1,1);
     if ($seq =~ /^(\w+)\*(\w+)/ || $start ne 'M' || $end ne '*' || $len < 50){
         $fake ->{$key} = "Y";
     }
   }

}

sub cds_process {

   my ($cds, $fake) = @_;

   for my $key(sort keys %{$cds}){
      if($cds ->{$key} =~ /N/){
         $fake ->{$key} = "Y";
         print "contain N genes: $key\n";
      }
   }

}
