#!/usr/bin/perl -w
use strict;
require("/home/caix/bin/all-sub.pl");
##- Xu Cai

my $gff3 = $ARGV[0];      ##-  gff3
my $genome = $ARGV[1];    ##-  genome
my $cleanGff3 = $gff3.".HQ"; ##-  clean gff3

    &process();
sub process {

    my $pep = "first.pep.fa";
    `rm $pep`, if(-e $pep);
    `sh /40t_1/caix/software/tools/bin/iTools.sh  $genome $gff3 first`;
    `gunzip first.pep.fa.gz`;

    ##- index pep
    my %id2seq = ();
       &readFasta($pep, \%id2seq);
   
    ##- detect stop and short genes (< 50)
    my %fakegeneid = ();
    my $fakeGenes = "fakegenes.list";
       &detect_fakeGenes(\%id2seq, \%fakegeneid, $fakeGenes);

    `rm first.cds.fa.gz`;

    ##- generate clean gff3
       &read_gff3($gff3, \%fakegeneid, $cleanGff3);

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
      if($temp[2] eq "gene" && ($temp[8] =~ /^ID=(\S+?);/ || ($temp[8] !~ /;/ && $temp[8] =~ /^ID=(\S+)$/))){
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
    print "Number of all genes: ", scalar(keys %gff3), "\n";

    ##- format fake genes
    my %cleanFake = ();
    for my $key(sort keys %{$fakegenes}){
        my $fakegeneID;
        if($key =~ /(\S+?)\.m/){   ##- it depends_
           $fakegeneID = $1;
           $cleanFake{$fakegeneID} = $fakegenes ->{$key};
        }
        else{
           die 'error 1', "\n";
        }
    }

    ##- ouput
    my $count = 0;
    open OUT1, ">$cleangff3";
    for my $key(sort keys %gff3){
       if(not exists $cleanFake{$key}){
          $count += 1;
          print OUT1 $gff3{$key}, "\n";
       }
    }
    close OUT1;
    print "Number of high quality genes: ", $count, "\n";
    print "Number of low quality genes: ", scalar(keys %{$fakegenes}), "\n";

}




sub detect_fakeGenes {

    my ($pep, $fake, $fakeGenesList) = @_;
    
    open OUT0, ">$fakeGenesList";
    for my $key(sort keys %{$pep}){
        my $seq   = $pep ->{$key};
        my $len   = length($seq);
        my $start = substr ($seq,0,1);
        my $end   = substr ($seq,$len-1,1);
        if($start ne 'M' || $seq =~ /(\w+)\*(\w+)/ || $end ne '*'  || $len < 50){
           $fake ->{$key} = "Y";
           print OUT0 $key, "\t", "not M start\n",     if($start ne 'M');
           print OUT0 $key, "\t", "not * end\n",       if($end ne '*');
           print OUT0 $key, "\t", "middle *\n",        if($seq =~ /(\w+)\*(\w+)/);
           print OUT0 $key, "\t", "length short 50\n", if($len < 50);
        }        
    }
    close OUT0;

}






