#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- fragmentsHits.list.sort
   &main();


sub main {

    my $out = $in0.".merge";
       &read_frags($in0, $out);

}

sub read_frags {

    my ($fragsFile, $out) = @_;
  
    my %index = ();
    my %queryLen = ();
    my %refLen = ();
    
    open IN0, $fragsFile;
    while(<IN0>){
      chomp;
      my @temp = split("\t", $_);
         $queryLen{$temp[0]} = $temp[1];
         $refLen{$temp[5]} = $temp[6];
      my ($start, $end) = ($temp[7], $temp[8]);
         ($start, $end) = ($temp[8], $temp[7]), if($temp[7] > $temp[8]);
         ($temp[7], $temp[8]) = ($start, $end);

      if(not exists $index{$temp[0]}{$temp[5]}){
         $index{$temp[0]}{$temp[5]}{1} = $temp[2]."\t".$temp[3];
         $index{$temp[0]}{$temp[5]}{2} = $temp[7]."\t".$temp[8];
         $index{$temp[0]}{$temp[5]}{3} = $temp[4];
      }
      else{
         $index{$temp[0]}{$temp[5]}{1} .= "\t".$temp[2]."\t".$temp[3];
         $index{$temp[0]}{$temp[5]}{2} .= "\t".$temp[7]."\t".$temp[8];
      }
    }   
    close IN0;

    my %record = ();
    open OUT0, ">$out";
    for my $key1(sort keys %index){
        for my $key2(keys %{$index{$key1}}){
            if($key1 eq $key2){
               $record{$key1} = "Y";
               print join("\t", $key1, $queryLen{$key1}, 0, $queryLen{$key1}, "+", $key2, $refLen{$key2},  0, $refLen{$key2}), "\n";
            }
            print OUT0 join("\t", $key1, $queryLen{$key1}, $index{$key1}{$key2}{1}, $index{$key1}{$key2}{3}, $key2, $refLen{$key2}, $index{$key1}{$key2}{2}), "\n";
        }
    }
    close OUT0;
 
    open OUT3, ">cannot.hits.log";
    open OUT4, ">cannot.hits.log.list";
    for my $key1(sort keys %index){
        for my $key2(keys %{$index{$key1}}){
            my @query = split("\t", $index{$key1}{$key2}{1});
            my @sortQuery = sort {$a<=>$b} @query;
            my @ref = split("\t", $index{$key1}{$key2}{2});
            my @sortRef = sort {$a<=>$b} @ref;
            print OUT3 join("\t", $key1, $queryLen{$key1}, $index{$key1}{$key2}{1}, $index{$key1}{$key2}{3}, $key2, $refLen{$key2}, $index{$key1}{$key2}{2}), "\n", if(not exists $record{$key1});
            print OUT4 join("\t", $key1, $queryLen{$key1}, $sortQuery[0], $sortQuery[-1],  $index{$key1}{$key2}{3}, $key2, $refLen{$key2}, $sortRef[0], $sortRef[-1]), "\n", if(not exists $record{$key1});
        }
    }
    close OUT3;
    close OUT4;

}
