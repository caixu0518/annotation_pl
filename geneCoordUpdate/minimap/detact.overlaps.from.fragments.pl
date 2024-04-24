#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- fragments merged (Temp and fragments file)

   &read_frags($in0);


sub read_frags {

    my ($fragsFile) = @_;
  
    my $tmp = $fragsFile.".sort";
    system("sort -k6,6 -k8,8n $fragsFile > $tmp");

    my %chrGenes = ();
    my %info = ();
    my %index2info = ();
    my $num = 0;
    my $count = 0;
    open IN0, $tmp;
    open OUT0, ">$fragsFile.addIndex";
    while(<IN0>){
      chomp;
      $num += 1;
      my @temp = split("\t", $_);
         $index2info{$num} = \@temp;
      print OUT0 join("\t", $num, @temp),"\n";
      push(@temp, $num);
      my ($start, $end) = ($temp[7], $temp[8]);
         ($start, $end) = ($temp[8], $temp[7]), if($temp[7] > $temp[8]);
         ($temp[7], $temp[8]) = ($start, $end);
      
      ##- count gene number in chromosomes
      if(not exists $chrGenes{$temp[5]}){
         $chrGenes{$temp[5]} = 1;
         $count = 1; 
         $info{$temp[5]}{$count} = \@temp;
      }
      else{
         $chrGenes{$temp[5]} += 1;
         $count += 1;
         $info{$temp[5]}{$count} = \@temp;
      }
    }   
    close IN0;
    close OUT0;
 
    open OUT, ">Map_Overlap.List";
    for my $key1(sort keys %info){
        for my $key2(sort {$a<=>$b} keys %{$info{$key1}}){
            my ($s1, $e1, $f_index) = ($info{$key1}{$key2} ->[7], $info{$key1}{$key2} ->[8],  $info{$key1}{$key2} ->[9]);
            my $link = $f_index;
            for(my $n=$key2+1; $n<= $chrGenes{$key1}; $n++){
                my ($s2, $e2, $s_index) = ($info{$key1}{$n} ->[7], $info{$key1}{$n} ->[8],  $info{$key1}{$n} ->[9]);
                
                if(($e2 >=$s1 && $e2 <= $e1) || ($s2 >= $s1 && $e2 <= $e1) || ($s2 >= $s1 && $s2 <= $e1) || ($s2 <= $s1 && $e2 >= $e1)){
                    $link .= ";".$s_index;
                }
            }
            my @mapinfo = ();
            if($link ne $f_index){
               my @link = split(/;/, $link);
               for(@link){
                   my ($S1, $S2, $Len) = ($index2info{$_} ->[2], $index2info{$_} ->[3], $index2info{$_} ->[1]);
                   #print join("\t", $_, $S2-$S1, $Len), "\n";
                   if(($S2 - $S1) == $Len){
                       $_ = "$_(C)";
                   }
                   else{
                       $_ = "$_(F)";
                   }
                   push(@mapinfo, $_);
               }     
               print OUT  join(";", @mapinfo), "\n";
            } 
        }
    }
    close OUT0;

}
