#!/usr/bin/perl -w
use strict;

my $in2 = $ARGV[0]; ##- genewise.gff3
my $in3 = $ARGV[1]; ##- final.evm.gff3.HQ
my $out = $ARGV[2]; ##-  

   &process();

sub process {
   
    my %gene2gff3 = ();
    my $gene = $in3.".gene.List";
       &read_gff3($in3, \%gene2gff3, $gene);

    my %genewise2gff3 = ();
    my $genewise = $in2.".gene.List";
       &read_gff3($in2, \%genewise2gff3, $genewise);

    #my %errorGenes = ();
    #   &index_errorGenes($in1, $in0, \%errorGenes);
    #my $num = scalar(keys %errorGenes);
    #print "Number of error genes: $num\n";
       
    ##- delete error genes
    #my $cleanGenes = $in3.".gene.List.clean";
    #   &del_Error($gene, \%errorGenes, $cleanGenes);

    ##- insert genewise results
    my $finalGff3 = $out;
    &insert($genewise, $gene, \%genewise2gff3, \%gene2gff3, $finalGff3);    

}


sub insert {
     
    my ($genewise, $cleangenes, $genewiseIndex, $geneIndex, $insertGff3) = @_;

    my %scfindex = ();
    my @insertGenes = ();
    open IN4, $genewise;
    while(<IN4>){
      chomp;
      my @temp = split(/\t/, $_);
      my ($start, $end) = ($temp[2], $temp[3]);
      if(not exists $scfindex{$temp[1]}){
         $scfindex{$temp[1]} = "Y";
      }

      if(not -e $temp[1]){
         `awk '{if(\$2~/$temp[1]/) print}' $cleangenes | sort -k2,2 -k3,3n > $temp[1]`;
      }
      ##- skip, if file is empty
      next, if(-z $temp[1]);

      my %scf2genes = ();
      my $num = 0;
      open IN1, $temp[1];
      while(my $line = <IN1>){
         chomp($line);
         $scf2genes{$num} = $line;
         $num += 1;         
      }
      close IN1;     
      $num = $num -1;      

      for(my $m = 0; $m< $num; $m++){
          my @first = split(/\t/, $scf2genes{$m});
          my @next  = split(/\t/, $scf2genes{$m+1});
          if($m == 0){
             if($end < $first[2]){
                push(@insertGenes, $temp[0]); 
                last;     
             }
          }
          if($m >0 && $m < $num){
             if($start > $first[3] && $end < $next[2]){
                push(@insertGenes, $temp[0]);
                last;
             }           
          }
      }
      my @lastLine = split(/\t/, $scf2genes{$num});
      if($start > $lastLine[3]){
         push(@insertGenes, $temp[0]);
      }
    }
    close IN4;
    
    open OUT, ">$insertGff3";
    open IN, $cleangenes;
    while(<IN>){
      chomp;
      my $geneID = (split(/\t/, $_))[0];
      if(exists $geneIndex ->{$geneID}){
         print OUT $geneIndex ->{$geneID}, "\n";
      }
      else{
         die "code 1.\n";
      }
    }
    close IN;
    
    print "Number of insert genes: ", scalar(@insertGenes), "\n";
    open OUT4, ">insertGenes.gff3";
    for my $element(@insertGenes){
        if(exists $genewiseIndex ->{$element}){
           print OUT $genewiseIndex ->{$element}, "\n";
           print OUT4 $genewiseIndex ->{$element}, "\n";
        }
    }
    close OUT;
    close OUT4;

    ##- clean 
    `rm $_`, for(keys %scfindex);
}



sub del_Error {

    my ($geneList, $errorGenes, $cleanGenes) = @_;

    open IN3, $geneList;
    open OUT3, ">$cleanGenes.Temp";
    while(<IN3>){
      chomp;
      my $gene = (split(/\t/, $_))[0];
      next, if(exists $errorGenes ->{$gene});
      print OUT3 $_, "\n";
    }
    close IN3;
    close OUT3;

    `sort -k2,2 -k3,3n $cleanGenes.Temp > $cleanGenes`;
    `rm $cleanGenes.Temp`;

}


sub index_errorGenes {

    my ($tandem, $error, $index) = @_;

    my %tandemGenes = ();
    open IN1, $tandem;
    while(<IN1>){
      chomp;
      my $info = (split(/\t/, $_))[1];
      my @temp = split(/,/, $info);
      my $representative = shift(@temp);      
         $tandemGenes{$representative} = \@temp;
    }
    close IN1;

    open IN2, $error;
    while(<IN2>){
      s/\.TA//g;
      chomp;
      my @temp = split(/;/, $_);
      for my $key(@temp){
          $index ->{$key} = "Y";
          if(exists $tandemGenes{$key}){
             my @info = split(/\t/, @{$tandemGenes{$key}});
             for my $tanGenes(@info){
                  $index ->{$tanGenes} = "Y";
             }
          }
      }
    }
    close IN2;

}


sub read_gff3 {

    my ($gff3, $index, $geneList) = @_;

    my $id;
    open IN0, $gff3;
    open OUT0, ">$geneList";
    while(<IN0>){
      next, if(/^\s+$/ ||  /^#/);
      chomp;
      if(/\tgene\t/ && (/ID=(\S+?);/ || /\tAUGUSTUS\t.*?ID=(\S+)$/)){
         my @temp = split(/\t/, $_);
         print OUT0 join("\t", $1, $temp[0], $temp[3], $temp[4], $temp[6]), "\n";
         $id = $1;
         $index ->{$id} = $_;
      }
      else{
         $index ->{$id} .= "\n".$_;
      }
    }
    close IN0;
    close OUT0;

} 
