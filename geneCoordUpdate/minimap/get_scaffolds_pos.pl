#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- aln.paf.sort
my $in1 = $ARGV[1]; ##= scaffolds sizes file

   &main();

sub main {
    
    ##- index scaffolds sizes
    my %scf2Sizes = ();
       &readSizes($in1, \%scf2Sizes);

    ##- read sorted paf file
    my $matchedCoords = "scfs_chr.coords.list";
       &read_paf(\%scf2Sizes, $in0, $matchedCoords);


}

sub read_paf {

    my ($index, $paf, $out) = @_;

    my $outTmp = $out.".Temp";
    open IN1, $paf;
    open OUT0, ">$outTmp";
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
      my ($s, $e) = ($temp[7], $temp[8]);
         ($s, $e) = ($temp[8], $temp[7]), if($temp[7] > $temp[8]);
         ($temp[7], $temp[8]) = ($s, $e);
      if((exists $index ->{$temp[0]}) && ($index ->{$temp[0]} == $temp[3] - $temp[2])  && $temp[9] == $temp[10] && $temp[11] == 60){
         my $tmp = $temp[8];
         if($temp[4] eq "-"){
            $temp[8] = $temp[7];
            $temp[7] = $tmp;
         }
         print OUT0 join("\t", $temp[0], $temp[1], $temp[2], $temp[3], $temp[4], $temp[5], $temp[6], $temp[7], $temp[8]), "\n";
      }
    }
    close IN1; 
    close OUT0;

    ##- filter multi-hits
    my %scfs2Chr = ();
    open IN1, $outTmp;
    while(<IN1>){
      chomp;
      my @temp =split(/\t/, $_);
      if(not exists $scfs2Chr{$temp[0]}){
         $scfs2Chr{$temp[0]} = $_;
      }
      else{
         $scfs2Chr{$temp[0]} = ";".$_;
      }
    }
    close IN1;

    my %notFind = ();
    for my $element(keys %{$index}){
        if(exists $scfs2Chr{$element}){
           if($scfs2Chr{$element} =~ /;/){
              print "One hit multi-locus: $element\n";
           }
        }
        else{
              print "cannot find hits: $element\n";
              $notFind{$element} = "Y";
        }
    }

    open IN, $paf;
    open OUT, ">fragmentsHits.list";
    while(<IN>){
      chomp;
      my @temp = split("\t", $_);
      my ($s, $e) = ($temp[7], $temp[8]);
         ($s, $e) = ($temp[8], $temp[7]), if($temp[7] > $temp[8]);
         ($temp[7], $temp[8]) = ($s, $e);
      if(exists $notFind{$temp[0]}){
         print OUT join("\t", $temp[0], $temp[1], $temp[2], $temp[3], $temp[4], $temp[5], $temp[6], $temp[7], $temp[8]), "\n";
      }
    }
    close IN;
    close OUT;

}


sub readSizes {

     my ($scafSizes, $index) = @_;

     open IN0, $scafSizes;
     while(<IN0>){
       chomp;
       my @temp = split(/\t/, $_);
       if(not exists $index ->{$temp[0]}){
          $index ->{$temp[0]} = $temp[1];
       }
       else{
          print "same scaffolds: $temp[0]\n";
       }
     }
     close IN0;

}
