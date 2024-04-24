#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; # query tandem file
my $in1 = $ARGV[1]; # reference tandem file
my $in2 = $ARGV[2]; # SynOrth file


my (%Qtandem, %Rtandem);
   &Tandem($in0, \%Qtandem);
   &Tandem($in1, \%Rtandem);

my $out = $in2.".Mul2One";
   &SynOrth($in2, \%Qtandem, \%Rtandem, $out);
  
 
##- - - 
sub SynOrth {
  my($SynOrth, $Qtandem, $Rtandem, $out) = @_;
  
  my (%synGene);
  open IN1, $SynOrth;
  while(<IN1>){
    chomp;
    my @temp = split("\t", $_);
    if(not exists $synGene{$temp[4]}){
      $synGene{$temp[4]} = $temp[0];
    }
    else{
      $synGene{$temp[4]} .= ";".$temp[0];
    }
  }
  close IN1;
  print "Number of reference Syntenic genes:", scalar(keys %synGene), "\n"; 

  open OUT0, ">$out";
  my ($newKey, $one2one);
  for my $key(sort keys %synGene){
    if(exists $Rtandem ->{$key}){
       $newKey = $key.".TA";
    }
    else{
       $newKey = $key;
    }
    
    if($synGene{$key} =~ /;/){
      my @info = split(/;/, $synGene{$key});
      for my $element(@info){
        $element = $element.".TA", if(exists $Qtandem ->{$element});
      }
      print OUT0 join("\t", $newKey, scalar(@info), join(";", @info)), "\n";
    }
    else{
      $one2one += 1;
    }
  } 
  close OUT0;
  print "Number of one to one syntenic genes:", $one2one, "\n";
}


sub Tandem {
  my ($file, $tandemG) = @_;

  open IN0, $file;
  while(<IN0>){
    chomp;
    my @temp = split("\t", $_);
    if($temp[1] =~ /^(\S+?),/ && ! exists $tandemG ->{$1}){
      $tandemG ->{$1} = $temp[0];
    }
  }
  close IN0;

}
