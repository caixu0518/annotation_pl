#!/usr/bin/perl -w
# Xu Cai
use strict;

my $in0 = $ARGV[0]; # Ref to At        eg. v15_to_At
my $in1 = $ARGV[1]; # variety to At    eg. v3.0_TAIR10
my $in2 = $ARGV[2]; # variety to Ref   eg. v3.0_v1.5
   
   &main($in0, $in1, $in2);

sub main{
  my($SynRef, $SynVariety, $Variety2Ref) = @_;
  
  my(%Refgene, %Varietygene, %Variety2Refgene); 
     &readSynOrth($SynRef, \%Refgene);
     &readSynOrth($SynVariety, \%Varietygene);
     &readSynOrth($Variety2Ref, \%Variety2Refgene);

##---------------------------------------------------------
  my %match = ();
     &hitAt(\%V15gene, \%V30gene, \%V30ToV15gene, \%match);  

  my $ReliableUniqGenes = "Ref_uniq.gene.list";
     &getReliableUniqGenes(\%V15gene, \%V30ToV15gene, $ReliableUniqGenes);

  open OUT0, ">matched.At.genes.list";
  for my $key1(sort keys %match){
    for my $key2(sort keys %{$match{$key1}}){
      print OUT0 join("\t", $key1, $key2, $match{$key1}{$key2}), "\n";
     }
  }
  close OUT0;
}

sub getReliableUniqGenes {
  my ($V15gene, $V30ToV15gene, $ReliableUniqGenes) = @_;

  open OUT1, ">$ReliableUniqGenes";
  
  my @V15matched = values %{$V30ToV15gene};
  my %matchedV15 = ();
     $matchedV15{$_} = "Y", for(@V15matched);
  
  for my $key1(sort keys %{$V15gene}){
    print OUT1 $key1, "\n", if(not exists $matchedV15{$key1});

  }
  close OUT1;

}


sub hitAt{
  my($V15gene, $V30gene, $V30ToV15gene, $match) = @_;
  
  for my $key(sort keys %{$V30ToV15gene}){
    if((exists $V30gene ->{$key}) && (exists $V15gene ->{($V30ToV15gene ->{$key})})){
      $match ->{$key} ->{$V30ToV15gene ->{$key}} = $V30gene ->{$key}, if(($V30gene ->{$key}) eq ($V15gene ->{($V30ToV15gene ->{$key})}));
    }
  }

}

sub readSynOrth {
  my ($SynOrth, $geneindex) = @_;
   
  open IN0, $SynOrth;
  while(<IN0>){
    chomp;
    my @temp = split("\t", $_);
    $geneindex ->{$temp[0]} = $temp[4];
  }
  close IN0;

}
