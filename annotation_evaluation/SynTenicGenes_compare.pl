#!/usr/bin/perl -w
use strict;
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; # reference protein sequences. eg. Brassica_rapa.20100225.pep
my $in1 = $ARGV[1]; # variety protein sequences.   eg. Brapa_v3.0_pep.fasta
my $in2 = $ARGV[2]; # Ref cds sequences fasta   eg. V1.5 cds file
my $in3 = $ARGV[3]; # variety cds sequences fasta  eg. V3.0 cds file
my $in4 = $ARGV[4]; # Ref to Variety

    &main();
sub main {

   ##- get pep cds sequences
   my (%pepRef2seq, %pepVariety2seq, %Refgene2seq, %Varietygene2seq);
       &readFasta($in0, \%pepRef2seq);
       &readFasta($in1, \%pepVariety2seq);
       &readFasta($in2, \%Refgene2seq);
       &readFasta($in4, \%Varietygene2seq);

   ##- read SynOrth results
   my %Variety2Ref_Syn = ();
      &readSynOrth($in4, \%Variety2Ref_Syn);
   
   ##- caculate similarity
   my $out = "Variety_Ref.SynGene.similarity.list";
      &output(\%Variety2Ref_Syn, \%pepRef2seq, \%pepVariety2seq, \%Refgene2seq, \%Varietygene2seq, $out);

}

sub output {

    my ($Variety2Ref_Syn, $pepRef2seq, $pepVariety2seq, $Refgene2seq, $Varietygene2seq, $out) = @_;
     
    open OUT0, ">$out";
    print OUT0 join("\t", 'Varietygene', 'Refgene', 'AA_similarity(Variety_to_Ref)', 'Nt_similarity(Variety_to_Ref)'), "\n";
    for my $key(sort keys %{$Variety2Ref_Syn}){
        
        my @temp = ($key, $Variety2Ref_Syn ->{$key});
        open OUT3, ">Variety_to_Ref.fa";
        $pepVariety2seq->{$temp[0]} =~ s/\*//g;
        $pepRef2seq->{$temp[1]} =~ s/\*//g;
        print OUT3 ">$temp[0]\n$pepVariety2seq->{$temp[0]}\n>$temp[1]\n$pepRef2seq->{$temp[1]}\n";
        close OUT3;        
        my $simVariety2Ref = &similarity('Variety_to_Ref.fa'); 

        open OUT4, ">Variety_to_Ref.gene.fa";
        print OUT4 ">$temp[0]\n$Varietygene2seq->{$temp[0]}\n>$temp[1]\n$Refgene2seq->{$temp[1]}\n";
        close OUT4;
        my $simVariety2RefGene = &similarity('Variety_to_Ref.gene.fa');
        print OUT0 join("\t", $key, $simVariety2Ref, $simVariety2RefGene), "\n";
        system("rm Variety_to_Ref.fa Variety_to_Ref.gene.fa");
    }
    close OUT0;

}

sub similarity {
  my($twoseq) = @_;
  my $muscleF = "Temp.fa";

  system("muscle3.8.31_i86linux64 -in $twoseq -quiet  -clw  -out $muscleF");
  my($len, $matched, $same);
  open IN0, $muscleF;
  <IN0>;
  while(<IN0>){
    next, if(/^\s+$/);
    chomp;
    my @temp = split(/\s+/, $_);
       $len += length($temp[1]);
       <IN0>;
    my $matched_line = <IN0>;
       $same = $matched_line =~ tr/\*/\*/;
       $matched += $same;
  }
  close IN0;
  my $sim = $matched/$len;
  return($sim);
  system("rm Temp.fa");

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
