#!/usr/bin/perl -w
# Xu Cai
use strict;
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; # TAIR10_GFF3_genes.gff.representative.pep
my $in1 = $ARGV[1]; # reference protein sequences. eg. Brassica_rapa.20100225.pep
my $in2 = $ARGV[2]; # variety protein sequences.   eg. Brapa_v3.0_pep.fasta
my $in3 = $ARGV[3]; # core genes (two genomes mapped to At, three column: variety geneid, refrence geneid, Atid)

my $in4 = $ARGV[4]; # Ref cds sequences fasta   eg. V1.5 cds file
my $in5 = $ARGV[5]; # variety cds sequences fasta  eg. V3.0 cds file

   &main($in0, $in1, $in2, $in3, $in4, $in5);

sub main {
  my ($pepAt, $pepRef, $pepVariety, $matchedGenes, $Refgene, $Varietygene) = @_;

  my (%pepAt2seq, %pepRef2seq, %pepVariety2seq, %Refgene2seq, %Varietygene2seq);
      &readFasta($pepAt, \%pepAt2seq);
      &readFasta($pepRef, \%pepRef2seq);
      &readFasta($pepVariety, \%pepVariety2seq);
      &readFasta($Refgene, \%Refgene2seq);
      &readFasta($Varietygene, \%Varietygene2seq);

  my $output = "Variety_Ref_At_similarity.genes.list";
     &output(\%pepAt2seq, \%pepRef2seq, \%pepVariety2seq, \%Refgene2seq, \%Varietygene2seq,  $output, $matchedGenes);
  
}

sub output{
 my ($pepAt2seq, $pepRef2seq, $pepVariety2seq, $Refgene2seq, $Varietygene2seq, $output, $matchedGenes) = @_;
 my ($simRef, $simVariety, $AA_Variety_to_Ref, $Nt_Variety_to_Ref) = (0, 0, 0);
 
 open OUT0, ">$output";
 print OUT0 join("\t", 'Varietygene', 'Refgene', 'Atgene', 'sim(Variety_to_At)', 'sim(Ref_to_At)', 'AA_similarity(Variety_to_Ref)', 'Nt_similarity(Variety_to_Ref)'), "\n";
 open IN1, $matchedGenes;
 while(<IN1>){
   chomp;
   my @temp = split("\t", $_);
  
     open OUT1, ">Variety_At.fa";
     $pepVariety2seq->{$temp[0]} =~ s/\*//g; 
     $pepAt2seq->{$temp[2]} =~ s/\*//g;
     print OUT1 ">$temp[0]\n$pepVariety2seq->{$temp[0]}\n>$temp[2]\n$pepAt2seq->{$temp[2]}\n";     
     close OUT1;
     $simV30 = &similarity('Variety_At.fa');     

     open OUT2, ">Ref_At.fa";
     $pepRef2seq->{$temp[1]} =~ s/\*//g;
     $pepAt2seq->{$temp[2]} =~ s/\*//g;
     print OUT2 ">$temp[1]\n$pepRef2seq->{$temp[1]}\n>$temp[2]\n$pepAt2seq->{$temp[2]}\n";
     close OUT2; 
     $simV15 = &similarity('Ref_At.fa');
   
     open OUT3, ">Variety_to_Ref.fa";     
     $pepVariety2seq->{$temp[0]} =~ s/\*//g;
     $pepRef2seq->{$temp[1]} =~ s/\*//g;
     print OUT3 ">$temp[0]\n$pepVariety2seq->{$temp[0]}\n>$temp[1]\n$pepRef2seq->{$temp[1]}\n";
     close OUT3;
     $simV15ToV30 = &similarity('Variety_to_Ref.fa');     

     open OUT4, ">Variety_to_Ref.gene.fa";     
     print OUT4 ">$temp[0]\n$Varietygene2seq->{$temp[0]}\n>$temp[1]\n$Refgene2seq->{$temp[1]}\n";
     close OUT4;
     $simV15ToV30gene = &similarity('Variety_to_Ref.gene.fa');

     print OUT0 join("\t", @temp, $simV30, $simV15, $simV15ToV30, $simV15ToV30gene), "\n";
     system("rm -rf Variety_At.fa Ref_At.fa Variety_to_Ref.fa Variety_to_Ref.gene.fa");
 }
 close IN1;
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
