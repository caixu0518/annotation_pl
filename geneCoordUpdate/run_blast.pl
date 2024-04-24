#!/usr/bin/perl -w
##- Xu Cai
use strict;
require("/home/caix/bin/all-sub.pl");

my $query = $ARGV[0]; # query file
my $subject = $ARGV[1]; # reference fasta genoem
   
   &main;

sub main {
 ##------Run blastn------------------------
  my $blastOut = "$query.blastn";
     &blastn($query, $subject, $blastOut);
 ##---------------------------------------- 
 
  my %gene2seq = ();
     &readFasta($query, \%gene2seq);
  
  my %matched_genes = ();
  my $output = "Matched_genes.list";
     &process_Blast($blastOut, \%gene2seq, \%matched_genes, $subject, $output);

  my $unmappedGenes = "Unmapped_genes.list";
  my $sameGenes = "One_hit_multi-locus.list";
     &printUnmappdGenes(\%gene2seq, $output, $unmappedGenes, $sameGenes);

}


sub printUnmappdGenes{
  my ($gene2seq, $output, $unmappedGenes, $sameGenes) = @_;
  
  my %mappedGenes =();
  open IN5, $output;
  while(<IN5>){
    chomp;
    my @temp = split("\t", $_);
    if(not exists $mappedGenes{$temp[0]}){
       $mappedGenes{$temp[0]} = $temp[1]."\t".$temp[2]."\t".$temp[3]."\t".$temp[4];
    }
    else{
       $mappedGenes{$temp[0]} .= ";".$temp[1]."\t".$temp[2]."\t".$temp[3]."\t".$temp[4];
    }
  }
  close IN5;
  
  open OUT3, ">$unmappedGenes";
  open OUT4, ">$sameGenes";
  for my $key1 (sort keys %{$gene2seq}){
    if(not exists $mappedGenes{$key1}){
       print OUT3 "Unmapped gene: $key1\n";
    }
    else{
       if($mappedGenes{$key1} =~ /;/){
          print OUT4 $key1, "\t", $mappedGenes{$key1}, "\n";
       }
    }
  }
  close OUT3;
  close OUT4;

}

sub process_Blast {
 my ($blast_out, $gene2seq, $matched_genes, $ref_genome, $output) = @_;
 my %chr2seq = ();
    &readFasta($ref_genome, \%chr2seq);  

 open IN0, $blast_out;
 open OUT0, ">$output.Temp";
 while(<IN0>){
   if(/^[^#]/){
     chomp;
     my @temp = split("\t", $_);
        if($temp[2] == 100 && $temp[3] == length($gene2seq ->{$temp[0]})){
           my ($seq, $stread, $s_start, $s_end);
           if($temp[8] < $temp[9]){
              $stread = "+";
              ($s_start, $s_end) = ($temp[8], $temp[9]);
           } 
           else{
              $stread = "-";
              ($s_start, $s_end) = ($temp[9], $temp[8]);
           }
           
           $seq = substr($chr2seq{$temp[1]}, $s_start - 1, $s_end - $s_start + 1);
           $seq =~ tr/ATCG/TAGC/ if($stread eq "-"); 
           $seq = reverse($seq) if($stread eq "-");

           if($gene2seq ->{$temp[0]} eq $seq){   
              print OUT0 join("\t", $temp[0], $temp[1], $temp[8], $temp[9], $stread), "\n";
           }
           else{
              print "Error: You need to check $temp[0]\n";
           }
        }
        else{
          next; # genes may be fragemented
        }      
     }
 } 
 close IN0;
 close OUT0;
 `sort -k2,2 -k3,3n $output.Temp > $output`;   
 `rm $output.Temp`;
}

sub blastn {
  my ($queryFile, $refFile, $blast_out) = @_;
  
  system("makeblastdb -in $refFile  -parse_seqids -hash_index  -out $refFile  -dbtype nucl");
  system("blastn -db $refFile  -query $queryFile  -out $blast_out  -evalue 1e-5 -num_threads 30   -outfmt 7 -perc_identity 100");

}

