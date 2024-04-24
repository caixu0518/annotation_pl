#!/usr/bin/perl -w
##- Xu Cai
use strict;
use threads;
require("/home/caix/bin/all-sub.pl");

my $query = $ARGV[0]; # query file
my $subject = $ARGV[1]; # reference fasta genome
   
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
  my %gene2seq = ();  
     &readFasta($queryFile, \%gene2seq);  
  my @gene = keys %gene2seq; 
  my $parallelNum = 10;  ##- it depends__
  my $eachNum = int (scalar(@gene)/$parallelNum);
  my @fileName = ();

  ##------------split fasta sequence-----------------------------------
  for(my $n = 1; $n <= $parallelNum; $n++){
      my $file = $queryFile.".".($n - 1);
      push(@fileName, $file);
      if($n <= $parallelNum -1){
         open OUT0, ">$file";
         for(my $m = ($n -1)*$eachNum; $m < $n*$eachNum; $m++){
             print OUT0 ">", $gene[$m], "\n", $gene2seq{$gene[$m]}, "\n";
         }
         close OUT0;
      }
      else{
         ##- print last file
         open OUT1, ">$file";
         for(my $m = ($n -1)*$eachNum; $m <= $#gene; $m++){
             print OUT1 ">", $gene[$m], "\n", $gene2seq{$gene[$m]}, "\n";
         }
         close OUT1;
      }
  }
  ##---use threads------------------------------------------------------------------
  my @thr;
   for(my $i=0; $i< $parallelNum; $i+=1){
       $thr[$i] = threads->create(\&run_blastn, $fileName[$i], $refFile, "$i.blastn.out");
   }
   for(my $i=0; $i< $parallelNum; $i+=1){
       $thr[$i]->join;
   }
   
   `rm $blast_out`, if(-e $blast_out);
   for(my $i=0; $i< $parallelNum; $i+=1){
       `cat "$i.blastn.out" >> $blast_out`;
       `rm "$i.blastn.out"`;
       `rm $query.".".$i`;
   }
   ##-------------------------------------------------------------------------------
}

sub run_blastn {
  my ($query, $ref, $out) = @_;
  system("blastn -db $ref  -query  $query  -out $out  -evalue 1e-5 -num_threads 4  -task blastn  -outfmt 7 -perc_identity 100");
}

