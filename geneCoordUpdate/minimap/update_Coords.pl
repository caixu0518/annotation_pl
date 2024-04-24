#!/usr/bin/perl -w
use strict;

my $geneCoords = $ARGV[0]; ##- genome.gff3.replaced.ft.form.sort.coords
my $scf2Chr = $ARGV[1];    ##- scfs_chr.coords.list.merge.addIndex.manual

    &main();
sub main {

    ##- index coordinates to genes
    my %Coord2Genes = ();
       &indexGenesCoord($geneCoords, \%Coord2Genes);

    ##- index scfCoords
    my %scfCoords = ();
       &index_scfCoord($scf2Chr, \%scfCoords);

    ##- update coords  
    my $out = "gene.matched.list"; 
       &updateCoord(\%Coord2Genes, \%scfCoords, $out);

}

sub updateCoord {

    my ($Coord2Genes, $scfCoords, $out) = @_;

    open OUT, ">$out";
    for my $key1(sort keys %{$Coord2Genes}){
        if(exists $scfCoords ->{$key1}){
           for my $key2(sort {$a<=>$b} keys %{$Coord2Genes ->{$key1}}){
                  for my $key3(keys %{$Coord2Genes ->{$key1} ->{$key2}}){
                      for my $key4(keys %{$Coord2Genes ->{$key1} ->{$key2} ->{$key3}}){
                          ##- a whole scaffold
                          my ($oldGenes, $matched, $newGene);
                          $oldGenes = join("\t", $Coord2Genes ->{$key1} ->{$key2} ->{$key3} ->{$key4}, $key1, $key2, $key3, $key4);
                          if(not $scfCoords ->{$key1} =~ /;/){
                             $matched = $scfCoords ->{$key1};
                          }
                          else{
                             my @multi = split(/;/, $scfCoords ->{$key1});
                             for my $element(@multi){
                                 my @element = split("\t", $element);
                                 my ($s, $e) = ($element[2]+1, $element[3]+1);
                                 if($key2 >= $s && $key2 <= $e && $key3 >= $s && $key3 <= $e){
                                    $matched = $element;
                                 }
                                 else{
                                    next;
                                 }
                             }
                          }
                          if(defined $oldGenes && defined $matched){
                             &update(\$oldGenes, \$matched, \$newGene);
                             print OUT $newGene, "\n";  
                          }
                          else{
                             print "Warning: cannot find ",  $Coord2Genes ->{$key1} ->{$key2} ->{$key3} ->{$key4}, " in matched list\n\n";
                          }
                         
                      }
                  }
           }
        }
        else{
           die "cannot find scaffolds: $key1\n";
        }

    }
    close OUT;

}

sub update {
   
    my ($old, $match, $new) = @_;
    
    my @old_array = split(/\s+/, $$old);
    my @match_array = split(/\s+/, $$match);
    print $$match, "\n", if(not defined $match_array[2] || not defined $match_array[3] || not defined $match_array[7] || not defined $match_array[8]);  
 
    ##  gene00001       jcf7180000011159_1      17638   19633   +
    ##  jcf7180000011159_1      799999  0       799999  +       C03     76625942        53177310        53977309
    my ($s1, $e1, $chr, $scfStrand, $newS, $newE)   = ($match_array[2]+1, $match_array[3], $match_array[5], $match_array[4], $match_array[7]+1, $match_array[8]);
    my ($geneid, $geneS, $geneE, $geneScf, $strand) = ($old_array[0], $old_array[2], $old_array[3], $old_array[1], $old_array[4]);
    if($geneScf eq $match_array[0]){
       my $geneLen = $geneE - $geneS + 1; 
       my ($newStart, $newEnd, $NewStrand);

       ##- gene + and scaffold +
       if($scfStrand eq "+" && $strand eq "+"){
          $newStart = $geneS - $s1 + $newS;
          $NewStrand = "+";
       }

       ##- gene - and scaffold +
       if($scfStrand eq "+" && $strand eq "-"){
          $newStart = $geneS - $s1 + $newS;
          $NewStrand = "-";
       }
       
       ##- gene + and scaffold -
       if($scfStrand eq "-" && $strand eq "+"){
          $newStart = $e1 - $geneE + $newS;
          $NewStrand = "-";
       }       

       ##- gene - and scaffold -
       if($scfStrand eq "-" && $strand eq "-"){
          $newStart = $e1 - $geneE + $newS;
          $NewStrand = "+";
       }
       $newEnd = $newStart + $geneLen -1;
       $$new = join("\t", $geneid, $chr, $newStart, $newEnd, $NewStrand);
    } 
    else{
       die "cannot match scaffolds: $geneScf and $match_array[0]\n";
    }   

}


sub index_scfCoord {

    my ($scfCoord, $index) = @_;

    open IN1, $scfCoord;
    while(<IN1>){
      chomp;
      my @temp = split(/\s+/, $_);
      my ($s, $e) = ($temp[7], $temp[8]);
         ($s, $e) = ($temp[8], $temp[7]), if($temp[7] > $temp[8]);
         ($temp[7], $temp[8]) = ($s, $e);
      my $line = join("\t", @temp);    
 
      if(not exists $index ->{$temp[0]}){
         $index ->{$temp[0]} = $line;
      }
      else{
         $index ->{$temp[0]} .= ";".$line;
      }
    }
    close IN1;

}


sub indexGenesCoord {

    my ($geneCoord, $index) = @_;

    open IN0, $geneCoord;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      if(not exists $index ->{$temp[1]} ->{$temp[2]} ->{$temp[3]} ->{$temp[4]}){
         $index ->{$temp[1]} ->{$temp[2]} ->{$temp[3]} ->{$temp[4]} = $temp[0];
      }
      else{
         die "Fatal error: one loci has two genes, please check gene: $temp[0]\n";
      }
    }
    close IN0;
    
}
