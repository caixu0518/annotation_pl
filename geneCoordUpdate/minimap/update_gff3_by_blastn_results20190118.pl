#!/usr/bin/perl -w
use strict;

my $gff3 = $ARGV[0]; # old gff3 file
my $matchFile = $ARGV[1]; # matched file
my $prefixname = 'Bra';
my $prefixnameScf = 'BraAnng';
   
   &main($gff3, $matchFile);
   
sub main {
  my ($gff3F, $matchedF) = @_;
  
  ##- index Gff3 file  
  my %indexGff3 = ();
     &readGff3($gff3F, \%indexGff3);

  ##- process matched file and give new name
  my $old_newGenes = "old_new_genes.list";
     &readMatched($matchFile, $old_newGenes);

  my $output = "$gff3.update";
     &output(\%indexGff3, $old_newGenes, $output);

}


sub output {

    my ($gff3, $old_newGenes, $output) = @_;
    
    open IN, $old_newGenes;
    open OUT, ">$output";
    while(<IN>){
      chomp;
      my $updateCoords = 'NA';
      my $oldgeneid = (split(/\t/, $_))[1];
      if(exists $gff3 ->{$oldgeneid}){
         my ($oldgeneInfo, $matchLine) = ($gff3 ->{$oldgeneid}, $_);
            &updateGeneCoords(\$oldgeneInfo, \$matchLine, \$updateCoords);
      }
      else{
            die "cannot find geneid: $oldgeneid.\n\n";
      }
      print OUT $updateCoords, "\n";
    }
    close IN;
    close OUT;

}


sub updateGeneCoords {

    my ($oldCoords, $match, $newCoords) = @_;
    
    ##- read old gff3 
    my ($gene, $mrna);
    my @exon = ();
    my @cds = ();
    my @gff3 = split("\n", $$oldCoords);
    for my $each(@gff3){
        $gene = $each,      if($each =~ /\tgene\t/);
        $mrna = $each,      if($each =~ /\tmRNA\t/);
        push(@cds, $each),  if($each =~ /\tCDS\t/);
    }
    die "Fatal error: cannot initial CDS info in gff3 file.\n", if(scalar(@cds) == 0);

    ##- matched coordinate information
    my @oldgene = split(/\t/, $gene);
    my ($oldStart, $oldStop, $oldStrand) = ($oldgene[3],  $oldgene[4], $oldgene[6]);
    my ($newgeneid, $geneid, $chr, $start, $end, $strand) = (split(/\t/, $$match))[0,1,2,3,4,5];
       $gene = join("\t", $chr, "EVM", "gene", $start, $end, ".", $strand, ".", "ID=$newgeneid;");
       $mrna = join("\t", $chr, "EVM", "mRNA", $start, $end, ".", $strand, ".", "ID=$newgeneid.m1;Parent=$newgeneid;");

    my @cdsForm = ();
    my $cdsLen;
    my $line; 
    my ($newStart, $newEnd);
    
    ##- process CDS
    ##- new + and old strand +
    if($strand eq "+"  && $oldStrand eq "+"){
       for my $eachCDS(@cds){
           my @CDSinfo = split(/\t/, $eachCDS);
              $cdsLen = $CDSinfo[4] - $CDSinfo[3] + 1;
              $newStart = $CDSinfo[3] - $oldStart + $start;
              $newEnd = $newStart + $cdsLen -1;
              $line = join("\t", $chr, "EVM", "CDS", $newStart, $newEnd, ".", $strand, $CDSinfo[7], "ID=$newgeneid.m1.cds;Parent=$newgeneid.m1;");
              push(@cdsForm, $line);
       }
       die "Fatal error: code 1.\n", if(scalar(@cdsForm) == 0);         
    }

    ##- new + and old strand -
    if($strand eq "+"  && $oldStrand eq "-"){
       for my $eachCDS(@cds){
           my @CDSinfo = split(/\t/, $eachCDS);
              $cdsLen = $CDSinfo[4] - $CDSinfo[3] + 1;
              $newStart = $oldStop - $CDSinfo[4] + $start;
              $newEnd = $newStart + $cdsLen -1;
              $line = join("\t", $chr, "EVM", "CDS", $newStart, $newEnd, ".", $strand, $CDSinfo[7], "ID=$newgeneid.m1.cds;Parent=$newgeneid.m1;");
              push(@cdsForm, $line);
       }
       die "Fatal error: code 2.\n", if(scalar(@cdsForm) == 0);
    }  

    ##- new - and old strand +
    if($strand eq "-"  && $oldStrand eq "+"){
       for my $eachCDS(@cds){
           my @CDSinfo = split(/\t/, $eachCDS);
              $cdsLen = $CDSinfo[4] - $CDSinfo[3] + 1;
              $newStart = $oldStop - $CDSinfo[4] + $start;
              $newEnd = $newStart + $cdsLen -1;       
              $line = join("\t", $chr, "EVM", "CDS", $newStart, $newEnd, ".", $strand, $CDSinfo[7], "ID=$newgeneid.m1.cds;Parent=$newgeneid.m1;");
              push(@cdsForm, $line);
       }
       die "Fatal error: code 3.\n", if(scalar(@cdsForm) == 0);
    }
    
    ##- new - and old strand -
    if($strand eq "-"  && $oldStrand eq "-"){
       for my $eachCDS(@cds){
           my @CDSinfo = split(/\t/, $eachCDS);
              $cdsLen = $CDSinfo[4] - $CDSinfo[3] + 1;
              $newStart = $CDSinfo[3] - $oldStart + $start;
              $newEnd = $newStart + $cdsLen -1;
              $line = join("\t", $chr, "EVM", "CDS", $newStart, $newEnd, ".", $strand, $CDSinfo[7], "ID=$newgeneid.m1.cds;Parent=$newgeneid.m1;");
              push(@cdsForm, $line);
       }
       die "Fatal error: code 4.\n", if(scalar(@cdsForm) == 0);
    }
    
    ##- store new coordinate information
    die "Fatal Error: $gene have no CDS, please check.\n", if(scalar(@cdsForm) == 0);
    $$newCoords = join("\n", $gene, $mrna, @cdsForm);

}

sub readMatched {

   my ($matchedFile, $newGenes) = @_;

   my $matchedFileTemp = $matchedFile.".Temp";
   open IN, $matchedFile;
   open OUT, ">$matchedFileTemp";
   while(<IN>){
     chomp;
     my @temp = split(/\t/, $_);
     my ($start, $end) = ($temp[2], $temp[3]);
        ($start, $end) = ($temp[3], $temp[2]), if($temp[2] > $temp[3]);
        ($temp[2], $temp[3]) = ($start, $end);
        print OUT join("\t", @temp), "\n";
   }
   close IN;
   close OUT;

   my $sortmatchedFile = $matchedFile.".sort";
   system("sort -k2,2 -k3,3n $matchedFileTemp > $sortmatchedFile");
   `rm $matchedFileTemp`;   

   my $newid;
   my %chr = ();
   my ($scfGenesCount, $chrGeneCount, $formatCount) = (0, 0, 0);
   open IN3, $sortmatchedFile;
   open OUT3, ">$newGenes";
   while(<IN3>) {
     chomp;
     my @temp = split(/\t/, $_);
     if($temp[1] =~ /Scaffold/ || $temp[1] =~ /scaffold/){

        ##- process  genes on scaffolds
        $scfGenesCount += 1;
        &formatNum($scfGenesCount, \$formatCount);
        $newid = $prefixnameScf.$formatCount.'.3C';
     }
     else{

        ##- process genes on chrs
        if(not exists $chr{$temp[1]}){
           $chrGeneCount = 1;
           $chr{$temp[1]} = "Y";
        }
        else{
           $chrGeneCount += 1;
        }
        &formatNum($chrGeneCount, \$formatCount);
        $newid = $prefixname.$temp[1]."g".$formatCount.'.3C';
     }
     print OUT3 join("\t", $newid, @temp), "\n";

   }  
   close IN3;
   close OUT3;

}


sub formatNum {
  my ($n, $num) = @_;

  $n = $n*10;
  $$num = "0000".$n  if($n >= 10 && $n < 100);
  $$num = "000".$n  if($n >= 100 && $n < 1000);
  $$num = "00".$n  if($n >= 1000 && $n < 10000);
  $$num = "0".$n  if($n >= 10000 && $n < 100000);
  $$num = $n  if($n >= 100000);

}


sub readGff3 {

  my ($gff3F, $indexGff3) = @_;
  
  open IN0, $gff3F;
  my $id;
  while(<IN0>){
      next, if(/^\s+$/ || /^#/ || /\texon\t/);
      chomp;
      my @temp = split("\t", $_);
      if($temp[2] =~ /gene/ && $temp[8] =~ /^ID=(\S+?);/){
         $id = $1;
         $indexGff3 ->{$id} = $_;
      } 
      else{
         $indexGff3 ->{$id} .= "\n".$_; # joined by "\n" elements
      }    
  }
  close IN0;

}
