#!/usr/bin/perl -w
use strict;

my $gff3 = $ARGV[0];   ##- gff3 file
my $prefix = $ARGV[1]; ##- prefix 
my $out = $gff3.'.rename.sort.gff3';

    &process();

sub process {

    my $geneList = 'gene.sort.list';
    my %indexGff3 = ();
    my %scfIndex = ();
       &read_gff3($gff3, $geneList, \%indexGff3, \%scfIndex);
     
    ##- rename
    #my $prefix = 'Bra';
    my $last =  $prefix; 
       &rename(\%scfIndex, $geneList, \%indexGff3, \$prefix, \$last, $out); 

}

sub rename {

    my ($scfIndex, $geneList, $indexGff3, $prefix, $last, $output) = @_;
 
    open OUT0, ">$output";
    for my $key(sort keys %{$scfIndex}){
        `awk '{if(\$2 == "$key") print}' $geneList  | sort -k2,2 -k3,3n > $key`, if(not -e $key);
         my $count = 0;
         open IN1, $key;
         while(<IN1>){
            chomp;
            my @temp = split(/\t/, $_);
            $count += 1;
            my $num = "";
            &formatNum($count, \$num);
            my $id = "$key.".$num.".$$last";
            if(exists $indexGff3 ->{$temp[0]}){
               my @info = split(/\n/, $indexGff3 ->{$temp[0]});
               for my $line(@info){
                   my @lineInfo = split(/\t/, $line);
                      $lineInfo[5] = '.';
                      if($lineInfo[2] eq 'gene'||  $lineInfo[2] eq 'mRNA' || $lineInfo[2] eq 'CDS'){ 
                         $lineInfo[8] = "ID=$id;", if($lineInfo[2] eq "gene");
                         $lineInfo[8] = "ID=$id.m1;Parent=$id;", if($lineInfo[2] eq "mRNA");
                         $lineInfo[8] = "ID=$id.m1.cds;Parent=$id.m1;", if($lineInfo[2] eq "CDS");
                         print OUT0 join("\t", @lineInfo), "\n";
                      }
                      else{
                         print "Wanning: $line will be ignore.\n";
                      }
               }
            }
            else{
               die "code 1.\n";
            }           

         }
         close IN1;
    }
    close OUT0;

    ##- clean
    `rm $_`, for(keys %{$scfIndex});

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

sub read_gff3 {

    my ($gff3, $geneList, $index, $scfIndex) = @_;
    
    my %coord2gene = ();
    my %indexgff3 = ();
    my $id;
    open IN0, $gff3;
    open OUT0, ">$geneList.Temp";
    while(<IN0>){
      next, if(/^\s+$/ ||  /^#/ || /\texon\t/ || /intron/ || /start_codon/ || /stop_codon/ || /five_prime_UTR/ || /three_prime_UTR/);
      s/\ttranscript\t/\tmRNA\t/;
      chomp;
      my @temp = split(/\t/, $_);
         $temp[1] = 'EVM';
      my $line = join("\t", @temp);
      if($temp[2] eq "gene" && ($temp[8] =~ /^ID=(\S+?);/ || ($temp[8] !~ /;/ && $temp[8] =~ /ID=(\S+)$/))){
         $id = $1;
         $indexgff3{$id} = $line;
         if(not exists $coord2gene{$temp[0]}{$temp[3]}{$temp[4]}){
            $coord2gene{$temp[0]}{$temp[3]}{$temp[4]} = $1;
            print OUT0 join("\t", $id, $temp[0], $temp[3], $temp[4], $temp[6]), "\n"; 
         }
         else{
            print "Two genes in the same locus: ", $1, "\t", $coord2gene{$temp[0]}{$temp[3]}{$temp[4]}, "\n";
         }

      }
      else{
         $indexgff3{$id} .= "\n".$line;
      }

      ##- index scaffolds
      $scfIndex ->{$temp[0]} = "Y", if(not exists $scfIndex ->{$temp[0]});

    }
    close IN0;
    close OUT0;
   
    `sort -k2,2 -k3,3n $geneList.Temp > $geneList`;
    `rm $geneList.Temp`;
 
     for my $key (sort keys %indexgff3){
         my $info = $indexgff3{$key};
         my $sortinfo = &sortgff($info);
         $index ->{$key} = $sortinfo;     
     }

}



sub sortgff {

        my $out;
        my $info = $_[0];
        my @line = split /\n/, $info;
        @_ = split /\t/, $line[0];
        my $strand = $_[6];
        my %sort;
        foreach (@line) {
                @_ = split /\t/;
                if ($_[2] eq "gene") {
                        $out .= "$_\n";
                }
                elsif ($_[2] eq "mRNA") {
                        $out .= "$_\n";
                }
                elsif ($_[2] eq "CDS" or $_[2] eq "exon" or $_[2] eq "five_prime_UTR" or $_[2] eq "three_prime_UTR") {
                        $sort{$_} = $_[3];
                }
        }

        if ($strand eq "+") {
                foreach (sort {$sort{$a} <=> $sort{$b} or $b cmp $a } keys %sort) {
                        $out .= "$_\n";
                }
        }
        elsif ($strand eq "-") {
                foreach (sort {$sort{$b} <=> $sort{$a} or $b cmp $a } keys %sort) {
                        $out .= "$_\n";
                }
        }

        return $out;
}
