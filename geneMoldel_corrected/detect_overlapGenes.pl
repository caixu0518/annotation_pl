#!/usr/bin/perl -w
##- Xu Cai
use strict;

my $in0 = $ARGV[0]; ##- gene gff3 file
my $out = $ARGV[1]; ##- overlap.genes.list
   $out = 'overlap.genes.list', if(not defined $out);

    &main();
sub main {

    my %gff3index = ();
    my $geneLists = "genes.list";
       &read_gff3($in0, \%gff3index, $geneLists);
     
    ##- detect overlap genes
    &detect_overlapGenes($geneLists, $out);


}

sub detect_overlapGenes {
   
    my ($geneList, $output) = @_;
    `sort -k2,2 -k3,3n $geneList > $geneList.sort`;

    my %geneinfo = ();
    my %scfs = ();
    my $count = 0;
    open IN1, "$geneList.sort";
    while(<IN1>){
       chomp;
       my @temp = split(/\t/, $_);
       if(not exists $scfs{$temp[1]}){
          $scfs{$temp[1]} = 1;
          $count += 1 ;
          $geneinfo{$temp[1]}{$count} = $_;
       }
       else{
          $scfs{$temp[1]} += 1;
          $count += 1;
          $geneinfo{$temp[1]}{$count} = $_;
       }
    }
    close IN1;    
 
    open OUT0, ">$output";
    for my $chr(sort keys %geneinfo){
        for my $coords(sort {$a<=>$b} keys %{$geneinfo{$chr}}){
            my ($id1, $s1, $e1) = (split(/\t/, $geneinfo{$chr}{$coords}))[0,2,3];
            for(my $n=$coords+1; $n<= $scfs{$chr}; $n++){
                my ($id2, $s2, $e2) = (split(/\t/, $geneinfo{$chr}{$n}))[0,2,3];
                if(($s2 >= $s1 && $s2 <= $e1) || ($e2 >= $s1 && $e2 <= $e1) || ($s2 <= $s1 && $e2 >= $e1)){
                    print OUT0 $id1, "-", $id2, "\n";;
                }
                last, if($s2 > $e1);
            }
        }
    }  
    close OUT0;

}

sub read_gff3 {

    my ($gff3File, $gff3index, $geneLists) = @_;

    open IN0, $gff3File;
    open OUT0, ">$geneLists";
    my $geneid;
    while(<IN0>){
       next if m/^#/;
       next if m/^\s*$/;
       my @temp = split(/\t/, $_);
       if (/\tgene\t.*ID=([^;\s]+)/) {
           $geneid = $1;
           $gff3index ->{$geneid} = $_;
           print OUT0 join("\t", $geneid, $temp[0], $temp[3], $temp[4], $temp[6]), "\n";
       }
       else{
           $gff3index ->{$geneid} .= $_;
       }
    }
    close IN0;
    close OUT0;

}

