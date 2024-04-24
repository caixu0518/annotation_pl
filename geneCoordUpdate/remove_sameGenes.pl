#!/usr/bin/perl -w
use strict;

my $in0 = $ARGV[0]; ##- Matched_genes.list
my $out = $in0.".ft"; ##- Matched_genes.list.ft

my %coord2genes = ();
   open OUT, ">MultiGeneshitOneloci.log";
   open IN0, $in0;
   while(<IN0>){
     chomp;
     my @temp = split(/\t/, $_);
     my $tmp = $temp[2];
     if($temp[4] eq "-"){
        $temp[2] = $temp[3];
        $temp[3] = $tmp;
     } 
     my $id = join("\t", $temp[1], $temp[2], $temp[3]);
     if(not exists $coord2genes{$id}){
        $coord2genes{$id} = join("\t", @temp);
     }
     else{
        print OUT  "Multi-genes hit one loci by blast, remove redundant genes: $temp[0]\n";
        next;
     }
   }
   close IN0;
   close OUT;

my $out1 = $in0.".ft1";
   open OUT0, ">$out1";
   for my $key1(keys %coord2genes){
       print OUT0 $coord2genes{$key1}, "\n"; 
   }
   close OUT0;
   `sort -k2,2 -k3,3n $out1 > $out`;
   `rm $out1`;

##- detect overlap genes 
    my $oupt = "OverlapGenes.list";
    &detect_overlapGenes($out, $oupt);

##-----------------------------------------------------
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
                    my $overlap = 0;
                    my $info = "NA";
                       $overlap = ($e1 - $s2) + 1, if($s2 >= $s1 && $s2 <= $e1);             
                       $overlap = ($e2 - $s1) + 1, if($e2 >= $s1 && $e2 <= $e1);
                       $overlap = ($e1 - $s1) + 1, if($s2 <= $s1 && $e2 >= $e1);
                       $info = "$overlap:($id1:$s1-$e1;$id2:$s2-$e2)";
                       print OUT0 $id1, "-", $id2, "\t", $info, "\n";;
                }
                last, if($s2 > $e1);
            }
        }
    }
    close OUT0;

}
