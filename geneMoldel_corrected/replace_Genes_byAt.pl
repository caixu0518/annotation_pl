#!/usr/bin/perl -w
##- Xu Cai
use strict;

my $in0 = $ARGV[0];  ##- SynOrth: genewise to At
my $in1 = $ARGV[1];  ##- SynOrth: evm gene to  At
my $in2 = $ARGV[2];  ##- genewise similarity
my $in3 = $ARGV[3];  ##- EVM genes similarity
my $in4 = $ARGV[4];  ##- overlap genes list


    &main();
sub main {

    ##- read synorth
    my (%at2genewise, %at2evmgenes);
        &read_SynOrth($in0, \%at2genewise);
        &read_SynOrth($in1, \%at2evmgenes);
    
    ##- index similarity
    my (%genewiseSim, %evmgeneSim);
        &read_sim($in2, \%genewiseSim);
        &read_sim($in3, \%evmgeneSim);
    
    ##- index overlap genes. ignore overlap genes
    my %overlapGenes = ();
        &index_overlapGenes($in4, \%overlapGenes);       

    ##- output 
    my $out = "replaced.genes.list";    
       &output(\%at2genewise, \%at2evmgenes, \%genewiseSim, \%evmgeneSim, \%overlapGenes, $out);    

}

sub output {
  
    my ($at2genewise, $at2evmgenes, $genewiseSim, $evmgeneSim, $overlapGenes, $out) = @_;

    open OUT0, ">$out";
    for my $key(sort keys %{$at2genewise}){
           ##- only coordinate overlap genes will be as candidate genes of replacing
           if(exists $at2evmgenes ->{$key}){
              my @genewise = split("\n", $at2genewise ->{$key});
              my @evm = split("\n", $at2evmgenes ->{$key});
              my ($flag, $results);
                 &get_pair(\@genewise, \@evm, \$flag, \$results);
                 if($flag > 0){
                    my @pairs = split(/;/, $results);
                    for my $eachPair(@pairs){
                        my ($geneid, $genewiseid);
                        if($eachPair =~ /(\S+?):(\S+)/){
                           ($geneid, $genewiseid) = ($1, $2);
                           if(exists $genewiseSim ->{$genewiseid} && exists $evmgeneSim ->{$geneid}){
                              if(($genewiseSim ->{$genewiseid})  > ($evmgeneSim ->{$geneid})){
                                 print OUT0 "$genewiseid\t$geneid\t", $genewiseSim ->{$genewiseid}, "\t", $evmgeneSim ->{$geneid}, "\n", if(not exists $overlapGenes ->{$genewiseid});
                              }
                           }
                           else{
                                 die "cannot find similarity index(genewise: $genewiseid; evm gene: $geneid)...\n";
                           }
                        }
                        else{
                            die "cannot find matched gene pairs.\n";
                        }

                    }

                 }
           }
    }
    close OUT0;

}

sub index_overlapGenes {

   my ($overlapList, $index) = @_;

   open IN, $overlapList;
   while(<IN>){
     chomp;
     my ($gene1, $gene2);
     if(/(\S+?)-(\S+)/){
        ($gene1, $gene2) = ($1, $2);
         $index ->{$gene1} = "Y";
         $index ->{$gene2} = "Y";
     }
   }
   close IN;

}

sub get_pair {

    my ($genewise, $evm, $flag, $results) = @_;

    $$flag = 0;
    my ($geneid, $chr, $start, $end);
    my ($geneid1, $chr1, $start1, $end1);
    for my $geneinfo(@{$evm}){
        ($geneid, $chr, $start, $end) = (split(/\t/, $geneinfo))[0,1,2,3];
        for my $genwiseinfo(@{$genewise}){
            if($genwiseinfo =~ /\t$chr\t/){
               ($geneid1, $chr1, $start1, $end1) = (split(/\t/, $genwiseinfo))[0,1,2,3];      
                if(($start1 >= $start && $start1 <= $end) || ($end1 >= $start && $end1 <= $end)){
                    ##- coordinate overlap genes
                    $$flag += 1;
                    if($$flag == 1){
                       $$results = $geneid.":".$geneid1;
                    }
                    elsif($$flag > 1){
                       $$results .= ";".$geneid.":".$geneid1;
                    }                    
                }
            }           
        }

    }

}


sub read_sim {

    my ($simFile, $indexsim) = @_;

    open IN1, $simFile;
    while(<IN1>){
      chomp;
      my @temp = split(/\t/, $_);
      if(not exists $indexsim ->{$temp[0]}){
         $indexsim ->{$temp[0]} = $temp[2];
      }
      else{
         die "same geneid: $temp[0]\n";
      }
    }
    close IN1;

}


sub read_SynOrth {

    my ($SynOrthFile, $at2other) = @_;
 
    open IN0, $SynOrthFile;
    while(<IN0>){
      chomp;
      my @temp = split(/\t/, $_);
      if(not exists $at2other ->{$temp[4]}){
         $at2other ->{$temp[4]} = $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3];
      }
      else{
         $at2other ->{$temp[4]} .= "\n".$temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3];
      }
    }
    close IN0;

}







