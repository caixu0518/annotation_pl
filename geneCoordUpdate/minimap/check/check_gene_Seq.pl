#!/usr/bin/perl -w
use strict;
##- Xu Cai
require("/home/caix/bin/all-sub.pl");

my $in0 = $ARGV[0]; ##- old_new_genes.list
my $in1 = $ARGV[1]; ##- final.evm.gff3.HQ.sort.coords
my $in2 = $ARGV[2]; ##- scaffolds fasta 
my $in3 = $ARGV[3]; ##- genome fasta

   my ($sweetheartOld,  $sweetheartNew) = ('scfGenes.regions', 'chrGenes.regions');
   `awk '{print \$2"\t"\$3"\t"\$4"\t"\$5"\t"\$1}' $in1 | sort -k1,1 -k2,2n > $sweetheartOld`;
   `awk '{print \$3"\t"\$4"\t"\$5"\t"\$6"\t"\$2}' $in0 | sort -k1,1 -k2,2n > $sweetheartNew`;

   my ($scfGenesSeq, $chrGeneSeq) = ($sweetheartOld.'.fa',  $sweetheartNew.'.fa');
       &getSeqs($in2, $sweetheartOld, $scfGenesSeq);
       &getSeqs($in3, $sweetheartNew, $chrGeneSeq);

   my %scfSeqindex;
   my %chrSeqindex;
      &readFasta($scfGenesSeq, \%scfSeqindex);
      &readFasta($chrGeneSeq, \%chrSeqindex);

   ##- generate check info
   &check(\%scfSeqindex, \%chrSeqindex);
   
   ##- remove tmp file
   `rm $sweetheartOld $sweetheartNew $scfGenesSeq $chrGeneSeq`;


sub check {

    my ($scfSeqindex, $chrSeqindex) = @_;

    for my $key(sort keys %{$chrSeqindex}){
        if(exists $chrSeqindex ->{$key} && exists $scfSeqindex ->{$key}){
           if($chrSeqindex ->{$key} eq $scfSeqindex ->{$key}){
            #  print $key, "\t", "Check OK!\n";
           }
           else{
              print $key, "\t", "Check Error!\n";
           }
        }
        else{
           die "cannot find index $key\n\n";
        }
    }

}










