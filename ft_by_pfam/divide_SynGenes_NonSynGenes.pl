#!/usr/bin/perl -w
##- Xu Cai
use strict;

my $gff3 = shift(@ARGV);
my $tandemList = shift(@ARGV);
my @SynFile = @ARGV;
   die "Please input SynOrth file ...\n\n", if(scalar @SynFile == 0);

    &main();
sub main {
   
   ##-- read gff3 
   my %geneList;
      &read_gff3($gff3, \%geneList);

   ##-- read tandem 
   my %tandem; 
      &read_tandem($tandemList, \%tandem);
  
   ##- read SynOrth file
   my $redundantSynGenes = "SynGenes.redundant";
      system("rm $redundantSynGenes"), if(-e $redundantSynGenes);
      for my $file(@SynFile){
          `cut -f 1 $file >> $redundantSynGenes`;
      }
   ##- output
   my $SynGenesList = "SynGenes.list";
   my $NonSynGeneList = "Non_SynGenes.list"; 
      &output(\%tandem, \%geneList, $redundantSynGenes, $SynGenesList, $NonSynGeneList); 
   `rm $redundantSynGenes`;

} 

sub output {

    my ($tandem, $gene, $redundantSynGenes, $SynGenesList, $NonSynGeneList) = @_;
 
    open OUT1, ">$SynGenesList";
    my %SynGenes = ();
    open IN2, $redundantSynGenes;
    while(<IN2>){
      if(/^(\S+)/ && not exists $SynGenes{$1}){
         $SynGenes{$1} = "Y";
         if(exists $tandem ->{$1}){
            my @TASynGenes = split(/,/, $tandem ->{$1});
            for my $key(@TASynGenes){
                $SynGenes{$key} = "Y";
                print OUT1 $key, "\n";
            }            
         }
         else{
            print OUT1 $1, "\n";
         }   
      }
      else{
         next;
      }
    }
    close OUT1;
    close IN2;

    ##- generate non_SynGenes
    open OUT2, ">$NonSynGeneList";
    for my $element(sort keys %{$gene}){
        print OUT2 $element, "\n", if(not exists $SynGenes{$element});
    }
    close OUT2;    

}


sub read_tandem {

     my ($tandemF, $tandem) = @_;
     open IN1, $tandemF;
     while(<IN1>) {
       chomp;
       my $TAgenes = (split(/\t/, $_))[1];
       my @TAgenes = split(/,/, $TAgenes);
       if(not exists $tandem ->{$TAgenes[0]}){
          $tandem ->{$TAgenes[0]} = $TAgenes;
       }
       else{
          print "Waning: same TA genes.\n";
       }
     }
     close IN1;


}

sub read_gff3 {
  
    my ($gff3, $geneList) = @_;
    
    open IN0, $gff3;
    while(<IN0>){
       chomp;
       my @temp = split(/\t/, $_);
       if($temp[2] =~ /gene/ && $temp[8] =~ /^ID=(\S+?);/){  ##- it depends_ mRNA or gene
          if(not exists $geneList ->{$1}){
             $geneList ->{$1} = "Y";
          }
          else{
             print "Waining: same gene id: $1\n";
          }
       }
    } 
    close IN0;

}
