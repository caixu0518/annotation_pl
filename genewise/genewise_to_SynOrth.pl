#!/usr/bin/perl -w
##- Xu Cai

my $in0 = $ARGV[0]; ##- genewise.gff3
my $genome = $ARGV[1]; ##- genome fasta
my $SynOrthResults = 'genewise_to_At'; ##- genewise_to_At
my $At_coord = '/40t_1/caix/Genomes/Genomes/Tair10/At10.gene.coords';
my $At_prot =  '/40t_1/caix/Genomes/Genomes/Tair10/TAIR10_GFF3_genes.gff.representative.pep';
`ln -s $At_coord At10.gene.coords .`, if(not -e "At10.gene.coords");
`ln -s $At_prot TAIR10_GFF3_genes.gff.representative.pep .`, if(not -e "TAIR10_GFF3_genes.gff.representative.pep");


my $gff3Filter = '/40t_1/caix/software/tools/geta/Filter_LQ.genes.pl';
my $gff3_to_protein = '/40t_1/caix/software/tools/geta/gff3_to_protein.pl';
my $gff3to_Coord = '/40t_1/caix/software/tools/geta/gff3ToSynOrthFromat.pl';
my $SynGeneSim = '/40t_1/caix/software/tools/geta/para_caculate_SynTenicGenes_Sim.pl';
my $cmdString = 'NA';

   ##- ft genewise.gff3 results
   my $genewise_ft = $in0.".HQ";
   
   $cmdString = "perl  $gff3Filter $in0  $genome $genewise_ft";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

   ##- clean gff3 to pep sequences
   $cmdString = "perl $gff3_to_protein  $genome $genewise_ft > proteins.fasta";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
   
   ##- trans gff3 to gene coordinates
   $cmdString = "perl $gff3to_Coord $genewise_ft";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n"; 
   
   ##- change pep format and make sure it matchs pep file
   $cmdString = "sed -i 's/\.mRNA//' proteins.fasta";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
   ##--run SynOrth
   
   $cmdString = "SynOrths -a proteins.fasta -b TAIR10_GFF3_genes.gff.representative.pep   -p $genewise_ft.sort.coords -q At10.gene.coords  -o $SynOrthResults";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
   system("rm -rf temp");

   ##- caculate_aa similarity
   $cmdString = "$SynGeneSim $At_prot  proteins.fasta $SynOrthResults";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
   
   print STDERR "All done.\n";   
