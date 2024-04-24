#!/usr/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../lib/";
use FASTA;
use ANNOTATION;

my $usage=<<USAGE;
	Usage:	$0 <exonerate.gff> [-o <output.gff>] [-p <program>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput_gff, $program);
GetOptions ("o:s"=>\$oput_gff,
            "p:s"=>\$program);
die $usage unless @ARGV;
my (@iput) = @ARGV;
$oput_gff ||= 'out.gff';
$program ||= 'exonerate';


my %gff;
foreach (@iput) {
	ANNOTATION::Read_Exonerate ($_, \%gff, $program);
}


my $count;
open my $OPUT, ">$oput_gff" or die;
foreach my $scaffold_id (ANNOTATION::sort_scaffold_id(\%gff)) {
	my $gff = $gff{$scaffold_id};
	ANNOTATION::sort_gff  ($gff);
	ANNOTATION::sort_gene ($gff);
	ANNOTATION::unique_gene_id ($gff, \$count, $scaffold_id);
	ANNOTATION::Print_GFF ($OPUT, $gff, $scaffold_id);
}
close $OPUT;


__END__
