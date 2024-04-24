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
	Usage:	$0 <scaffold.masked.fa> <snap.gff> [-o <output.gff>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput_gff);
GetOptions ("o:s"=>\$oput_gff);
die $usage unless (@ARGV >= 2);
my ($iput, @iput) = @ARGV;
$oput_gff ||= 'out.gff';


my %gff;
foreach (@iput) {
	ANNOTATION::Read_GFF ($_, \%gff);
}


my $count;
open my $OPUT, ">$oput_gff" or die;
open my $IPUT, "<$iput" or die;
while (!eof $IPUT) {
	my $scaffold = FASTA::read_fasta_one_seq ($IPUT);
	my $gff = $gff{$$scaffold{id}};
	$gff or next;
	motify ($_) foreach (@$gff);
	ANNOTATION::sort_gff ($gff);
	ANNOTATION::sort_gene ($gff);
	ANNOTATION::extract_fa ($gff, $scaffold);
#	ANNOTATION::filter_gap ($gff, 1000);
	ANNOTATION::unique_gene_id ($gff, \$count, $$scaffold{id}, 'SNAP');
	ANNOTATION::Print_GFF ($OPUT, $gff, $$scaffold{id});
}
close $OPUT;
close $IPUT;


sub motify {
	my $g = shift;
	$$_{ann} =~ s/Name=[^;]+\;?// foreach (@$g);
}

__END__
