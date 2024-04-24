#!/usr/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;

#use FindBin qw($Bin);
use lib "/40t_1/caix/software/tools/PerlLib/";
use FASTA;
use ANNOTATION;

my $usage=<<USAGE;
	Usage:	$0 <evm.gff> [-o <output.gff>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput_gff);
GetOptions ("o:s"=>\$oput_gff);
die $usage unless @ARGV;
my (@iput) = @ARGV;
$oput_gff ||= 'out.gff';


my %gff;
foreach (@iput) {
	ANNOTATION::Read_GFF ($_, \%gff);
}


my $count;
open my $OPUT, ">$oput_gff" or die;
foreach my $scaffold_id (ANNOTATION::sort_scaffold_id(\%gff)) {
	my $gff = $gff{$scaffold_id};
	motify ($_) foreach (@$gff);
	ANNOTATION::sort_gff ($gff);
	ANNOTATION::unique_gene_id ($gff, \$count, $scaffold_id, 'EVM');
	ANNOTATION::Print_GFF ($OPUT, $gff, $scaffold_id);
}
close $OPUT;


sub motify {
	my $g = shift;
	$$_{ann} =~ s/Name=[^;]+\;?// foreach (@$g);
}

__END__
