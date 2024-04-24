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
	Usage:	$0 <augustus.gff> [-o <output.gff>]
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
	ANNOTATION::sort_gff ($gff);
	motify ($_) foreach (@$gff);
	ANNOTATION::unique_gene_id ($gff, \$count, $scaffold_id);
	ANNOTATION::Print_GFF ($OPUT, $gff, $scaffold_id);
}
close $OPUT;


sub motify {
	my $g = shift;
	my @g;
	foreach (my $i=0;$i<@$g;$i++) {
		if ($$g[$i]{type} eq 'CDS' && $$g[$i-1]{type} ne 'exon') {
			my $j = @g;
			$g[$j]{program} = $$g[$i]{program};
			$g[$j]{type}    = 'exon';
			$g[$j]{start}   = $$g[$i]{start};
			$g[$j]{end}     = $$g[$i]{end};
			$g[$j]{score}   = '.';
			$g[$j]{strand}  = $$g[$i]{strand};
			$g[$j]{phase}   = '.';
			$g[$j]{parent}  = $$g[$i]{parent};
		}
		push @g, $$g[$i] if ($$g[$i]{type} ne 'intron');
	}
	$$_{ann} = '' foreach (@g);
	@$g = @g;
}

__END__
