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
	Usage:	$0 <scaffold.fa>  <pasa.gff> [-o <output_prefix>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput);
GetOptions ("o:s"=>\$oput);
die $usage unless (@ARGV>=2);
my ($iput, $pasa) = @ARGV;
$oput ||= 'cds';


my $oput_cds = "$oput.fa";
my $oput_lst = "$oput.lst";


my %gff;
ANNOTATION::Read_GFF ($pasa, \%gff);


my $i=0;
open my $OCDS, ">$oput_cds" or die;
open my $OLST, ">$oput_lst" or die;
open my $IPUT, "<$iput" or die;
while (!eof $IPUT) {
	my $scaffold = FASTA::read_fasta_one_seq ($IPUT);
	my $gff = $gff{$$scaffold{id}};
	$gff or next;
	ANNOTATION::extract_cds ($gff);
	foreach (@$gff) {
		my $start = $$_[0]{start}-1;
		my $end   = $$_[0]{end};
		my $len   = $end - $start;
		my %cds;
		$cds{id} = 'seq'.++$i;
		$cds{seq} = substr ($$scaffold{seq}, $start, $len);
		$cds{len} = $len;
		next if($cds{seq} =~ /NNNN/i); #fix by zhanghk
		FASTA::write_fasta_one_seq ($OCDS, \%cds);
		for (my $j=1;$j<=$#$_;$j++) {
			if ($$_[0]{strand} eq '+') {
				print $OLST join(" ", $cds{id},$$_[$j]{start}-$start,$$_[$j]{end}-$start), "\n";
			}
			else {
				my $k = $#$_ - $j + 1;
				print $OLST join(" ", $cds{id},$$_[$k]{end}-$start,$$_[$k]{start}-$start), "\n";
			}
		}
		print $OLST "\n";
	}
}
close $IPUT;
close $OCDS;
close $OLST;

__END__
