#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw(mesh);
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../lib/";
use FASTA;
use TOOLS;
use MAKEFILE;

my $usage=<<USAGE;
        Usage:  $0 [-s <sample.tab>] [-r <ref.tab>] [-t <transcript.tab>] [-fungi <1|0>] [-o <output_dir>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");
######################################################

my ($sample_tab,$ref_tab,$transcript_tab,$fungi,$oput_dir);
GetOptions (
	"s:s"=>\$sample_tab,
	"r:s"=>\$ref_tab,
	"t:s"=>\$transcript_tab,
	"fungi:i"=>\$fungi,
	"o:s"=>\$oput_dir);
die $usage if (!$sample_tab || !$ref_tab || !$transcript_tab);
$oput_dir ||="./";
$fungi ||="0";

my $config_ini  = "$oput_dir/protein.ini";
if ($fungi == 0){
	-e $config_ini or `cp  $Bin/protein.ini  $config_ini`;}
else{
	-e $config_ini or `cp  $Bin/protein_fungus.ini  $config_ini`;}

my $config  = TOOLS::read_ini ($config_ini);
$$config{BIN} = $Bin;

my ($key1,$sample) = TOOLS::read_table ($sample_tab);
my ($key2,$ref) = TOOLS::read_table ($ref_tab);
my ($key3,$transcript) = TOOLS::read_table($transcript_tab);
my $homolog_protein = $$ref[-1]{HOMOLOG_PROTEIN_PEP};

my $splitn1 = $$config{FASTA_SPLIT} || 100;
my $splitn2 = $$config{FASTA_SPLIT} || 100;

my @make;
my $j=0;

foreach (my $i=0;$i<@$sample;$i++){
	my $smp = $$sample[$i];
	my $tran = $$transcript[$i];
	die "ERROR: the file $sample_tab doesn't match $transcript_tab" if ($$smp{sample} ne $$tran{sample});
	
	my $smp_dir = "$oput_dir/$$smp{sample}";
	`mkdir -p $smp_dir/fa1`;
	`mkdir -p $smp_dir/fa2`;
	my (@fa1,@fa1_pre,@fa2,@fa2_pre);
	@fa1 = FASTA::split_fasta_by_length ($$smp{scaffold}, $splitn1, "$smp_dir/fa1/fa1");
	for (my $i=0;$i<@fa1;$i++) {
        	($fa1[$i], $fa1_pre[$i]) = (TOOLS::absolute_path ($fa1[$i]))[0,1];
	}	
	@fa2 = FASTA::split_fasta_by_length ($$smp{scaffold_marked}, $splitn2, "$smp_dir/fa2/fa2");
	for (my $i=0;$i<@fa2;$i++) {
        ($fa2[$i], $fa2_pre[$i]) = (TOOLS::absolute_path ($fa2[$i]))[0,1];
	}	
	
	$$config{SPECIES} = $$smp{species};
	
	my (@pro_gff,@gene_gff);
	#exonerate($$smp{sample},$$smp{scaffold},$ref,\@pro_gff);
	genemoma($$smp{sample},$$smp{scaffold},$ref,\@pro_gff, $$tran{bam});
	die "Please input ref homolog_protein for CUFFLINKS_PASA_FILTER " unless ($$ref[-1]{NAME} eq "REF" or @$ref==1);
	cufflinks_pasa($$smp{sample},$$smp{scaffold},$$tran{fa},$$ref[-1]{HOMOLOG_PROTEIN_PEP});
	augustus_train($$smp{sample},$$smp{scaffold});
	augustus_rna($$smp{sample},\@fa1,\@fa1_pre,$$tran{bam},\@gene_gff);
##	augustus_run($$smp{sample},\@fa2,\@fa2_pre);
		
	snap_train ($$smp{sample},$$smp{scaffold});
	snap_run ($$smp{sample},$$smp{scaffold_marked},\@gene_gff);

	genemarkes ($$smp{sample},$$tran{bam},$$smp{scaffold},$$smp{scaffold_marked},\@gene_gff);

	glimmerhmm_train ($$smp{sample},$$smp{scaffold});
	glimmerhmm_run ($$smp{sample},\@fa2,\@fa2_pre,$$smp{scaffold_marked},\@gene_gff);
	
	
	evm ($$smp{sample},$$smp{scaffold},\@pro_gff,\@gene_gff);
	evm_pasa ($$smp{sample},$$smp{scaffold},$$tran{fa});	
	
	gene_stat($$smp{sample},\@pro_gff,\@gene_gff);
}

MAKEFILE::write_makefile ($oput_dir, 'protein', \@make);

#############################################
sub exonerate {
	my ($sample,$iput,$ref,$pro_gff)=@_;
	my $M = 90;
		
	foreach (my $N=0;$N<@$ref;$N++){
		my @tmp = ();
                my $tmp = "";
		my $NAME = $$ref[$N]{NAME};
		next if ($NAME eq "REF");
		my $PEP = $$ref[$N]{HOMOLOG_PROTEIN_PEP};
		for (my $i=1;$i<=$M;$i++) {
                        $make[$j]{I} = "";
			$make[$j]{D} = "$sample/01exonerate/$N$NAME/tmp/$i";
			$make[$j]{O} = "out.$i";
			$make[$j]{C} = TOOLS::cmd ($config, 'EXONERATE_CMD1', $PEP, $iput, $i, $M, "$make[$j]{O}.tmp");
                        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
                        push @tmp, "$make[$j]{D}/$make[$j]{O}";
                        $tmp .= " tmp/$i/$make[$j]{O}";
                        $j++;
		}
		$make[$j]{I} = join ' ', @tmp;
                $make[$j]{D} = "$sample/01exonerate/$N$NAME";
                $make[$j]{O} = "$NAME.exonerate.gff";
                $make[$j]{C} = TOOLS::cmd ($config, 'EXONERATE_CMD2', $tmp, "$NAME.exonerate.out", "exonerate.$NAME", "$make[$j]{O}.tmp");
                $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
		push (@$pro_gff,"01exonerate/$N$NAME/$NAME.exonerate.gff");
                $j++;
	}
#	$make[$j]{I} = join " ",(map{"$sample/01exonerate/$_"}@$pro_gff;
#	$make[$j]{D} = "$sample";
#	$make[$j]{O} = "exonerate.stat.xls";
#	$make[$j]{C} = TOOLS::cmd ($config, 'EXONERATE_CMD3',(join ",",@$pro_gff),"$make[$j]{O}.tmp");
#	$j++;
}

sub genemoma{
	my ($sample, $iput, $ref, $pro_gff, $bam) = @_;
	foreach (my $N=0;$N<@$ref;$N++){
		my $NAME = $$ref[$N]{NAME};
		next if ($NAME eq "REF");
		my ($gff, $fasta) = split /,/, $$ref[$N]{HOMOLOG_PROTEIN_PEP};
		$make[$j]{I} = "";
		$make[$j]{D} = "$sample/01genemoma/$N$NAME";
		$make[$j]{O} = "$NAME.genemoma.gff";
		#$make[$j]{C} = TOOLS::cmd ($config, 'GENEMOMA_CMD1', $gff, $fasta, $iput, "$make[$j]{O}.tmp");
		$make[$j]{C} = TOOLS::cmd ($config, 'GENEMOMA_CMD2', $gff, $fasta, $iput,$bam, "$make[$j]{O}.tmp");
		$make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '7G', 'scr.q');
		push (@$pro_gff,"01genemoma/$N$NAME/$NAME.genemoma.gff");
		$j++;
	}
}

sub cufflinks_pasa {
	my ($sample,$iput,$cufflinks_fa,$homolog_protein)=@_;
		
        $make[$j]{I} = '';
        $make[$j]{D} = "$sample/03cufflinks_pasa";
        $make[$j]{O} = 'pasa.1.end';
        $make[$j]{C} = TOOLS::cmd ($config, 'CUFFLINKS_PASA', $iput, $cufflinks_fa,"$make[$j]{O}.tmp");
#        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '20G', 'scr.q');
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '20G', 'scr.q');
	$j++;
        $make[$j]{I} = "$sample/03cufflinks_pasa/pasa.1.end";
        $make[$j]{D} = "$sample/03cufflinks_pasa";
        $make[$j]{O} = 'pasa.2.end';
        $make[$j]{C} = TOOLS::cmd ($config, 'CUFFLINKS_PASA_FILTER',$homolog_protein,"$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '20G', 'scr.q');
        $j++;
} 

sub augustus_train {
	my ($sample,$iput)=@_;
        $make[$j]{I} = "$sample/03cufflinks_pasa/pasa.2.end";
        $make[$j]{D} = "$sample/04augustus/train";
        $make[$j]{O} = 'augustus.train.end';
        $make[$j]{C} = TOOLS::cmd ($config, 'AUGUSTUS_TRAIN', $iput, "../../03cufflinks_pasa/pasa.train.gff", "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'scr.q');
        $j++;
}

sub augustus_rna {
	my ($sample,$fa,$fa_pre,$bam,$gene_gff)=@_;

        $make[$j]{I} = "$sample/04augustus/train/augustus.train.end";
        $make[$j]{D} = "$sample/04augustus";
        $make[$j]{O} = 'hints.gff';
        $make[$j]{C} = TOOLS::cmd ($config, 'AUGUSTUS_CMD3',$bam, "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
        $j++;
        my @tmp = ();
        my $tmp = "";
        for (my $i=0;$i<@$fa;$i++) {
                $make[$j]{I} = "$sample/04augustus/hints.gff";
                $make[$j]{D} = "$sample/04augustus/tmp/$$fa_pre[$i]";
                $make[$j]{O} = "$$fa_pre[$i].gff";
                $make[$j]{C} = TOOLS::cmd ($config, 'AUGUSTUS_CMD4', $$fa[$i], '../../hints.gff', "$make[$j]{O}.tmp");
                $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
                push @tmp, "$make[$j]{D}/$make[$j]{O}";
                $tmp .= " tmp/$$fa_pre[$i]/$make[$j]{O}";
                $j++;
        }
        $make[$j]{I} = join ' ', @tmp;
        $make[$j]{D} = "$sample/04augustus";
        $make[$j]{O} = 'augustus.gff';
        $make[$j]{C} = TOOLS::cmd ($config, 'AUGUSTUS_CMD2', $tmp, "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
	push (@$gene_gff,"04augustus/augustus.gff");
        $j++;
}

sub augustus_run {
	my ($sample,$fa,$fa_pre,$gene_gff)=@_;
        my @tmp = ();
        my $tmp = "";
        for (my $i=0;$i<@$fa;$i++) {
                $make[$j]{I} = "$sample/04augustus/train/augustus.train.end";
                $make[$j]{D} = "$sample/04augustus/tmp/$$fa_pre[$i]";
                $make[$j]{O} = "$$fa_pre[$i].gff";
                $make[$j]{C} = TOOLS::cmd ($config, 'AUGUSTUS_CMD1', $$fa[$i], "$make[$j]{O}.tmp");
                $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
                push @tmp, "$make[$j]{D}/$make[$j]{O}";
                $tmp .= " tmp/$$fa_pre[$i]/$make[$j]{O}";
                $j++;
        }
        $make[$j]{I} = join ' ', @tmp;
        $make[$j]{D} = "$sample/04augustus";
        $make[$j]{O} = 'augustus.gff';
        $make[$j]{C} = TOOLS::cmd ($config, 'AUGUSTUS_CMD2', $tmp, "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
	push (@$gene_gff,"04augustus/augustus.gff");
        $j++;
}

sub snap_train {
	my ($sample,$iput)=@_;
        $make[$j]{I} = "$sample/03cufflinks_pasa/pasa.2.end";
        $make[$j]{D} = "$sample/05snap/train";
        $make[$j]{O} = "$$config{SPECIES}.hmm";
        $make[$j]{C} = TOOLS::cmd ($config, 'SNAP_TRAIN', $iput, "../../03cufflinks_pasa/pasa.train.gff", "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
        $j++;
}


sub snap_run {
	my ($sample,$iput,$gene_gff)=@_;
        $make[$j]{I} = "$sample/05snap/train/$$config{SPECIES}.hmm";
        $make[$j]{D} = "$sample/05snap";
        $make[$j]{O} = "snap.gff";
        $make[$j]{C} = TOOLS::cmd ($config, 'SNAP_CMD1', $iput, "./train/$$config{SPECIES}.hmm", "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
	push (@$gene_gff,"05snap/snap.gff");
        $j++;
}


sub genemarkes {
	my ($sample,$bam,$iput1,$iput2,$gene_gff)=@_;
        $make[$j]{I} = '';
        $make[$j]{D} = "$sample/06genemarkes";
        $make[$j]{O} = 'genemarkes.gff';
        $make[$j]{C} = TOOLS::cmd ($config, 'GENEMARKES_CMD',$bam, $iput1, $iput2, "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '20G', 'scr.q');
	push(@$gene_gff,"06genemarkes/genemarkes.gff");
        $j++;
}

sub glimmerhmm_train {
	my ($sample,$iput)=@_;
        $make[$j]{I} = "$sample/03cufflinks_pasa/pasa.2.end";
        $make[$j]{D} = "$sample/07glimmerhmm/train";
        $make[$j]{O} = 'glimmer.train.end';
        $make[$j]{C} = TOOLS::cmd ($config, 'GLIMMERHMM_TRAIN', $iput, "../../03cufflinks_pasa/pasa.train.gff", "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
        $j++;
}


sub glimmerhmm_run {
	my ($sample,$fa,$fa_pre,$iput,$gene_gff)=@_;
        my @tmp = ();
        my $tmp = "";
        for (my $i=0;$i<@$fa;$i++) {
                $make[$j]{I} = "$sample/07glimmerhmm/train/glimmer.train.end";
                $make[$j]{D} = "$sample/07glimmerhmm/tmp/$$fa_pre[$i]";
                $make[$j]{O} = "$$fa_pre[$i].gff";
                $make[$j]{C} = TOOLS::cmd ($config, 'GLIMMERHMM_CMD1', $$fa[$i], "../../train/new", "$make[$j]{O}.tmp");
                $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
                push @tmp, "$make[$j]{D}/$make[$j]{O}";
                $tmp .= " tmp/$$fa_pre[$i]/$make[$j]{O}";
                $j++;
        }
        $make[$j]{I} = join ' ', @tmp;
        $make[$j]{D} = "$sample/07glimmerhmm";
        $make[$j]{O} = 'glimmerhmm.gff';
        $make[$j]{C} = TOOLS::cmd ($config, 'GLIMMERHMM_CMD2', $iput, $tmp, "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'all.q');
	push(@$gene_gff,"07glimmerhmm/glimmerhmm.gff");
        $j++;
}


sub evm {
	my ($sample,$iput,$p_gff,$g_gff)=@_;
	my $pro_gff = join " ",(map{"../$_"}@$p_gff);
	my $tran_gff = '../03cufflinks_pasa/*assemblies.fasta.transdecoder.genome.gff3'; 
	my $gene_gff = join " ",(map{"../$_"}@$g_gff);
        $make[$j]{I} = join " ",(map{"$sample/$_"}(@$p_gff,@$g_gff),"03cufflinks_pasa/pasa.1.end");
	$make[$j]{D} = "$sample/08evm";
        $make[$j]{O} = "evm.gff";
        $make[$j]{C} = TOOLS::cmd ($config, 'EVM_CMD1', $iput, $pro_gff, $tran_gff, $gene_gff, "$make[$j]{O}.tmp");
        $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '20G', 'scr.q');
        $j++;
	
}


sub evm_pasa {
	my ($sample,$iput,$cufflinks_fa)=@_;
        $make[$j]{I} = "$sample/08evm/evm.gff";
        $make[$j]{D} = "$sample/09evm_pasa";
        $make[$j]{O} = "evm_pasa.end";
        $make[$j]{C} = TOOLS::cmd ($config, 'EVM_PASA', $iput, $sample, "../08evm/evm.gff",$cufflinks_fa, "$make[$j]{O}.tmp");
#       $make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G', 'scr.q -l h=compute-0-0');
	$make[$j]{Q} = TOOLS::cmd ($config, 'QSUB', '3G','scr.q');
        $j++;
}

sub gene_stat {
	my ($sample,$p_gff,$g_gff)=@_;
	$make[$j]{I} = join " ",(map{"$sample/$_"}(@$p_gff,@$g_gff,"03cufflinks_pasa/pasa.1.end"));
	$make[$j]{D} = "$sample";
	$make[$j]{O} = "stat";
	$make[$j]{C} = TOOLS::cmd ($config, 'STAT_CMD',(join ",",@$p_gff),(join ",",@$g_gff),"03cufflinks_pasa/*assemblies.fasta.transdecoder.genome.gff3",$sample,"$make[$j]{O}.tmp");
	$j++;
}
__END__
