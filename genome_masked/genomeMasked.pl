#!/usr/bin/perl -w
use strict;

my $genome = $ARGV[0];
`ln -sf $genome genome.fasta`;
my $pwd = `pwd`;
chomp $pwd;
$genome = "$pwd/genome.fasta";

my $scripts= '/40t_1/caix/software/tools/geta';
my $REPEATMASKER = '/disk02/caix/src/RepeatMasker/RepeatMasker';
my $RepeatModeler = '/disk02/caix/src/RepeatModeler-open-1.0.10/RepeatModeler';
my $BuildDatabase = '/disk02/caix/src/RepeatModeler-open-1.0.10/BuildDatabase';
my $cpu = 30;
my $RM_species = "arabidopsis";
my $RM_lib;
my $cmdString = 'NA';
my $pwd = "NA";


# Step 0: RepeatMasker and RepeatModeler
print STDERR "\n============================================\n";
print STDERR "Step 0: RepeatMasker and RepeatModeler\n";
mkdir "0.RepeatMasker" unless -e "0.RepeatMasker";
unless (-e "0.RepeatMasker.ok") {
    chdir "0.RepeatMasker";
    $pwd = `pwd`; print STDERR "PWD: $pwd";
    
    # 进行RepeatMasker分析
    mkdir "repeatMasker" unless -e "repeatMasker";
    $cmdString = "$REPEATMASKER -e ncbi -gff  -pa $cpu -species $RM_species -dir repeatMasker/ $genome &> repeatmasker.log";
    unless (-e "RepeatMasker.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "RepeatMasker.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 进行RepeatModeler分析
    mkdir "repeatModeler" unless -e "repeatModeler";
    chdir "repeatModeler";
    $pwd = `pwd`; print STDERR "PWD: $pwd";
    unless ($RM_lib) {
        $cmdString = "$BuildDatabase  -name species -engine ncbi $genome";
        unless (-e "BuildDatabase.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "BuildDatabase.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }
        $cmdString = "$RepeatModeler  -engine ncbi -pa $cpu -database species &> RepeatModeler.log";
        unless (-e "RepeatModeler.ok") {
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            open OUT, ">", "RepeatModeler.ok" or die $!; close OUT;
        }
        else {
            print STDERR "CMD(Skipped): $cmdString\n";
        }
        $cmdString = "$REPEATMASKER  -e ncbi -gff  -pa $cpu -lib RM_\*/\*.classified -dir ./ $genome &> repeatmasker.log";
    }
    else {
        $cmdString = "$REPEATMASKER  -e ncbi -gff   -pa $cpu -lib $RM_lib -dir ./ $genome &> repeatmasker.log";
    }
    unless (-e "RepeatMasker.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "RepeatMasker.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }
    chdir "../";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    # 合并RepeatMasker和RepeatModeler的结果
    $cmdString = "$scripts/merge_repeatMasker_out.pl repeatMasker/genome.fasta.out repeatModeler/genome.fasta.out > genome.repeat.stats";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    $cmdString = "$scripts/maskedByGff.pl genome.repeat.gff3 $genome > genome.masked.fasta";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    $cmdString = "$scripts/maskedByGff.pl --mask_type softmask genome.repeat.gff3 $genome > genome.softmask.fasta";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
    chdir "../";
    open OUT, ">", "0.RepeatMasker.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 0 for the file 0.RepeatMasker.ok exists\n";
}
