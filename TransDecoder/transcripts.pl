#!/usr/bin/perl -w
use strict;

my $sortedBam = $ARGV[0]; ##- star sorted bam 
my $genome = $ARGV[1]; ##- masked genome
my $cpu = 20;
my $binPath = '/40t_1/caix/software/tools/geta';
my $TransDecoder_Predict = '/40t_1/caix/software/gene_prediction/TransDecoder-2.0.1/TransDecoder.Predict';
my $TransDecoder_LongOrfs = '/40t_1/caix/software/gene_prediction/TransDecoder-2.0.1/TransDecoder.LongOrfs';
my $strand_specific;

my ($cmdString, $pwd);
   $cmdString = "samtools view -h $sortedBam > star.sorted.sam";
   print STDERR (localtime) . ": CMD: $cmdString\n";
   system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

##- Step 3: Transcript
print STDERR "\n============================================\n";
print STDERR "Step 3: Transcript\n";
mkdir "3.transcript" unless -e "3.transcript";
unless (-e "3.transcript.ok") {
    chdir "3.transcript";
    $pwd = `pwd`; print STDERR "PWD: $pwd";

    # 将SAM比对将诶过分成一个个比对区域
    $cmdString = "$binPath/split_sam_from_non_aligned_region ../star.sorted.sam splited_sam_out 10 > splited_sam_files.list";
    unless (-e "split_sam_from_non_aligned_region.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "split_sam_from_non_aligned_region.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 每个比对区域生成一个sam2transfrag命令
    open IN, "splited_sam_files.list" or die $!;
    open CMD, ">", "command.sam2transfrag.list" or die $!;
    my $no_strand_specific = "--no_strand_specific";
    $no_strand_specific = "" if $strand_specific;
    while (<IN>) {
        s/\.sam\n//;
		my $sam2transfrag_opt = '--fraction 0.05 --min_expressed_base_depth 2 --max_expressed_base_depth 50 --min_junction_depth 2 --max_junction_depth 50 --min_fragment_count_per_transfrags 10 --min_intron_length 20';
        print CMD "$binPath/sam2transfrag $sam2transfrag_opt $no_strand_specific --intron_info_out $_.intron $_.sam > $_.gtf\n";
    }
    close CMD;
    close IN;

    # 批量并行化进行transcripts计算
    $cmdString = "ParaFly -c command.sam2transfrag.list -CPU $cpu &> /dev/null";
    print STDERR (localtime) . ": CMD: $cmdString\n";
    system("$cmdString") == 0 or die "failed to execute: $cmdString\n";

    # 整合所有并行化结果，得到GTF文件和可信的Intron信息
    open IN, "splited_sam_files.list" or die $!;
    open OUT1, ">", "transfrag.gtf" or die $!;
    open OUT2, ">", "intron.txt"  or die $!;
    while (<IN>) {
        s/\.sam\n//;
        open IN1, "$_.gtf" or die "Cannot open file $_.gtf, $!\n";
        print OUT1 <IN1>;
        close IN1;
        if (-e "$_.intron") {
            open IN1, "$_.intron" or die $!;
            print OUT2 <IN1>;
            close IN1;
        }
    }
    close OUT2; close OUT1; close IN;

    # 将GTF文件转换成GFF3文件和transcripts序列
    # 若是非链特异性测序，则得到的single exon转录本序列是没有方向的。
    unless (-e "transfragDump.ok") {
        $cmdString = "$binPath/transfragDump transfrag.gtf ../$genome 2> transfrag.stats";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "transfragDump.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    # 对transcripts序列使用Transdecoder进行ORF分析
    unless (-e "TransDecoder.ok") {
        $cmdString = "$TransDecoder_LongOrfs  -m 100 -G universal -t transfrag.strand.fasta -S &> /dev/null";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString = "$TransDecoder_Predict  --retain_long_orfs 900  -t transfrag.strand.fasta &> /dev/null";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        $cmdString = "cp transfrag.strand.fasta.transdecoder.gff3 transfrag.transdecoder.gff3";
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        unless ($strand_specific) {
            $cmdString = "$TransDecoder_LongOrfs  -m 100 -G universal -t transfrag.noStrand.fasta &> /dev/null";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $cmdString = "$TransDecoder_Predict  --retain_long_orfs 900 -t transfrag.noStrand.fasta --train transfrag.strand.fasta.transdecoder_dir/longest_orfs.cds.top_500_longest &> /dev/null";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
            $cmdString = "cat transfrag.noStrand.fasta.transdecoder.gff3 >> transfrag.transdecoder.gff3";
            print STDERR (localtime) . ": CMD: $cmdString\n";
            system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        }
        open OUT, ">", "TransDecoder.ok" or die $!; close OUT;
    }
    else {
        $cmdString = "$TransDecoder_LongOrfs  -m 100 -G universal -t transfrag.strand.fasta -S &> /dev/null";
        print STDERR "CMD(Skipped): $cmdString\n";
        $cmdString = "$TransDecoder_Predict  --retain_long_orfs 900 -t transfrag.strand.fasta &> /dev/null";
        print STDERR "CMD(Skipped): $cmdString\n";
        $cmdString = "cp transfrag.strand.fasta.transdecoder.gff3 transfrag.transdecoder.gff3";
        print STDERR "CMD(Skipped): $cmdString\n";
        unless ($strand_specific) {
            $cmdString = "$TransDecoder_LongOrfs  -m 100 -G universal -t transfrag.noStrand.fasta &> /dev/null";
            print STDERR "CMD(Skipped): $cmdString\n";
            $cmdString = "$TransDecoder_Predict  --retain_long_orfs 900 -t transfrag.noStrand.fasta --train transfrag.strand.fasta.transdecoder_dir/longest_orfs.cds.top_500_longest &> /dev/null";
            print STDERR "CMD(Skipped): $cmdString\n";
            $cmdString = "cat transfrag.noStrand.fasta.transdecoder.gff3 >> transfrag.transdecoder.gff3";
            print STDERR "CMD(Skipped): $cmdString\n";
        }
    }

    # 将transcripts的ORF预测结果映射到基因组序列上，得到transcripts的基因预测结果： transfrag.genome.gff3
    $cmdString = "$binPath/transdecoder2ORF --out_protein proteins.fasta transfrag.gtf transfrag.transdecoder.gff3 ../$genome > transfrag.genome.gff3";
    unless (-e "transdecoder2ORF.ok") {
        print STDERR (localtime) . ": CMD: $cmdString\n";
        system("$cmdString") == 0 or die "failed to execute: $cmdString\n";
        open OUT, ">", "transdecoder2ORF.ok" or die $!; close OUT;
    }
    else {
        print STDERR "CMD(Skipped): $cmdString\n";
    }

    chdir "../";
    open OUT, ">", "3.transcript.ok" or die $!; close OUT;
}
else {
    print STDERR "Skip Step 3 for the file 3.transcript.ok exists\n";
}
