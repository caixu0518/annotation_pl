#!/usr/bin/perl -w
use strict;

my $gff = shift;
my $sample = shift;

my $gene = "";
my $gene_flag = 1; #������Ŀ
my @mRNA;  #��¼ת¼����Ŀ����GeMoMA��ƥ��prediction
my @CDS;   #��¼ĳ���򹲰������ٸ�CDS������CDS��ID��ӦmRNA
my @exon;  #��¼ĳ�����������exon��Ŀ
my @count_cds; #��¼ת¼����Ӧ����ʼ��ֹλ��
my $flag = 0;  #��¼ÿ��ת¼��������exon��CDS����Ŀ

my ($name, $other) = split(/\./, $sample, 2);

open GFF, "<$gff" or die "$!";
while(my $line = <GFF>){
	chomp $line;
	next if($line=~/^#/);
	my @eles = split /\t/, $line;
	my @temp = splice(@eles, 8);
	$eles[1] = "GeMoMa.$name";
	if($eles[2] eq "gene"){
		my $gene = join("\t", @eles);
		push @count_cds, $flag;
		print join("\t", @eles, "ID=gene_$name\_$gene_flag;Name=GeMoMa.gene_$name\_$gene_flag;")."\n";
		for(my $k=0; $k<=$#mRNA; $k++){
			print join("\t",$mRNA[$k], "ID=gene_$name\_$gene_flag.m$k;Parent=gene_$name\_$gene_flag;")."\n";
			for(my $i=$count_cds[$k]; $i<$count_cds[$k+1]; $i++){
				print join("\t", $exon[$i], "ID=$name\_$gene_flag.exon.$i;Parent=gene_$name\_$gene_flag.m$k;")."\n";
				print join("\t", $CDS[$i], "ID=cds_$name\_$gene_flag.$k;Parent=gene_$name\_$gene_flag.m$k;")."\n";
			}
		}
		@mRNA = ();
		@count_cds = ();
		@exon = ();
		@CDS = ();
		$flag = 0;
		$gene_flag++;
	}else{
		if($eles[2] eq "prediction"){
			$eles[2] = "mRNA";
			push @mRNA, join("\t", @eles);
			push @count_cds, $flag;
		}
		if($eles[2] eq "CDS"){
			$flag++;
			push @CDS, join("\t", @eles);
			$eles[2] = "exon";
			push @exon, join("\t",@eles);
		}
	}
}
close GFF;
