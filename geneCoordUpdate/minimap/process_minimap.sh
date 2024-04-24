/share/mg6t/caix/src/Minimap2/minimap2-2.10_x64-linux/minimap2 -cx asm5 -t 30 Brapa_genome_v3.0_1801.fasta Brapa_scaffold_v3.fasta > aln.paf 
sort -k6,6 -k8,8n aln.paf > aln.paf.sort
