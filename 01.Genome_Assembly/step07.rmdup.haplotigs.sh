/share/backup/xujiabao/Software/BWA/bwa-0.7.16a/bwa index  Pilon_genome.fasta
/share/backup/xujiabao/Software/BWA/bwa-0.7.16a/bwa mem -t 40 Pilon_genome.fasta  input//Lib_270bp_1.fq.gz input//Lib_270bp_2.fq.gz|/share/backup/xujiabao/Software/Samtools/samtools-1.3.1/samtools view -b -S - |/share/backup/xujiabao/Software/Samtools/samtools-1.3.1/samtools sort - -o Lib_270.sort.bam
/share/backup/xujiabao/Software/Samtools/samtools-1.3.1/samtools index Lib_270.sort.bam  
export PATH="/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Deal_Assemble/purge_haplotigs/purge_haplotigs/bin:/zfssz5/BC_PUB/Software/03.Soft_ALL/R-3.4.1/bin:/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Deal_Assemble/purge_haplotigs/usr/bin/:/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Bin/ncbi-blast-2.7.1+/bin/:/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Bin/lastz/:/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Bin/samtools/bin/:$PATH" 
/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Deal_Assemble/purge_haplotigs/purge_haplotigs/bin/purge_haplotigs readhist  Lib_270.sort.bam 
/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Deal_Assemble/purge_haplotigs/purge_haplotigs/bin/purge_haplotigs contigcov -i Lib_270.sort.bam.genecov -l 20 -m 45 -h 80
/hwfssz4/BC_PUB/pipeline/DNA/DNA_Denovo/Deal_Assemble/purge_haplotigs/purge_haplotigs/bin/purge_haplotigs purge -g  Pilon_genome.fasta -c coverage_stats.csv  -b Lib_270.sort.bam  -t 16  
