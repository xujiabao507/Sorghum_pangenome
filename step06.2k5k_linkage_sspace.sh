/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa mem -t 16 02.db/IS3614.genome.fa HiSeq/2k/IS3614-2/180101_I137_FCHM3HCBBXX_L7_WHSORkpzDCAADWAAPEI-52_1.fq.gz.clean.gz HiSeq/2k/IS3614-2/180101_I137_FCHM3HCBBXX_L7_WHSORkpzDCAADWAAPEI-52_2.fq.gz.clean.gz|/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h - -b -o IS3614/IS3614.2K_L7_WHSORkpzDCAADWAAPEI-52.pri.bam
perl S30.Hic/08.2k_5k_linkage/bin/script/bwa_mem_uniq.pl IS3614/IS3614.2K_L7_WHSORkpzDCAADWAAPEI-52.pri.bam IS3614/IS3614.2K_L7_WHSORkpzDCAADWAAPEI-52.uniq.bam
perl S30.Hic/08.2k_5k_linkage/bin/script/uniq_reads.pl S30.Hic/08.2k_5k_linkage/Result/03.bwa/S273/S273.2K.uniq.bam Input/ipa/Lib_2000bp_1.fq.gz Input/ipa/Lib_2000bp_2.fq.gz S30.Hic/08.2k_5k_linkage/Result/04.uniq/S273_2K
perl bin/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l  SSPACE.AusTRCF317961.lib -s  02.db/AusTRCF317961.genome.fa -x 0 -k 5 -a 0.7 -n 15 -z 0 -T 20 -b AusTRCF317961




