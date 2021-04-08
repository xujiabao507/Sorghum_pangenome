/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa bwasw -t 16 11.Map_linkage/db/Sbicolor_454_v3.0.1.fa  step05.sspace/AusTRCF317961/AusTRCF317961/AusTRCF317961.final.scaffolds.fasta |/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h -  -b -o 11.Map_linkage/03.bwa/AusTRCF317961.bwasw.contig.bam

perl /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/11.Map_linkage/bin/scripts/contig_split.pl step05.sspace/AusTRCF317961/AusTRCF317961/AusTRCF317961.final.scaffolds.fasta |gzip >AusTRCF317961.fq.gz

/zfssz5/BC_PUB/Software/03.Soft_ALL/bwa-0.7.17/bwa mem -t 16 11.Map_linkage/db/Sbicolor_454_v3.0.1.fa 11.Map_linkage//02.fqs/AusTRCF317961.fq.gz|/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h - -b -o 11.Map_linkage//03.bwa/AusTRCF317961.bwa.mem.bam

perl /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/11.Map_linkage/bin/scripts/bwa_best.pl 11.Map_linkage//02.fqs/AusTRCF317961.fq.gz 11.Map_linkage//03.bwa/AusTRCF317961.bwa.mem.bam 11.Map_linkage//03.bwa/AusTRCF317961.best.gz


perl /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/11.Map_linkage/bin/scripts/bwa_best_scaffold_linkage.pl  11.Map_linkage//03.bwa/AusTRCF317961.best.gz >11.Map_linkage//04.linkage/AusTRCF317961.best.list

perl /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/11.Map_linkage/bin/scripts/unmap_scf.pl 04.linkage/AusTRCF317961.best.list /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/08.2k_5k_linkage/Result/01.shell/step05.sspace/AusTRCF317961/AusTRCF317961/AusTRCF317961.final.scaffolds.fasta  04.linkage/AusTRCF317961.unmap.fa


/opt/bio/ncbi/bin/formatdb -i 04.linkage/AusTRCF317961.unmap.fa -p F 
/opt/bio/ncbi/bin/blastall -p blastn -m 8 -a 4 -e 1e-10 -d 04.linkage/AusTRCF317961.unmap.fa -i /hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/DB/Sbicolor_454_v3.1.1.cds.fa -o 04.linkage/AusTRCF317961.unmap.m8

perl /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/11.Map_linkage/bin/scripts/blast_gene_order_same.pl 04.linkage/AusTRCF317961.unmap.gene.order.list 04.linkage/AusTRCF317961.unmap.gene.order.list.index
perl /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/11.Map_linkage/bin/scripts/lianj.pl 04.linkage/AusTRCF317961.best.list /zfssz4/BC_COM_P4/F17FTSNCKF1543/SORipaD/S30.Hic/08.2k_5k_linkage/Result/01.shell/step05.sspace/AusTRCF317961/AusTRCF317961/AusTRCF317961.final.scaffolds.fasta 04.linkage/AusTRCF317961.unmap.gene.order.list.index.2 /04.linkage/AusTRCF317961.map.fa