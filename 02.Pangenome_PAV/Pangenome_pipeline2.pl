use strict;
use Getopt::Long;
use File::Basename;


my $bin_file=shift;
my $project_file=shift;
my $fqs_file=shift;
my $sequel_file=shift;


###########  载入配置文件
my %bincf;
my %project;
my %genomes;
my %genomes2;
my %fqlist;
my %sequel;
&bin_load("$bin_file");
&project_load("$project_file");
&fqs_load("$fqs_file");
&sequel_load("$sequel_file");
################################################### outdir 设置
my $outdir=$project{"OutDir"};
my $shell_dir="$outdir/S01.Shell";
&creatdir($shell_dir);
my $db_dir="$outdir/S02.DB";
&creatdir($db_dir);
my $Mummer_dir="$outdir/S03.Mummer";
&creatdir($Mummer_dir);
my $refgenome=$project{ref};
my $Mummerfilt_dir="$outdir/S04.MumFilt_dir";
&creatdir($Mummerfilt_dir);
my $bwa_dir="$outdir/S05.BWA";
&creatdir($bwa_dir);
#my $vcf_dir="$outdir/S06.VCF";
#&creatdir($vcf_dir);
#my $gatk_dir="$outdir/S07.GATK";
#&creatdir($gatk_dir);
#my $cnv_dir="$outdir/S08.CNV";
#&creatdir($cnv_dir);
#my $pav_dir="$outdir/S09.PAV";
#&creatdir($pav_dir);
my $pav_dir2="$outdir/S10.PAV";
&creatdir($pav_dir2);

############################# Start
#&creatdir("$shell_dir/S01.split/");
#################  S01 split

my $sam1=0;
&creatdir("$shell_dir/S01.index_DB");
print "$shell_dir/S01.index_DB\n";
foreach my $sam (sort keys %genomes) {
	#next if($sam ne "S273" && $sam ne "S353");
	$sam1=$sam;
	&creatdir("$db_dir/$sam");
	open my $SH01a,">$shell_dir/S01.index_DB/step01.$sam.index.sh" or die $!;
	my $genome=$genomes{$sam};
	if(exists $genomes2{$sam}){
		$genome=$genomes2{$sam};
	}
	print $SH01a "perl $bincf{bindir}/script/step01.genome_split.pl $genome $sam $db_dir \n";
	print $SH01a "ln -s $genomes{$sam} $db_dir/$sam.genome.fa\n";
	print $SH01a "$bincf{bwa}/bwa index $db_dir/$sam.genome.fa\n";
	print $SH01a "$bincf{samtools}/samtools faidx $db_dir/$sam.genome.fa\n";
	print $SH01a "$bincf{java} -jar $bincf{picard}/CreateSequenceDictionary.jar R=$db_dir/$sam.genome.fa  OUTPUT=$db_dir/$sam.genome.dict\n";
	close $SH01a;
}
my $base=(split/\//,$refgenome)[-1];
open my $SH01b,">$shell_dir/S01.index_DB/step01.ref.index.sh" or die $!;
print $SH01b "ln -s $refgenome $db_dir/$base \n";
$refgenome="$db_dir/$base";
print $SH01b "$bincf{bwa}/bwa index $db_dir/$base \n";
print $SH01b "$bincf{samtools}/samtools faidx $db_dir/$base\n";
if($base=~/fa$/ or $base=~/fasta$/){
	my @base=split/\./,$base;
	my $numbase=scalar(@base);
	my $base2=join ".",@base[0..$numbase-2];
	print $SH01b "$bincf{java} -jar $bincf{picard}/CreateSequenceDictionary.jar R=$db_dir/$base OUTPUT=$db_dir/$base2.dict\n";
	#print   "$bincf{java} -jar $bincf{picard}/CreateSequenceDictionary.jar R=$db_dir/$base OUTPUT=$db_dir/$base2.dict\n";
}else{
	print $SH01b "$bincf{java} -jar $bincf{picard}/CreateSequenceDictionary.jar R=$db_dir/$base OUTPUT=$db_dir/$base.dict\n";
	#print "$bincf{java} -jar $bincf{picard}/CreateSequenceDictionary.jar R=$db_dir/$base OUTPUT=$db_dir/$base.dict\n";
}

close $SH01b;

################   S02 nucmer
my @fapre=glob "$db_dir/$sam1/*fasta";
if(scalar(@fapre)==0){
	print "please run step01\n";
	exit;
}

#### 比对
foreach my $sam (sort keys %genomes) {
	#next if($sam ne "S273" && $sam ne "S353");
	&creatdir("$shell_dir/S02.nucmer/$sam");
	print "$shell_dir/S02.nucmer/$sam\n";
	&creatdir("$Mummer_dir/$sam");
	my @fa=glob "$db_dir/$sam/*fasta";
	my $fanum=scalar(@fa);
	print "$sam\t$fanum\t$db_dir/$sam\n";
	foreach my $fa (@fa) {
		my $base=(split/\//,$fa)[-1];
		my $output_nucmer="$Mummer_dir/$sam/$base";
		open my $SH02,">$shell_dir/S02.nucmer/$sam/step02.$base.sh" or die $!;
		#print "$shell_dir/S02.nucmer/$sam/step02.$base.sh\n";
		print $SH02 "$bincf{mummer}/nucmer --maxmatch -c 50 -l 40  -t 6 $refgenome $fa  -p $output_nucmer\n";
		close $SH02;
	}
}
#### 过滤
foreach my $sam (sort keys %genomes) {
	#next if($sam ne "S273" && $sam ne "S353");
	&creatdir("$shell_dir/S03.filter/$sam");
	print "$shell_dir/S03.filter/$sam\n";
	my @fa=glob "$db_dir/$sam/*fasta";
	foreach my $fa (@fa) {
		my $base=(split/\//,$fa)[-1];
		my $output_nucmer="$Mummer_dir/$sam/$base.delta";
		my $output_filt="$Mummer_dir/$sam/$base.filt";
		my $output_coords="$Mummer_dir/$sam/$base.coords";
		my $output_filt_all="$Mummer_dir/$sam/$base.filt.all";
		open my $SH03,">$shell_dir/S03.filter/$sam/step03.$base.sh" or die $!;
		print $SH03 "$bincf{mummer}/delta-filter -1 -i 50 -l 50 -u 10 -o 100 $output_nucmer >$output_filt\n";
		print $SH03 "$bincf{mummer}/delta-filter -i 50 -l 50 -u 10 -o 100 $output_nucmer >$output_filt_all\n";
		#print $SH03 "$bincf{mummer}/show-coords -c -d -I 0.5 -l -L 1000 -r $output_filt >$output_coords\n";
		close $SH03;
	}
}
#### merge coords
print "$shell_dir/S04.Merge_coord\n";
foreach my $sam (sort keys %genomes) {
	#next if($sam ne "S273" && $sam ne "S353");
	&creatdir("$shell_dir/S04.Merge_coord");
	&creatdir("$Mummerfilt_dir/$sam");
	open my $SH04,">$shell_dir/S04.Merge_coord/step04.$sam.merge.sh" or die $!;
	my $delta_onebyone_merge="$Mummerfilt_dir/$sam/$sam.filt.pre";
	my $coords_onebyone_merge="$Mummerfilt_dir/$sam/$sam.filt.pre.coords";
	my $delta_onebyone_ok="$Mummerfilt_dir/$sam/$sam.obo_ok";
	my $snps_onebyone_merge="$Mummerfilt_dir/$sam/$sam.filt.snps";
	print $SH04 "perl $bincf{bindir}/script/step04.merge_filt.pl $Mummer_dir/$sam $refgenome $genomes2{$sam} $delta_onebyone_merge\n";
	print $SH04 "$bincf{mummer}/show-coords -c -d -I 0.5 -l -L 100 -r $delta_onebyone_merge >$coords_onebyone_merge\n";
	print $SH04 "perl $bincf{bindir}/script/step05.stat_obo.pl $coords_onebyone_merge  $coords_onebyone_merge.select\n";
	print $SH04 "perl $bincf{bindir}/script/step06.obo_filt_delta.pl   $coords_onebyone_merge.select $delta_onebyone_merge $delta_onebyone_ok\n";
	print $SH04 "$bincf{mummer}/show-snps -r -T -C $delta_onebyone_ok >$snps_onebyone_merge \n";
	print $SH04 "perl $bincf{bindir}/script/step07.stat.pl $coords_onebyone_merge $coords_onebyone_merge.stat\n";
	print $SH04 "perl $bincf{bindir}/script/step07.stat.pl $coords_onebyone_merge.select $coords_onebyone_merge.select.stat\n";
	print $SH04 "perl $bincf{bindir}/script/step06.1.snp_indel_out.pl  $snps_onebyone_merge\n";
	print $SH04 "";
	close $SH04;
}
######## SV inversion + translocation
#
#foreach my $sam (sort keys %genomes) {
#	open my $SH051,">$shell_dir/S04.Merge_coord/step05.$sam.Inversion_Transloaction.sh" or die $!;
#	print $SH051 "perl $bincf{bindir}/script/step07.1.inversion_translocation.pl $Mummerfilt_dir/$sam/$sam.filt.pre.coords.select $Mummerfilt_dir/$sam/$sam.SV\n";
#	print $SH051 "perl $bincf{bindir}/script/step07.1.1.inversion_translocation_filt.pl $Mummerfilt_dir/$sam/$sam.SV.TransLocation\n";
#	print $SH051 "perl $bincf{bindir}/script/step07.1.1.inversion_translocation_filt.pl $Mummerfilt_dir/$sam/$sam.SV.Inversion\n";
#	close $SH051;
#}

exit;


######### short read alignment

print "$shell_dir/S05.BWA\n";
foreach my $sam (sort keys %genomes) {
	&creatdir("$shell_dir/S05.BWA");
	&creatdir("$bwa_dir/$sam");
	if(!exists  $fqlist{$sam}{1}){print "$sam\n";next;};
	open my $SH05,">$shell_dir/S05.BWA/step05.$sam.bwa.self.sh" or die $!;
	my $self_genome="$db_dir/$sam.genome.fa";
	my $bam_pre_self="$bwa_dir/$sam/$sam.self.pre.bam";
	my $bam_sort_self="$bwa_dir/$sam/$sam.self.sort.bam";
	my $bam_rmdup_self="$bwa_dir/$sam/$sam.self.rmdup.bam";
	print $SH05 "$bincf{bwa}/bwa mem -t 8 -M  -R \"\@RG\\tID:$sam\\tSM:$sam\\tCN:BGI\\tPL:illumina\" $self_genome $fqlist{$sam}{1} $fqlist{$sam}{2}|$bincf{samtools}/samtools view -Sb  -o $bam_pre_self \n";
	print $SH05 "$bincf{samtools}/samtools  sort -m 4000M $bam_pre_self -o $bam_sort_self \n";
	&creatdir("$bwa_dir/$sam/temp1/");
	print $SH05 "$bincf{java} -Xmx3g -Djava.io.tmpdir=$bwa_dir/$sam/temp1/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{picard}/MarkDuplicates.jar MAX_FILE_HANDLES=1000 REMOVE_DUPLICATES=false TMP_DIR=$bwa_dir/$sam/temp1/ M=$bwa_dir/$sam/$sam.rmdup.metrics I=$bam_sort_self O=$bam_rmdup_self \n";
	print $SH05 "$bincf{java} -Xmx2G -jar $bincf{picard}/BuildBamIndex.jar I=$bam_rmdup_self\n";
	open my $SH052,">$shell_dir/S05.BWA/step05.$sam.bwa_hq.self.sh" or die $!;
	print $SH052 "$bincf{samtools}/samtools view -q 20 -h -b  $bam_rmdup_self -o $bwa_dir/$sam/$sam.self.hq.bam\n";
	print $SH052 "$bincf{samtools}/samtools index  $bwa_dir/$sam/$sam.self.hq.bam\n";
	close $SH052;

	close $SH05;
	open my $SH06,">$shell_dir/S05.BWA/step06.$sam.bwa.ref.sh" or die $!;
	my $bam_pre_ref="$bwa_dir/$sam/$sam.ref.pre.bam";
	my $bam_sort_ref="$bwa_dir/$sam/$sam.ref.sort.bam";
	my $bam_rmdup_ref="$bwa_dir/$sam/$sam.ref.rmdup.bam";
	print $SH06 "#$bincf{bwa}/bwa mem -t 8 -M -R \"\@RG\\tID:$sam\\tSM:$sam\\tCN:BGI\\tPL:illumina\" $refgenome $fqlist{$sam}{1} $fqlist{$sam}{2}|$bincf{samtools}/samtools view -Sb    -o $bam_pre_ref \n";
	print $SH06 "#$bincf{samtools}/samtools  sort -m 4000M $bam_pre_ref -o $bam_sort_ref \n";
	&creatdir("$bwa_dir/$sam/temp2/");
	print $SH06 "$bincf{java} -Xmx3g -Djava.io.tmpdir=$bwa_dir/$sam/temp2/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{picard}/MarkDuplicates.jar MAX_FILE_HANDLES=1000 REMOVE_DUPLICATES=false TMP_DIR=$bwa_dir/$sam/temp2/ M=$bwa_dir/$sam/$sam.rmdup.metrics I=$bam_sort_ref O=$bam_rmdup_ref \n";
	print $SH06 "$bincf{java} -Xmx2G -jar $bincf{picard}/BuildBamIndex.jar I=$bam_rmdup_ref \n";
	close $SH06;
}


#foreach my $sam (sort keys %genomes) {
#	if(exists $sequel{$sam}){
#		open my $SH05,">$shell_dir/S05.BWA/step05.$sam.bwa.Sequel.sh" or die $!;
#		my $bam_pre_sequel_ref="$bwa_dir/$sam/$sam.sequel.ref.pre.bam";
#		my $bam_sort_sequel_ref="$bwa_dir/$sam/$sam.sequel.sort.bam";
#		print $SH05 "$bincf{bwa}/bwa mem -t 8 -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 $refgenome $sequel{$sam} |$bincf{samtools}/samtools view -Sb  -o $bam_pre_sequel_ref\n";
#		print $SH05 "$bincf{samtools}/samtools sort -m 4000M  $bam_pre_sequel_ref -o $bam_sort_sequel_ref\n";
#		print $SH05 "perl $bincf{bindir}/script/len.pl $sequel{$sam} >$bwa_dir/$sam/$sam.sequel.len\n";
#		print $SH05 "perl $bincf{bindir}/script/match2.pl $bwa_dir/$sam/$sam.sequel.len  $bam_sort_sequel_ref $bam_sort_sequel_ref.split\n";
#		
#		my $self_genome="$db_dir/$sam.genome.fa";
#		my $bam_pre_sequel_self="$bwa_dir/$sam/$sam.sequel.self.pre.bam";
#		my $bam_sort_sequel_self="$bwa_dir/$sam/$sam.sequel.self.sort.bam";
#		my $depthfile="$bam_sort_sequel_self.split.map.depth.gz";
#		my $read_cov_file="$bam_sort_sequel_self.split.map.read_cov.gz";
#		my $read_cov_reg_file="$bam_sort_sequel_self.split.map.read_cov_Reg.gz";
#		print $SH05 "$bincf{bwa}/bwa mem -t 8 -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 $self_genome  $sequel{$sam}  |$bincf{samtools}/samtools view -Sb  -o  $bam_pre_sequel_self\n";
#		print $SH05 "$bincf{samtools}/samtools sort -m 4000M $bam_pre_sequel_self  -o $bam_sort_sequel_self\n";
#		print $SH05 "perl $bincf{bindir}/script/match2.pl $bwa_dir/$sam/$sam.sequel.len  $bam_sort_sequel_self  $bam_sort_sequel_self.split\n";
#		print $SH05 "$bincf{samtools}/samtools depth $bam_sort_sequel_self.split.map.bam |gzip >$depthfile\n";
#		print $SH05 "$bincf{samtools}/samtools mpileup -O --output-QNAME -a -V -f $self_genome $bam_sort_sequel_self.split.map.bam |gzip >$read_cov_file\n";
#		print $SH05 "perl $bincf{bindir}/script/step05.2.read_cov.pl $read_cov_file $read_cov_reg_file\n";
#		print $SH05 "perl $bincf{bindir}/script/step07.2.inversion_translocation_cov.pl $read_cov_reg_file $Mummerfilt_dir/$sam/$sam.SV.Inversion $Mummerfilt_dir/$sam/$sam.SV.TransLocation  $Mummerfilt_dir/$sam/$sam.SV.read_support\n";
#		print $SH05 "perl $bincf{bindir}/script/step07.3.inversion_translocation_stat.pl $Mummerfilt_dir/$sam/$sam.SV.read_support $Mummerfilt_dir/$sam/$sam.SV.read_support.stat\n";
#		close $SH05;
#	}
#}

############################ samtools VCF
#print "$shell_dir/S06.VCF\n";
#foreach my $sam (sort keys %genomes) {
#	&creatdir("$shell_dir/S06.VCF");
#	&creatdir("$vcf_dir/$sam");
#	open my $SH06,">$shell_dir/S06.VCF/step06.$sam.vcf.ref.sh" or die $!;
#	my $bam_sort_ref="$bwa_dir/$sam/$sam.ref.sort.bam";
#	my $bam_qual_ref="$bwa_dir/$sam/$sam.ref.qual.bam";
#	my $vcf_ref_pre="$vcf_dir/$sam/$sam.pre.vcf.gz";
#	print $SH06 "perl $bincf{bindir}/script/step08.bam_quality.pl $bam_sort_ref $bincf{samtools}/samtools|$bincf{samtools}/samtools view -Sb - -o $bam_qual_ref\n";
#	print $SH06 "$bincf{samtools}/samtools mpileup -g -d 1000 -q 20 -Q 15 -t AD -f $refgenome $bam_qual_ref |$bincf{bcftools}/bcftools call -c -O z - -o $vcf_ref_pre \n";
#	close $SH06;
#}
#foreach my $sam (sort keys %genomes) {
#	&creatdir("$shell_dir/S06.VCF");
#	&creatdir("$vcf_dir/$sam");
#	open my $SH06,">$shell_dir/S06.VCF/step06.$sam.vcf.self.sh" or die $!;
#	my $bam_sort_self="$bwa_dir/$sam/$sam.self.sort.bam";
#	my $bam_qual_self="$bwa_dir/$sam/$sam.self.qual.bam";
#	my $vcf_sort_self="$vcf_dir/$sam/$sam.self.vcf.gz";
#	my $self_genome="$db_dir/$sam.genome.fa";
#	print $SH06 "perl $bincf{bindir}/script/step08.bam_quality.pl $bam_sort_self $bincf{samtools}/samtools|$bincf{samtools}/samtools view -Sb - -o $bam_qual_self\n";
#	print $SH06 "$bincf{samtools}/samtools mpileup -g -d 1000 -q 20 -Q 15 -t AD -f $self_genome $bam_qual_self |$bincf{bcftools}/bcftools call -c -O z - -o $vcf_sort_self \n";
#	close $SH06;
#}
###############################   samtools vcf filter
#foreach my $sam (sort keys %genomes) {
#	&creatdir("$shell_dir/S06.VCF");
#	&creatdir("$vcf_dir/$sam");
#	open my $SH06,">$shell_dir/S06.VCF/step07.$sam.vcf_filter.self.sh" or die $!;
#	print "$shell_dir/S06.VCF/step07.$sam.vcf_filter.self.sh\n";
#	my $vcf_ref_pre="$vcf_dir/$sam/$sam.self.vcf.gz";
#	my $vcf_ref_filt="$vcf_dir/$sam/$sam.var_filt.vcf.gz";
#	print $SH06 "zcat $vcf_ref_pre|perl /zfssz5/BC_PUB/Software/03.Soft_ALL/bcftools-1.5/bin/vcfutils.pl varFilter -Q 30 -d 10 -D 500 -a 1 -w 5  |gzip >$vcf_ref_filt \n";
#	close $SH06;
#}

################################ GATK 


#foreach my $sam (sort keys %genomes) {
#	&creatdir("$shell_dir/S07.GATK/$sam");
#	print "$shell_dir/S07.GATK/$sam\n";
#	open my $SH07,">$shell_dir/S07.GATK/$sam/step05.$sam.self.sh" or die $!;
#	my $self_genome="$db_dir/$sam.genome.fa";
##	my $bam_rmdup_self="$bwa_dir/$sam/$sam.self.rmdup.bam";
#	my $bam_rmdup_self="$bwa_dir/$sam/$sam.self.hq.bam";
#	my $gvcf_pri="$gatk_dir/$sam/$sam.pri.gvcf.gz";
#	my $vcf_pri="$gatk_dir/$sam/$sam.pri.vcf.gz";
#	my $indel_pri="$gatk_dir/$sam/$sam.indel.pri.vcf.gz";
#	my $intervals="$gatk_dir/$sam/$sam.indel.intervals";
#	my $bam_realign="$gatk_dir/$sam/$sam.realignment.bam";
#	my $gvcf_realignment="$gatk_dir/$sam/$sam.realign.gvcf.gz";
#	my $vcf_realignment="$gatk_dir/$sam/$sam.realign.vcf.gz";
#	my $snp_raw="$gatk_dir/$sam/$sam.snp_raw.vcf.gz";
#	my $snp_filt="$gatk_dir/$sam/$sam.snp_filt.vcf.gz";
#	my $snp_result="$gatk_dir/$sam/$sam.snp_result.vcf.gz";
#	my $indel_raw="$gatk_dir/$sam/$sam.indel_raw.vcf.gz";
#	my $indel_filt="$gatk_dir/$sam/$sam.indel_filt.vcf.gz";
#	my $indel_result="$gatk_dir/$sam/$sam.indel_result.vcf.gz";
#	&creatdir("$gatk_dir/$sam/temp/");
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T HaplotypeCaller -R $self_genome -I $bam_rmdup_self --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $gvcf_pri -nct 16 \n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $self_genome --variant $gvcf_pri -o $vcf_pri -stand_call_conf 30 -allSites\n";
#	print $SH07  "time /hwfssz4/BC_PUB/Software/03.Soft_ALL/jdk1.8.0_202/bin/java -XX:-UseGCOverheadLimit -Xmx4G -Djava.io.tmpdir=$gatk_dir/$sam/temp/  -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar  HaplotypeCaller -R  $self_genome -I $bam_rmdup_self --sample-ploidy 2 --standard-min-confidence-threshold-for-calling 30 --output-mode EMIT_ALL_SITES --native-pair-hmm-threads 16 -ERC GVCF -O $gvcf_pri\n";
#	print $SH07 "time /hwfssz4/BC_PUB/Software/03.Soft_ALL/jdk1.8.0_202/bin/java -XX:-UseGCOverheadLimit -Xmx4G -Djava.io.tmpdir=$gatk_dir/$sam/temp/  -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar  GenotypeGVCFs -R $self_genome -V  $gvcf_pri -O $vcf_pri\n";
#	print $SH07 "perl  /hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S07.GATK/stat/gatk-snp_het_pri.pl $sam\n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T SelectVariants  -selectType INDEL --excludeNonVariants  -R $self_genome -V $vcf_pri -o $indel_pri \n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $self_genome -I $bam_rmdup_self -known $indel_pri -o $intervals \n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T IndelRealigner  -R $self_genome -I $bam_rmdup_self -known $indel_pri -targetIntervals $intervals -o $bam_realign\n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T HaplotypeCaller -R $self_genome -I $bam_realign --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $gvcf_realignment -nct 2 \n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $self_genome --variant $gvcf_realignment -o $vcf_realignment -stand_call_conf 30 -allSites\n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T SelectVariants  -selectType SNP --excludeNonVariants  -R $self_genome -V $vcf_realignment -o $snp_raw\n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T VariantFiltration --filterExpression  \"QUAL>30.0 && QD > 20.0 && FS < 20.000 && MQ > 40.0 \" --filterName \"NEXTOK\" -R $self_genome -V $snp_raw -o $snp_filt\n";
#	#print $SH07 "zgrep NEXTOK $snp_filt |gzip >$snp_result\n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T SelectVariants  -selectType INDEL --excludeNonVariants  -R $self_genome -V $vcf_realignment -o $indel_raw\n";
#	#print $SH07 "$bincf{java} -Xmx4g -Djava.io.tmpdir=$gatk_dir/$sam/temp/ -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $bincf{gatk}/GenomeAnalysisTK.jar -T VariantFiltration --filterExpression  \"QUAL>30.0 && QD > 20.0 && FS < 20.000 && MQ > 40.0 \" --filterName \"NEXTOK\" -R $self_genome -V $indel_raw -o $indel_filt\n";
#	#print $SH07 "zgrep NEXTOK $indel_filt |gzip >$indel_result\n";
#	#print $SH07 "perl $bincf{bindir}/script/step06.2.snp_indel_stat.pl $outdir $sam $outdir/S07.GATK/$sam/$sam.small_var_stat.xls\n";
#	close $SH07;
#}



############################################   S08 CNV
##my $cnv_dir="";
#foreach my $sam (sort keys %genomes) {
#	#&creatdir("$shell_dir/S08.CNV/$sam");
#	print "$shell_dir/S08.CNV/\n";
#	&creatdir("$cnv_dir/$sam");
#	my $bam_rmdup_ref="$bwa_dir/$sam/$sam.ref.rmdup.bam";
#	my $depth_ref="$cnv_dir/$sam/$sam.ref.depth.gz";
#	my $genedepth_ref="$cnv_dir/$sam/$sam.gene.depth";
#	my $genedepth_ref_detail="$cnv_dir/$sam/$sam.gene.depth.detail.gz";
#	my $ntdepth_stat_ref="$cnv_dir/$sam/$sam.nt.depth.stat";
#	my $genedepth_stat_ref="$cnv_dir/$sam/$sam.gene.depth.stat";
#	my $bam_rmdup_self="$bwa_dir/$sam/$sam.self.rmdup.bam";
#	my $depth_self="$cnv_dir/$sam/$sam.self.depth.gz";
#	open my $SH081,">$shell_dir/S08.CNV/step081.$sam.depth.ref.sh" or die $!;
#	print $SH081 "#$bincf{samtools}/samtools depth -a $bam_rmdup_ref |gzip >$depth_ref \n";
#	print $SH081 "perl $bincf{bindir}/script/step11.gene_depth.pl $project{refgff}  $depth_ref $genedepth_ref\n";
#	print $SH081 "perl $bincf{bindir}/script/step11.gene_depth_detail.pl $project{refgff}  $depth_ref $genedepth_ref_detail\n";
#	print $SH081 "perl $bincf{bindir}/script/step12.stat.pl $depth_ref $ntdepth_stat_ref\n";
#	print $SH081 "perl $bincf{bindir}/script/step12.gene_stat.pl  $genedepth_ref_detail  $genedepth_stat_ref\n";
#	print $SH081 "perl $bincf{bindir}/script/step13.gene_p.pl $genedepth_ref_detail $genedepth_ref\n";
#	close $SH081;
#	#open my $SH082,">$shell_dir/S08.CNV/$sam/step082.depth.self.sh" or die $!;
#	#print $SH082 "$bincf{samtools}/samtools depth -Q 20 -a $bam_rmdup_self |gzip >$depth_self \n";
#	#close $SH082;
#}

###########################################   S09 PAV1
#foreach my $sam (sort keys %genomes) {
#	last;
#	&creatdir("$shell_dir/S09.PAV/$sam");
#	&creatdir("$pav_dir/$sam");
#	print "$shell_dir/S09.PAV\n";
#	open my $SH09,">$shell_dir/S09.PAV/$sam/step09.unmap.$sam.sh" or die $!;
#	my $map="$Mummerfilt_dir/$sam/$sam.filt.pre.coords.select";
#	my $unmap="$pav_dir/$sam/step01.$sam.unmap.fa";
#	&creatdir("$pav_dir/$sam/split_unmap");
#	print $SH09 "perl $bincf{bindir}/script/step14.unmap_read.pl $map $genomes{$sam} $unmap\n";
#	print $SH09 "$bincf{blast}/formatdb -i $refgenome -p F \n";
#	print $SH09 "perl $bincf{bindir}/script/step01.genome_split.pl $unmap $sam $pav_dir/$sam \n";
#	close $SH09;
#
#	my $blast="$pav_dir/$sam/step01.$sam.unmap.blast";
#	my $blast_filt="$pav_dir/$sam/step01.$sam.unmap.blast.filt";
#	
#	open my $SH10a,">$shell_dir/S09.PAV/$sam/step10.blast.$sam.sh" or die $!;
#	print $SH10a "#$bincf{blast}/blastall -p blastn -d $refgenome -i $unmap -a 4 -e 10.0 -m 8 -o $blast \n";
#	print $SH10a "perl $bincf{bindir}/script/step15.map_scale.pl $map $blast $blast_filt \n";
#	close $SH10a;
#
#	&creatdir("$pav_dir/$sam/blast");
#	my @subfa=glob "$pav_dir/$sam/$sam/*fasta";
#	my $cat="cat ";
#	foreach my $subfa (@subfa) {
#		my $base=(split/\//,$subfa)[-1];
#		&creatdir("$shell_dir/S09.PAV/$sam/step10.blast");
#		open my $SH10,">$shell_dir/S09.PAV/$sam/step10.blast/step10.$base.sh" or die $!;
#		print $SH10 "#$bincf{blast}/blastall -p blastn -d $refgenome -i $subfa -a 4 -e 10.0 -m 8 -o $pav_dir/$sam/blast/$base.m8 \n";
#		print $SH10 "perl $bincf{bindir}/script/step15.map_scale.pl $map $pav_dir/$sam/blast/$base.m8 $pav_dir/$sam/blast/$base.m8.filt \n";
#		$cat.=" $pav_dir/$sam/blast/$base.m8.filt ";
#		close $SH10;
#	}
#	open my $SH11,">$shell_dir/S09.PAV/$sam/step11.cat_filt.sh" or die $!;
#	my $blast_merge="$pav_dir/$sam/blast/$sam.merge.filt";
#	my $blast_map_largest="$pav_dir/$sam/blast/$sam.map_largest.xls";
#	my $blast_rate="$pav_dir/$sam/blast/$sam.blast_rate.xls";
#	my $blast_unmap="$pav_dir/$sam/blast/$sam.blast_unmap.xls";
#	print $SH11 "$cat >$blast_merge \n";
#	close $SH11;
#
#	open my $SH12,">$shell_dir/S09.PAV/$sam/step12.map_result.sh" or die $!;
#	print $SH12 "perl $bincf{bindir}/script/step16.blast_map_filt.pl $blast_merge  $blast_map_largest \n";
#	print $SH12 "perl $bincf{bindir}/script/step17.unmap2.pl $unmap $blast_map_largest $blast_rate\n";
#	#print $SH12 "awk \'\$4<0.6\' $blast_rate >$blast_unmap\n";
#	my $mapresult="$pav_dir/$sam/blast/$sam.map_result";
#	print $SH12 "perl $bincf{bindir}/script/step18.allmap_sample.pl $genomes{$sam} $map $blast_map_largest $blast_rate $mapresult\n";
#	my $contig_unmap="$pav_dir/$sam/blast/$sam.contig.unmap";
#	my $chr_unmap="$pav_dir/$sam/blast/$sam.chr.unmap";
#	print $SH12 "perl $bincf{bindir}/script/step19.specific-map.pl  $genomes{$sam} $mapresult $contig_unmap \n";
#	print $SH12 "perl $bincf{bindir}/script/step20.specific-map.pl  $refgenome $mapresult $chr_unmap \n";
#	print $SH12 "perl $bincf{bindir}/script/step20.unmap_merge.pl $chr_unmap $chr_unmap.merge\n";
#	print $SH12 "perl $bincf{bindir}/script/step20.unmap_merge.pl $contig_unmap $contig_unmap.merge\n";
#	close $SH12;
#}

###########################################   S10 PAV2

&creatdir("$shell_dir/S10.PAV");
&creatdir("$pav_dir2/db");

open my $SH13,">$shell_dir/S10.PAV/step13.windows_fasta.sh" or die $!;
print "$shell_dir/S10.PAV\n";
print $SH13 "perl $bincf{bindir}/script/step21.window_fa.pl $refgenome $pav_dir2/db/ref.winddow_fa.gz\n";
foreach my $sam (sort keys %genomes) {
	&creatdir("$pav_dir2/$sam");
	&creatdir("$shell_dir/S10.PAV/$sam");
	print $SH13 "perl $bincf{bindir}/script/step21.window_fa.pl $genomes{$sam} $pav_dir2/db/$sam.winddow_fa.gz\n";
}
close $SH13;

foreach my $sam (sort keys %genomes) {
	&creatdir("$pav_dir2/$sam");
	&creatdir("$shell_dir/S10.PAV/$sam");
	open my $SH14,">$shell_dir/S10.PAV/$sam/step14.$sam.bwa.sh" or die $!;
	my $fq="$pav_dir2/db/$sam.winddow_fa.gz";
	my $prebam="$pav_dir2/$sam/$sam.pre.bam";
	my $sortbam="$pav_dir2/$sam/$sam.sort.bam";
	print $SH14 "$bincf{bwa}/bwa mem -t 16 -w 500 -M $refgenome $fq |$bincf{samtools}/samtools view -Sb -o $prebam\n";
	print $SH14 "$bincf{samtools}/samtools sort -m 4000M  $prebam -o $sortbam\n";
	close $SH14;
	open my $SH15,">$shell_dir/S10.PAV/$sam/step15.$sam.ref_bwa.sh" or die $!;
	my $sam_genome="$db_dir/$sam.genome.fa";
	my $ref_fq="$pav_dir2/db/ref.winddow_fa.gz";
	my $pre_refbam="$pav_dir2/$sam/$sam.ref.pre.bam";
	my $sort_refbam="$pav_dir2/$sam/$sam.ref.sort.bam";
	my $sort_mapbam="$pav_dir2/$sam/$sam.map.bam";
	my $sort_refmapbam="$pav_dir2/$sam/$sam.refmap.bam";
	print $SH15 "$bincf{bwa}/bwa mem -t 16 -w 500 -M $sam_genome $ref_fq |$bincf{samtools}/samtools view -Sb -o $pre_refbam \n";
	print $SH15 "$bincf{samtools}/samtools sort -m 4000M  $pre_refbam -o $sort_refbam\n";
	close $SH15;
	open my $SH16,">$shell_dir/S10.PAV/$sam/step16.$sam.map2unmap.sh" or die $!;
	print "$shell_dir/S10.PAV/$sam/step16.$sam.map2unmap.sh\n";
	print $SH16 perl $bincf{bindir}/script/step22.filt2.pl $fq $sortbam $sort_mapbam $pav_dir2/$sam/$sam.map.reads.gz  $pav_dir2/$sam/$sam.unmap.reads\n";
	print $SH16 "$bincf{samtools}/samtools index $sort_mapbam\n";
	print $SH16 "perl $bincf{bindir}/script/step09.map_regoin_merge.pl $pav_dir2/$sam/$sam.map.reads.gz  $pav_dir2/$sam/$sam.mapregoin.merge.xls\n";
	print $SH16 "perl $bincf{bindir}/script/step10.map2unmap.pl  $outdir/S02.DB/$sam.genome.fa.fai $pav_dir2/$sam/$sam.mapregoin.merge.xls $pav_dir2/$sam/$sam.unmapregoin.xls\n";
	print $SH16 "perl $bincf{bindir}/script/step11.unmap_unmapread.pl  $pav_dir2/$sam/$sam.unmapregoin.xls $pav_dir2/$sam/$sam.unmap.reads $pav_dir2/$sam/step11.$sam.unmap_reads_merge.xls\n";
	close $SH16;
	open my $SH17,">$shell_dir/S10.PAV/$sam/step17.$sam.ref_depth.sh" or die $!;
	print "$shell_dir/S10.PAV/$sam/step17.$sam.ref_depth.sh\n";
	print $SH17 "$bincf{samtools}/samtools depth -a $sort_mapbam >$pav_dir2/$sam/step01.$sam.depth\n";
	print $SH17 "perl $bincf{bindir}/script/cov_depth.pl $pav_dir2/$sam/step01.$sam.depth $pav_dir2/$sam/step02.$sam.cov \n";
	print $SH17 "perl $bincf{bindir}/script/step12.window_fa_N.pl $fq $pav_dir2/$sam/step02.$sam.winddow_fa.list \n";
	print $SH17 "perl $bincf{bindir}/script/step13.window_edge.pl $pav_dir2/$sam/step02.$sam.winddow_fa.list $pav_dir2/$sam/step03.$sam.gap.xls\n";
	print $SH17 "perl $bincf{bindir}/script/step14.maploc.pl $pav_dir2/$sam/$sam.map.reads.gz $pav_dir2/$sam/step03.$sam.gap.xls $pav_dir2/$sam/step04.$sam.maplocation.xls\n";
	print $SH17 "perl $bincf{bindir}/script/step15.cluster_inter.pl $pav_dir2/$sam/step04.$sam.maplocation.xls $pav_dir2/$sam/step05.$sam.gap.xls  $pav_dir2/$sam/step05.$sam.map.contig.xls\n";
	print $SH17 "perl $bincf{bindir}/script/map.pl $pav_dir2/$sam/step05.$sam.gap.xls  $pav_dir2/$sam/step05.$sam.gap_real.xls\n";
	print $SH17 "perl $bincf{bindir}/script/step16.ref_specific.pl $pav_dir2/$sam/step05.$sam.gap_real.xls $pav_dir2/$sam/step02.$sam.cov  $pav_dir2/$sam/step07.$sam.unmap.xls\n";
	print $SH17 "perl $bincf{bindir}/script/cov_depth_merge.pl $pav_dir2/$sam/step07.$sam.unmap.xls >$pav_dir2/$sam/step08.$sam.unmapmerge.xls\n";
	print $SH17 "perl $bincf{bindir}/script/step22.filt2.pl $ref_fq $sort_refbam $sort_refmapbam $pav_dir2/$sam/$sam.refmap.reads.gz $pav_dir2/$sam/$sam.ref2unmap.reads\n";
	print $SH17 "perl $bincf{bindir}/script/step11.unmap_unmapread.pl $pav_dir2/$sam/step08.$sam.unmapmerge.xls $pav_dir2/$sam/$sam.ref2unmap.reads $pav_dir2/$sam/step08.$sam.unmapmerge_lowmap.xls\n";
	print $SH17 "perl $bincf{bindir}/script/step22.ref_map_bam_filt.pl  $bwa_dir/$sam/$sam.ref.rmdup.bam  $pav_dir2/$sam/$sam.ref.lib500.bam \n";
	print $SH17 "$bincf{samtools}/samtools  index $pav_dir2/$sam/$sam.ref.lib500.bam \n";
	print $SH17 "$bincf{samtools}/samtools  depth $pav_dir2/$sam/$sam.ref.lib500.bam >$pav_dir2/$sam/step09.$sam.reflib500.depth \n";
	print $SH17 "perl $bincf{bindir}/script/step23.pav_NGS_deletion.pl $pav_dir2/$sam/step08.$sam.unmapmerge_lowmap.xls $pav_dir2/$sam/step09.$sam.reflib500.depth $pav_dir2/$sam/step10.$sam.pav_cov.xls\n";
	print $SH17 "perl $bincf{bindir}/script/step23.pav_seq.pl $pav_dir2/$sam/step10.$sam.pav_cov.xls  $pav_dir2/$sam/step11.$sam.pav_result.xls \n";
 	close $SH17;
	open my $SH18,">$shell_dir/S10.PAV/$sam/step18.$sam.map_det.sh" or die $!;
	print "$shell_dir/S10.PAV/$sam/step18.$sam.map_det.sh\n";
	print $SH18 "perl $bincf{bindir}/script/step17.map_read.pl  $pav_dir2/$sam/$sam.map.reads.gz  $pav_dir2/$sam/step09.$sam.mapstat.xls\n";
	print $SH18 "perl $bincf{bindir}/script/map2file_merge.pl $pav_dir2/$sam/step09.$sam.mapstat.xls $pav_dir2/$sam/step09.$sam.mapstat_merge.xls\n";
	print $SH18 "perl $bincf{bindir}/script/map2file_merge_other.pl $pav_dir2/$sam/step09.$sam.mapstat_merge.xls $pav_dir2/$sam/step09.$sam.mapstat.xls $pav_dir2/$sam/step10.$sam.map_detail.xls\n";
	print $SH18 "perl $bincf{bindir}/script/step18.pav_out.pl $pav_dir2/$sam/step10.$sam.map_detail.xls $pav_dir2/$sam/step11.$sam.mapstat.xls\n";
	print $SH18 "perl $bincf{bindir}/script/step19.order.pl $pav_dir2/$sam/step11.$sam.mapstat.xls $pav_dir2/$sam/step12.$sam.order.xls  $pav_dir2/$sam/step13.$sam.SVs.xls\n";
	print $SH18 "perl $bincf{bindir}/script/step23.pav_detect.pl $pav_dir2/$sam/$sam.map.reads.gz  $pav_dir2/$sam/step12.$sam.order.xls  $pav_dir2/$sam/step14.$sam.SV_read.xls\n";
	print $SH18 "perl $bincf{bindir}/script/step24.pav_detect_cluster.pl  $pav_dir2/$sam/step14.$sam.SV_read.xls   $pav_dir2/$sam/step15.$sam.svs.xls $pav_dir2/$sam/step15.$sam.svs_detail.xls\n"; 
	print $SH18 "perl $bincf{bindir}/script/step23.sv_fa.pl  $outdir/S02.DB/Sbicolor_454_v3.0.1.fa $pav_dir2/$sam/step15.$sam.svs.xls $pav_dir2/$sam/step20.$sam.read.fa\n";
	print $SH18 "perl $bincf{bindir}/script/step23.sv_fq.pl  $outdir/S02.DB/Sbicolor_454_v3.0.1.fa $pav_dir2/$sam/step15.$sam.svs.xls $pav_dir2/$sam/step20.$sam.read.fq\n";
	print $SH18 "$bincf{bwa}/bwa mem -t 32 $outdir/S02.DB/$sam.genome.fa   $pav_dir2/$sam/step20.$sam.read.fq  > $pav_dir2/$sam/step22.$sam.read.bwa_mem.sam\n";
	print $SH18 "perl $bincf{bindir}/script/step24.map_uniq.pl step20.$sam.read.fa  step22.$sam.read.bwa_mem.sam ste22.$sam.map_uniq.xls\n";
	close $SH18;

}



#/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S05.BWA


############################# End










################## 子函数
sub sequel_load{
	my ($file)=@_;
	if(!-e $file){print  "sequel file is un-exists\n";return 0};
	open my $FL,$file or die $!;
	while (<$FL>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		next if($_!~/\w/);
		my @a=split/\s+/;
		my $sam=$a[0];
		my $seq=$a[1];
		$sequel{$sam}=$seq;
	}
	close $FL;
}

sub fqs_load{
	my ($file)=@_;
	open my $FL,$file or die $!;
	while (<$FL>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		next if($_!~/\w/);
		my @a=split/\s+/;
		my $sam=$a[0];
		my $fq1=$a[1];
		my $fq2=$a[2];
		my $insertsize=$a[3];
		$fqlist{$sam}{1}=$fq1;
		$fqlist{$sam}{2}=$fq2;
		$fqlist{$sam}{3}=$insertsize;
		print "$sam\t$fq1\t$fq2\n";
	}
	close $FL;
}
sub bin_load{
	my ($file)=@_;
	open my $FL,$file or die $!;
	while (<$FL>) {
		chomp;
		$_=~s/\s+//g;
		my @a=split/\=/;
		$bincf{$a[0]}=$a[1];
	}
	close $FL;
}


sub project_load{
	my ($file)=@_;
	open my $FL,$file or die $!;
	while (<$FL>) {
		chomp;
		next if($_=~/^#/);
		$_=~s/(^\s+|\s+$)//g;
		my @a=split/\=/;
		if($a[0]=~/(.*)\_genome/){
			my $sam=$1;
			$a[1]=~s/(^\s+|\s+$)//g;
			my @g=split/\s+/,$a[1];
			if(scalar(@g)==1){
				$genomes{$sam}=$g[0];
			}else{
				$genomes{$sam}=$g[0];
				$genomes2{$sam}=$g[1];
			}
			print "1\t$sam\t@g\n";
			#print "$1\n";
			next;
		}
		$a[0]=~s/\s+//g;
		$project{$a[0]}=$a[1];
	}
	close $FL;
}
 
sub creatdir{
	my ($dir)=@_;
	if(!-e $dir){
		`mkdir -p $dir`;
		print "$dir\n";
	}
}





