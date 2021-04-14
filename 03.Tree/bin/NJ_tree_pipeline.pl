use strict;


 
my $vcfin=shift;
my $outdir=shift;
my $bindir=shift;
my $bootstrap=shift;
 
$bootstrap||=100;
&creatdir("$outdir/Tree_shell");

my $snpnum=`wc -l  $vcfin`;
$snpnum=$snpnum-1;
print "$outdir/Tree_shell/step02.randsnp/step02.snprand snprand\n";
&creatdir("$outdir/S02.SNP");
for (my $bt=1;$bt<=$bootstrap ;$bt++) {
	&creatdir("$outdir/S02.SNP/Bootstrap_$bt");
	&creatdir("$outdir/Tree_shell/step02.randsnp");
	open  CH02,">$outdir/Tree_shell/step02.randsnp/step02.snprand.$bt.sh";
	print CH02 "perl $bindir/S02.bootstrap.pl $snpnum $vcfin $outdir/S02.SNP/Bootstrap_$bt\n";
	close CH02;
}

print "$outdir/Tree_shell/step04.p_distance.sh\n";
open SH04,">$outdir/Tree_shell/step04.p_distance.sh" or die $!;
&creatdir("$outdir/S03.Distance");
&creatdir("$outdir/S04.Distance_merge");
&creatdir("$outdir/Tree_shell/step03.distancelist");
for (my $bt=1;$bt<=$bootstrap ;$bt++) {
	&creatdir("$outdir/S03.Distance/Bootstrap_$bt");
	open  CH03,">$outdir/Tree_shell/step03.distancelist/step03.$bt.p_distance.sh";
	print CH03 "perl $bindir/PAV_tree_reseq.pl  $outdir/S02.SNP/Bootstrap_$bt/PAV44.vcf.gz $outdir/S03.Distance/Bootstrap_$bt/PAV44.dismatrix\n";
	close CH03;
}
close SH04;

############## S05 PHYLIP

print "$outdir/Tree_shell/step05.phylip.sh\n";
open SH05,">$outdir/Tree_shell/step05.phylip.sh" or die $!;
print SH05 "cd $outdir/S05.PHYLIP/\n";
&creatdir("$outdir/S05.PHYLIP");
for (my $bt=1;$bt<=$bootstrap ;$bt++) {
	open CONF,">$outdir/S05.PHYLIP/phylip.$bt.conf" or die $!;
	print CONF "$outdir/S03.Distance/Bootstrap_$bt/PAV44.dismatrix\n";
	print CONF "N\n";
	print CONF "Y\n";
	close CONF;
	print SH05 "/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S05.SL/SNP/02.PopGentic/bin/Tree/phylip-3.696/exe/neighbor <$outdir/S05.PHYLIP/phylip.$bt.conf\n";
	print SH05 "mv outfile nj.$bt.outfile\n";
	print SH05 "mv outtree nj.$bt.outtree\n";
}
close SH05;


print "$outdir/Tree_shell/step06.cat_NJ.sh\n";
&creatdir("$outdir/S06.Consense_Tree");
open SH06,">$outdir/Tree_shell/step06.cat_NJ.sh" or die $!;
my $cmd="cat ";
for (my $bt=1;$bt<=$bootstrap ;$bt++) {
	$cmd.=" $outdir/S05.PHYLIP/nj.$bt.outtree ";
}
print SH06 "$cmd >$outdir/S06.Consense_Tree/NJ.merge.outtree\n";
close SH06;

open CONF,">$outdir/S06.Consense_Tree/conf" or die $!;
print CONF "$outdir/S06.Consense_Tree/NJ.merge.outtree\nY\n";
close CONF;
#
#
print "$outdir/Tree_shell/step07.Consense_NJ.sh\n";
open SH07,">$outdir/Tree_shell/step07.Consense_NJ.sh" or die $!;
print SH07 "cd $outdir/S06.Consense_Tree/\n";
print SH07 "/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S05.SL/SNP/02.PopGentic/bin/Tree/phylip-3.696/exe/consense <$outdir/S06.Consense_Tree/conf\n";
print SH07 "mv outfile consense.outfile\n";
print SH07 "mv outtree consense.tre\n";
close SH07;
#
#
#print "$outdir/Tree_shell/step08.normal_name.sh\n";
#open SH08,">$outdir/Tree_shell/step08.normal_name.sh" or die $!;
#print SH08 "perl $bindir/NJ/S05.back2name.pl $outdir/S04.Distance_merge/Bootstrap_1.dismatrix.namespair  $outdir/S06.Consense_Tree/consense.tre  $outdir/S06.Consense_Tree/NormalName.tre\n";
#close SH08;
 
##############################################
sub creatdir{
	my ($dir)=@_;
	if(!-e $dir){
		`mkdir -p $dir`;
		print "$dir\n";
	}
}

