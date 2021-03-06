use strict;
 
use lib "/hwfssz4/BC_PUB/Software/03.Soft_ALL/PerlInfo/share/perl5/";
use Statistics::R;

my $stat=shift;
my $genefile=shift;
my ($mu,$sdev);
#open STAT,"$stat" or die $!;
#while (<STAT>) {
#	chomp;
#	next if($_!~/Summary/);
#	my @a=split/\s+/;
#	$mu=$a[3];
#	$sdev=$a[5];
#}
#close STAT;
#print "$mu\t$sdev\n";

my $Rbin = '/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/bin/R';
my $R = Statistics::R->new( r_bin => $Rbin);

my $sum=0;
my $num=0;
open my $FL0,"zcat $stat|" or die $!;
while(<$FL0>){
	chomp;
	my @a=split/\s+/;
	my $mean_g=$a[2];
	$sum+=$mean_g;
	$num++;
}
close $FL0;
my $mu=$sum/$num;
print "$sum\t$num\t$mu\n";
#####################################
#####################################
open my $FL,"$genefile" or die $!;
my $num=0;
open my $OUT,">$genefile.copynum.xls" or die $!;
while (<$FL>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $gene=$a[1];
	my $mean_g=$a[2];
	my $sdev_g=$a[3];
	my $samn_g=$a[4];
	if($mean_g<$mu*0.2 or $mean_g>$mu*1.8){
		#my $t=($mean_g-$mu)/($sdev/$samn_g**0.5);
		#my $df=$samn-1;
		#$R->set('t',"$t");
		#$R->set('df',"$df");
		#$R->run(q`dt(t,df)->pval`);
		#$R->run(q`pnorm(t,0,1)->pval`);
		#my $pval = $R->get('pval');
		my $copynum=int($mean_g/$mu+0.5);
		#print "$chr\t$gene\t$mu\t$mean\t$sdev\t$samn\t$pval\n";
		#if($pval>0.5){
		#	$pval=1-$pval;
		#}
		#if($pval<0.00001){
			$mean_g=int($mean_g+0.5);
			my $mean_all=int($mu+0.5);
			print $OUT "$chr\t$gene\t$mean_g\t$mean_all\t$copynum\n";
		#}
	}
}
close $FL;




