#!/usr/bin/perl -w
use strict;
my $snpnum=shift;
my $vcfin=shift;
my $outdir=shift;
my $max=$snpnum;

print "chrlist has been load\n";
my %snpturn;
my $turn=0;
while ($turn<$max) {
	$turn++;
	my $rand=int(rand($snpnum));
	$snpturn{$rand}++;
	if($turn%100000 ==0){
		print "$turn\n";
	}
}
print "the number $turn has been load\n";


my $num=0;
print "$outdir/PAV44.vcf.gz\n";
open VCF,$vcfin or die $!;
open OUT,"|gzip >$outdir/PAV44.vcf.gz" or die $!;
my $head=<VCF>;
print OUT $head;
while (<VCF>) {
	$num++;
	if(exists $snpturn{$num}){
		my $num2=$snpturn{$num};
		for (my $i=1;$i<=$num2 ;$i++) {
			print OUT "$_";
		}
	}
}
close VCF;
close OUT;
