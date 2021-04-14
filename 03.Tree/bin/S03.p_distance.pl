#!/usr/bin/perl -w
 
use strict;

my %chrlist;
my $snpfile=shift;
my $dismatrix=shift;

####################################### database
my %dismat;
& dismatrix();
foreach my $nt1 (sort keys %dismat) {
	foreach my $nt2 (sort keys %{$dismat{$nt1}}) {
		print "$nt1\t$nt2\t$dismat{$nt1}{$nt2}\n";
	}
}

######################################### output
print "$snpfile\n";
print "$dismatrix\n";
open my $OUT,">$dismatrix" or die $!;

#########################################
open OUT ,">$dismatrix" or die $!;
my (%dis,%cns);
print "$snpfile start\n";
my $line=0;
open VCF,$snpfile=~/gz$/?"zcat $snpfile|":"<$snpfile" or die $!;
my @sams;
while (<VCF>) {
	chomp;
	if($_=~/^#/){
		if($_=~/\#CHROM/){
			my @tmp=split/\s+/,$_,10;
			@sams=split/\s+/,$tmp[9];
		}
		next;
	}
	
	my @a=split/\s+/,$_,10;
	my $chr=$a[0];
	my $pos=$a[1];
	my $ref=$a[3];
	my $alt=$a[4];
	my @nts=split/\s+/,$a[9];
	my %sam2nt;
	my %freq;
	for (my $i=0 ;$i<@nts ;$i++) {
		my $sam=$sams[$i];
		my $nts=(split/\:/,$nts[$i])[0];
		$nts=~s/\|/\//g;
		$sam2nt{$sam}=$nts;
		if($nts eq "./."){next};
		my @base1=split/\//,$nts;
		$freq{$base1[0]}++;
		$freq{$base1[1]}++;
		#print "$sam\t$nts\n";
	}
	my @allelenum=sort keys %freq;
	#if(scalar(@allelenum)!=2){next};
#	print "$chr\t$pos\t$ref\t$alt\t@allelenum\n";
	my $alles=join "",@allelenum;
	if($alles=~/3/){next};
	$line++;
	for (my $i=0;$i<@sams ;$i++) {
		my $sam1=$sams[$i];
		my $nt1=$sam2nt{$sam1};
		next if($nt1 =~ /\./);
		for (my $j=$i+1;$j<@sams ;$j++) {
			my $sam2=$sams[$j];
			my $nt2=$sam2nt{$sam2};
			next if($nt2 =~ /\./);
		#	print "$sam1\t$sam2\t$chr\t$pos\t$nt1\t$nt2\t$dismat{$nt1}{$nt2}\n";
		#	print "$sam1\t$sam2\t$chr\t$pos\n";
			if(!exists $dismat{$nt1}{$nt2}){
				print "unexists\t$sam1\t$sam2\t$chr\t$pos\t1:$nt1\t2:$nt2\n";
				next;
			}
			$dis{$sam1}{$sam2}+=$dismat{$nt1}{$nt2};
			$cns{$sam1}{$sam2}++;
		}
	}
	if($line%1000==0){
		print "$line\t$chr\t$pos\t".localtime()."\n";
	}
#	last if($line >5);
}
close VCF;

foreach my $sam1 (sort keys %dis) {
	foreach my $sam2 (sort keys %{$dis{$sam1}}) {
		print OUT "$sam1\t$sam2\t$dis{$sam1}{$sam2}\t$cns{$sam1}{$sam2}\n";
	}
}
close OUT;

############################# sub function
sub dismatrix{
	my %allele;
	$allele{"0\t0"}=0;
	$allele{"0\t1"}=0;
	$allele{"0\t2"}=0;
	$allele{"1\t1"}=0;
	$allele{"1\t2"}=0;
	$allele{"2\t2"}=0;
	foreach my $alle1 (sort keys %allele) {
		foreach my $alle2 (sort keys %allele) {
			my @alle1=split/\s+/,$alle1;
			my @alle2=split/\s+/,$alle2;
			my $dis=0;
			if($alle1[0] eq $alle1[1]){
				my @temp;
				@temp=@alle1;
				@alle1=@alle2;
				@alle2=@temp;
			}
			if($alle1[0] ne $alle2[0] and $alle1[0] ne $alle2[1]){
				$dis+=0.5;
			}
			if($alle1[1] ne $alle2[0] and $alle1[1] ne $alle2[1]){
				$dis+=0.5;
			}
			$dismat{"$alle1[0]/$alle1[1]"}{"$alle2[0]/$alle2[1]"}=$dis;
			$dismat{"$alle1[0]/$alle1[1]"}{"$alle2[1]/$alle2[0]"}=$dis;
			$dismat{"$alle1[1]/$alle1[0]"}{"$alle2[0]/$alle2[1]"}=$dis;
			$dismat{"$alle1[1]/$alle1[0]"}{"$alle2[1]/$alle2[0]"}=$dis;
			$dismat{"$alle2[0]/$alle2[1]"}{"$alle1[0]/$alle1[1]"}=$dis;
			$dismat{"$alle2[0]/$alle2[1]"}{"$alle1[1]/$alle1[0]"}=$dis;
			$dismat{"$alle2[1]/$alle2[0]"}{"$alle1[0]/$alle1[1]"}=$dis;
			$dismat{"$alle2[1]/$alle2[0]"}{"$alle1[1]/$alle1[0]"}=$dis;
		}
	}
}
 
