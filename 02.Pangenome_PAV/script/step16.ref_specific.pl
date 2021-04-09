use strict;

my $inter=shift;
my $map=shift;
my $out=shift;

my %svinter;
open INT,$inter or die $!;
while (<INT>) {
	chomp;
	if($_=~/NA/){next};
	my @a=split/\s+/;
	my $chr=$a[0];
	my $st=$a[1];
	my $en=$a[2];
	for (my $loc=$st;$loc<=$en ;$loc++) {
		$svinter{$chr}{$loc}=0;
	}
}
close INT;

open OUT,">$out" or die  $!;
open MAP,$map or die $!;
while (<MAP>) {
	chomp;
	if($_=~/NA/){next};
	my @a=split/\s+/;
	my $chr=$a[0];
	my $st=$a[1];
	my $en=$a[2];
	my %unmap;
	for (my $loc=$st;$loc<=$en ;$loc++) {
		if(!exists $svinter{$chr}{$loc}){
			$unmap{$loc}=0;
		}
	}
	my @loc=sort {$a<=>$b} keys %unmap;
	my %scale;
	my $st=0;
	my $en=0;
	foreach my $loc (@loc) {
		my $loc1=$loc-1;
		my $loc2=$loc;
		my $loc3=$loc+1;
		if(!exists $unmap{$loc1} and exists $unmap{$loc2} and exists $unmap{$loc3}){
			$st=$loc2;
			$scale{$st}=0;
		}
		if(exists $unmap{$loc1} and exists $unmap{$loc2} and !exists $unmap{$loc3}){
			$en=$loc2;
			$scale{$st}=$en;
		}
	}
	my $sum=0;
	foreach my $st (sort {$a<=>$b} keys %scale) {
		my $len=($scale{$st}-$st)+1;
		$sum+=$len;
		print OUT "$chr\t$st\t$scale{$st}\t$len\t$sum\n";
	}
}
close MAP;
close OUT;







