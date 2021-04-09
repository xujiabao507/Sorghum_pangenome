use strict;
my $depthfile=shift;
my $output=shift;

open my $OUT,">$output" or die $!;
my %depmat;
my $sum=0;
my $num=0;
open my $DPH,$depthfile=~/gz$/?"zcat $depthfile|":"<$depthfile" or die $!;
while (<$DPH>) {
	chomp;
	my @a=split/\s+/;
	my $depth=int($a[2]+0.5);
#	if($depth>0 and $depth<400){
		$depmat{$depth}++;
		$num++;
		$sum+=$depth;
#	}
}
close $DPH;

my $mean=$sum/$num;

my $square=0;
my $plussub=0;
foreach my $depth (sort {$a<=>$b} keys %depmat) {
	next if($depth>2*$depth);
	$square+=((($depth-$mean)**2)*$depmat{$depth});
	my $rate_den=$depmat{$depth}/$sum;
	$plussub+=$depmat{$depth};
	my $rate_plus=$plussub/$sum;
	print $OUT "Depth:\t$depth\t$rate_den\t$rate_plus\t$depmat{$depth}\t$plussub\t$sum\n";
}
my $var=$square/$num;
my $sdev=$var**(0.5);
print $OUT "Summary:\t$num\t$sum\t$mean\t$var\t$sdev\n";
close $OUT;



