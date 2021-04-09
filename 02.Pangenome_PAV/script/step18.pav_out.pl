use strict;

my $map=shift;
my $out=shift;
open OUT,">$out" or die  $!;
print "$map start\n";
my %readmap;
my %scale;
open MAP,$map or die $!;
while (<MAP>) {
	chomp;
	my @a=split/\s+/;
	my $qchr=$a[0];
	my $qst=$a[1];
	my $tchr=$a[4];
	my $tst=$a[5];
	my $stall=$a[6];
	my $enall=$a[7];
	if($qchr ne $tchr){next};
	$readmap{$qchr}{$stall}{$qst}="$tst";
	$scale{$qchr}{$stall}=$enall;
#	last if($qchr ne "Chr01");
}
close MAP;

foreach my $chr (%scale) {
	foreach my $stall (sort {$a<=>$b} keys %{$scale{$chr}}) {
		my @loc=sort {$a<=>$b} keys %{$readmap{$chr}{$stall}};
		my %main100k;
		my %main100kcont;
		foreach my $loc(@loc) {
			my $tloc=$readmap{$chr}{$stall}{$loc};
			my $tk100=int($tloc/100000)*100000;
			$main100k{$tk100}++;
			$main100kcont{$tk100}{$tloc}=0;
		}
		my @k100=sort {$main100k{$b}<=>$main100k{$a}} keys %main100k;
		my $main100k=$k100[0];
		#my $num=0;
		#foreach my $k100 (@k100) {
		#	$num++;
		#	print "$chr\t$stall\t$scale{$chr}{$stall}\t$num\t$k100\t$main100k{$k100}\n";
		#}
		my $sum=0;
		my $num=0;
		foreach my $tlocs (sort {$a<=>$b} keys %{$main100kcont{$main100k}}) {
			$sum+=$tlocs;
			$num++;
		}
		my $mean=$sum/$num;
		my @loc2;
		my $len=abs($loc[-1]-$loc[0])+500;
		foreach my $loc (@loc) {
			if(abs(($readmap{$chr}{$stall}{$loc}-$mean))<$len){
				push @loc2,$loc;
			}
		}
		my $minx=$loc2[0];
		my $maxx=$loc2[-1];
		my $maxx2=$maxx+500;
		my $miny=$readmap{$chr}{$stall}{$minx};
		my $maxy=$readmap{$chr}{$stall}{$maxx}+500;
		my $len1=$maxx2-$minx;
		my $len2=$maxy-$miny;
		if($miny=~/\d/ and $maxy=~/\d/ ){
			print OUT "$chr\t$stall\t$len1\t$len2\t$minx\t$maxx2\t$miny\t$maxy\n";
		}
		#foreach my $loc (@loc2) {
		#	my $dis=$readmap{$chr}{$stall}{$loc}-$mean;
		#	print OUT "$chr\t$loc\t$readmap{$chr}{$stall}{$loc}\t$minx\t$maxx2\t$miny\t$maxy\n";
		#}
	}
}
close OUT;

