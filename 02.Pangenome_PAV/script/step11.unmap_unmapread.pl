use strict;
my $unmapreg=shift;
my $unmapreads=shift;
my $mapout=shift;

print "$unmapreg\n";
my %pavs;
open UMR,$unmapreg or die $!;
while (<UMR>) {
	chomp;
	my @a=split/\s+/,$_;
	my $chr1=$a[0];
	my $st1=$a[1];
	my $en1=$a[2];
	$pavs{$chr1}{$st1}=$en1;
}
close UMR; 

print "regoin caculation\n";
my %depth;
foreach my $chr (sort keys %pavs) {
	foreach my $st (sort keys %{$pavs{$chr}}) {
		for (my $loc=$st;$loc<=$pavs{$chr}{$st} ;$loc++) {
			$depth{$chr}{$loc}=0;
		}
	}
}


print "$unmapreads load\n";
open OUT,">$mapout" or die $!;
open READ,$unmapreads or die $!;
my %scale;
my %unmapread;
while (<READ>) {
	chomp;
	my @a=split/\s+/;
	my @b=split/\_/,$a[0];
	my ($chr1,$st1,$en1);
	if($a[0]=~/super/){
		$chr1="$b[0]\_$b[1]";
		$st1=$b[2];
		$en1=$b[3];
	}
	if($a[0]=~/Chr/){
		$chr1=$b[0];
		$st1=$b[1];
		$en1=$b[2];
	}
	my $rate=$a[1];
	if($rate<0.25){
		for (my $loc=$st1;$loc<=$en1 ;$loc++) {
			if(exists $depth{$chr1}{$loc}){
				$depth{$chr1}{$loc}=1;
			}
		}
	}
}
close READ; 


print "Output caculation\n";
my %covs;
foreach my $chr (sort keys %pavs) {
	foreach my $st (sort {$a<=>$b} keys %{$pavs{$chr}}) {
		my $sum=0;
		my $end=$pavs{$chr}{$st};
		for (my $loc=$st;$loc<=$end ;$loc++) {
			if($depth{$chr}{$loc}>0){
				$sum++;
			}
		}
		my $len=$end-$st+1;
		#my $rate=int(($sum/$len)*100000+0.5)/10000;
		my $rate=$sum/$len;
		print OUT "$chr\t$st\t$end\t$len\t$rate\t$sum\n";
	}
}
close OUT;









