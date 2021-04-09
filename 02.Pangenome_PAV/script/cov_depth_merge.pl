use strict;


my $uncov=shift;

my %scale;
open UC,$uncov or die $!;
while (<UC>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $st=$a[1];
	my $en=$a[2];
	$scale{$chr}{$st}=$en;
}
close UC;

foreach my $chr (sort keys %scale) {
	my %cluster;
	my @loc=sort {$a<=>$b} keys %{$scale{$chr}};
	for (my $i=0;$i<@loc ;$i++) {
		my $st1=$loc[$i];
		my $en1=$scale{$chr}{$st1};
		my $len=$en1-$st1+1;
		my $lensum=$len;
		$cluster{$st1}="$en1\t$len\t$len\t1";
		for (my $j=$i+1;$j<@loc ;$j++) {
			my $st2=$loc[$j];
			my $en2=$scale{$chr}{$st2};
			$lensum+=($en2-$st2+1);
			my $allsum=$en2-$st1+1;
			my $rate=$lensum/$allsum;
			if($rate>0.9){
				$cluster{$st1}="$en2\t$lensum\t$allsum\t$rate";
				$i=$j;
			}else{
				last;
			}
		}
	}
	foreach  my $st(sort {$a<=>$b} keys %cluster) {
		print "$chr\t$st\t$cluster{$st}\n";
	}
}













