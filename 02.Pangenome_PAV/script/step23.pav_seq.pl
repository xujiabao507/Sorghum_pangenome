use strict;

my $pav_dep=shift;
my $out=shift;
open OUT,">$out" or die $!;

print "$pav_dep is loading\n";
my %dep;
my $sum=0;
open PAV,$pav_dep or die $!;
while (<PAV>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	#last if($chr ne "Chr01");
	my $loc=$a[1];
	my $dep=$a[4];
	if($dep==0){
		$dep{$chr}{$loc}=$dep;
		$sum++;
	#	print "$chr\t$loc\t$dep\n";
	}
	#print "$chr\t$loc\t$dep\n";
}
close PAV;
print  "Sum of ref-specific $sum\n";

my %pav;
foreach my $chr (sort keys %dep) {
	my @loc=sort {$a<=>$b} keys %{$dep{$chr}};
	my $st;
	my $en;
	for (my $i=0;$i<(scalar(@loc)-1) ;$i++) {
		my $loc2=$loc[$i];
		my $loc1=$loc2-1;
		my $loc3=$loc2+1;
		if(!exists $dep{$chr}{$loc1} and exists $dep{$chr}{$loc2}){
			$st=$loc2;
		}
		if(exists $dep{$chr}{$loc2} and !exists $dep{$chr}{$loc3}){
			$en=$loc2;
			$pav{$chr}{$st}=$en;
		}
	}
}
foreach my $chr (sort keys %pav) {
	foreach my $st (sort {$a<=>$b} keys %{$pav{$chr}}) {
		my $en=$pav{$chr}{$st};
		my $len=$en-$st+1;
		print OUT "$chr\t$st\t$en\t$len\n";
	}
}
close OUT;












