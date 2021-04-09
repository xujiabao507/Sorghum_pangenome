use strict;

my $pav=shift;
my $depth=shift;
my $out=shift;

print "$pav is loading\n";
open OUT,">$out" or die $!;
my %pavs;
my %pavscale;
open PAV,$pav or die $!;
while (<PAV>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
#	last if($chr ne "Chr01");
	my $st=$a[1];
	my $en=$a[2];
	my $rate=$a[4];
	next if($rate<0.5);
	for (my $loc=$st;$loc<=$en ;$loc++) {
		$pavs{$chr}{$loc}=0;
		$pavscale{$chr}{$st}=$en;
		#print "$chr\t$loc\t$st\t$en\n";
	}
}
close PAV;

print "$depth is loading\n";
my %pavdel;
open DEP,$depth or die $!;
while (<DEP>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
#	last if($chr ne "Chr01");
	my $loc=$a[1];
	my $dep=$a[2];
	next if($dep==0);
	if(exists $pavs{$chr}{$loc}){
		$pavs{$chr}{$loc}=$dep;
	}
}
close DEP;

print "$out doing\n";
foreach my $chr(sort keys %pavscale) {
	foreach my $st (sort {$a<=>$b} keys %{$pavscale{$chr}}) {
		my $en=$pavscale{$chr}{$st};
		my $len=$en-$st+1;
		my $covlen=0;
		for (my $loc=$st;$loc<=$en ;$loc++) {
			if($pavs{$chr}{$loc}>0){
				#$covlen++;
				print OUT  "$chr\t$loc\t$st\t$en\t$pavs{$chr}{$loc}\n";
			}else{
				print OUT  "$chr\t$loc\t$st\t$en\t0\n";
			}
		}
	}
}
close OUT;















