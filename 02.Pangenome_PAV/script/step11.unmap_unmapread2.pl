use strict;
my $unmapreg=shift;
my $unmapreads=shift;
my $mapout=shift;

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

my %allread;
foreach my $chr (sort keys %pavs) {
	foreach my $st (sort keys %{$pavs{$chr}}) {
		for (my $loc=$st;$loc<=$pavs{$chr}{$st} ;$loc++) {
			if($loc%100==0){
				$allread{$chr}{$loc}=$st;
			}
		}
	}
}

open OUT,">$mapout" or die $!;
open READ,$unmapreads or die $!;
my %scale;
my %unmapread;
while (<READ>) {
	chomp;
	my @a=split/\s+/;
	my @b=split/\_/,$a[0];
	my $rate=$a[1];
	my ($chr1,$st1,$en1);
	if($a[0]=~/Chr/){
		$chr1=$b[0];
		$st1=$b[1];
		$en1=$b[2];
	}else{
		$chr1="$b[0]\_$b[1]";
		$st1=$b[2];
		$en1=$b[3];
	}
	if($rate<0.25){
		if(exists $allread{$chr1}{$st1}){
			my $loc=$allread{$chr1}{$st1};
			$scale{$chr1}{$loc}++;
			$unmapread{$chr1}{$loc}{$st1}=0;
			$unmapread{$chr1}{$loc}{$en1}=0;
		}
	}
}
close READ; 


foreach my $chr (sort keys %scale) {
	my @st=sort {$a<=>$b} keys %{$scale{$chr}};
	foreach my $st (@st) {
		if($scale{$chr}{$st}>0){
			my @loc=sort {$a<=>$b} keys %{$unmapread{$chr}{$st}};
			my $en=$pavs{$chr}{$st};
			my $len=$en-$st+1;
			print OUT "$chr\t$st\t$en\t$len\t$scale{$chr}{$st}\t$loc[0]\t$loc[-1]\n";
		}
	}
}

close OUT;
