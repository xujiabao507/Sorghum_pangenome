use strict;
my $map =shift;
my $mapout=shift;
open OUT ,">$mapout" or die $!;
my %pavs;
open MAP,$map or die $!;
while (<MAP>) {
	chomp;
	my @a=split/\s+/;
	my $chr1=$a[0];
	my $st1=$a[1];
	my $en1=$a[2];
	my $chr2=$a[4];
	my $st2=$a[5];
	my $cigar=$a[3];
	#if($chr1 ne "Chr01"){next};
	if($chr1 eq $chr2){
		$pavs{$chr1}{$st1}=$en1;
	}
}
close MAP; 


open OUT,">$mapout" or die $!;
my %scale;
foreach my $chr (sort keys %pavs) {
	my @st=sort {$a<=>$b} keys %{$pavs{$chr}};
	for (my $i=0;$i<@st ;$i++) {
		my $st1=$st[$i];
		my $en1=$pavs{$chr}{$st1};
		$scale{$chr}{$st1}=$en1;
		for (my $j=$i+1;$j<@st ;$j++) {
			my $st2=$st[$j];
			my $en2=$pavs{$chr}{$st2};
		#	print "0\t$st1\t$en1\t$st2\t$en2\n";
			if($st2<=($en1+1)){
				if($en1<=$en2){
		#			print "1\t$st1\t$en1\t$st2\t$en2\n";
					$scale{$chr}{$st1}=$en2;
					$en1=$scale{$chr}{$st1};
		#			print "2\t$st1\t$en1\t$st2\t$en2\n";
					$i=$j;
				}else{
					$i=$j;
				}
			}else{
				$i=$j-1;
				last;
			}
		}
	}
	my $sum=0;
	foreach my $st (sort {$a<=>$b} keys %{$scale{$chr}}) {
		my $len=$scale{$chr}{$st}-$st+1;
		$sum+=$len;
		print OUT "$chr\t$st\t$scale{$chr}{$st}\t$len\t$sum\n";
	}
}
close OUT;



 
