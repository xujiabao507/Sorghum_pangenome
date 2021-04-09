use strict;

my $mapread=shift;
my $out=shift;
open OUT,">$out" or die $!;
my %qst;
my %gaps;
open MR,$mapread or die $!;
while (<MR>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[1];
	my $st=$a[2];
	my $seq=$a[5];
	my $num=0;
	while ($seq=~/(N+)/) {
		my $seq1=$`;
		my $seq2=$1;
		my $seq3=$';
		$num+=length($seq2);
		$seq=$seq3;
	}
	$qst{$chr}{$st}=$num;
	if($num>50){
		$gaps{$chr}{$st}="Gap";
	}
}
close MR;

foreach my $chr (sort keys %qst) {
	my $ind=1;
	my @st=sort {$a<=>$b} keys %{$qst{$chr}};
	for (my $i=0;$i<scalar(@st)-1;$i++) {
		my $loc1=$st[$i];
		my $loc2=$st[$i+1];
		if(exists $gaps{$chr}{$loc1} and !exists $gaps{$chr}{$loc2}  ){
			$ind++;
		}
		if(!exists $gaps{$chr}{$loc1}){
			print OUT "$chr\t$loc1\t$ind\t$qst{$chr}{$loc1}\n";
		}else{
			print OUT "$chr\t$loc1\tGap\t$qst{$chr}{$loc1}\n";
		}
	}
	print OUT "$chr\t$st[-1]\t$ind\t$qst{$chr}{$st[-1]}\n";
}
close OUT;


















