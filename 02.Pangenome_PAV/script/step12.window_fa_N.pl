use strict;

my $window_fa=shift;
my $nseq=shift;

print "$window_fa\n";
my %segfa;
my %segst;
open FA,"zcat $window_fa|" or die $!;
while (my $name=<FA>) {
	my $seq=<FA>;
	chomp $name;
	chomp $seq;
	<FA>;<FA>;
	$name=reverse $name;
	chop $name;
	$name=reverse $name;
	my ($chr,$st,$en)=(split/\_/,$name)[0,1,2];
	$segfa{$chr}{$st}=$seq;
	$segst{$chr}{$st}=$en;
}
close FA;
print "outputing\n";
open SEQ,">$nseq" or die $!;
foreach my $chr (sort keys %segfa) {
	my $ind=0;
	foreach my $st (sort {$a<=>$b} keys %{$segfa{$chr}}) {
		$ind++;
		my $seq=$segfa{$chr}{$st};
		my $len=0;
		if($seq=~/(N+)/){
			$len=length($1);
		}
		print SEQ "$ind\t$chr\t$st\t$segst{$chr}{$st}\t$len\t$seq\n";
	}
}
close SEQ;
























