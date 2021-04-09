use strict;
my $readmap =shift;
my $order=shift;
my $pav_out=shift;

open OUT,">$pav_out" or die $!;
print "$readmap\n";
my $readnum=0;
my %mapread;
open MR, $readmap =~/gz$/ ? "zcat $readmap|" : "<$readmap" or die $!;
while(<MR>){
	chomp;
	my @a=split/\s+/;
	my $read=$a[0];
	my @b=split/\_/,$read;
	my $chr=$b[0];
	my $st=$b[1];
	my $en=$b[2];
	my $tarchr=$a[2];
	my $tarloc=$a[3];
	$mapread{$chr}{$st}="$tarchr\t$tarloc";
	$readnum++;
	if($readnum%100000==0){
		print "$readnum\t$chr\t$st\t$en\t".localtime()."\n";
	}
}
close MR;


print "$order\n";
my %qst;
my %tst;
my %qst2tst;
open ORDER,$order or die $!;
while (<ORDER>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $qst=$a[1];
	my $qen=$a[2];
	my $tst=$a[3];
	my $ten=$a[4];
	$qst{$chr}{$qst}=$qen;
	$tst{$chr}{$tst}=$ten;
	$qst2tst{$chr}{$qst}=$tst;
}
close ORDER;


print "Transfomration\n";
my %gapread;
foreach my $chr (sort keys %mapread) {
	my @binst=sort  {$a<=>$b} keys %{$qst{$chr}};
	foreach my $qst (sort {$a<=>$b} keys %{$mapread{$chr}}) {
		my $qen=$qst+500;
		my $tag=0;
		for (my $i=0;$i<(scalar(@binst)-1) ;$i++) {
			my $binst1=$binst[$i];
			my $binen1=$qst{$chr}{$binst1};
			if(($qst>=$binst1 and $qst<=$binen1) or ($qen>=$binst1 and $qen<=$binen1)){
				$tag=1;
				last;
			}
		}
		if($tag==0){
			$gapread{$chr}{$qst}=$qen;
		}
	}
}
my %reads;
foreach my $chr (sort keys %gapread) {
	foreach my $loc (sort {$a<=>$b} keys %{$gapread{$chr}}) {
		$qst{$chr}{$loc}=0;
	}
}

print "Output\n";
foreach my $chr (sort keys %qst) {
	foreach my $qst (sort {$a<=>$b} keys %{$qst{$chr}}) {
		if(exists $gapread{$chr}{$qst}){
			print OUT "SV\t$chr\t$qst\t$gapread{$chr}{$qst}\t$mapread{$chr}{$qst}\n";
		}else{
			my $tst=$qst2tst{$chr}{$qst};
			my $ten=$tst{$chr}{$tst};
			print OUT  "Order\t$chr\t$qst\t$qst{$chr}{$qst}\t$tst\t$ten\n";
		}
	}
}
close OUT;






















