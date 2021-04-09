use strict;
my $fastafile=shift;
my $mapfile=shift;
my $output=shift;
 
open OUT,">$output" or die $!;


my %fastalen;
open CONTIG, $fastafile =~/gz$/ ? "zcat $fastafile|" : "<$fastafile" or die $!;
$/=">";
<CONTIG>;
$/="\n";
while(<CONTIG>){
	chomp;
	my $name=(split/\s+/,$_)[0];
	$/=">";
	my $seq=<CONTIG>;
	$/="\n";
	$seq=~s/\s+//g;
	$seq=~s/\>//g;
	my $len=length($seq);
	$fastalen{$name}=$len;
}
close CONTIG;


open MATCH, $mapfile =~/gz$/ ? "zcat $mapfile|" : "<$mapfile" or die $!;
my %map;
while (<MATCH>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $st=$a[1];
	my $en=$a[2];
	next if($chr=~/sum/);
	$chr=~s/chr\_//g;
	$map{$chr}{$st}=$en;
}
close MATCH;

foreach my $chr (sort keys %map) {
	if(!exists $fastalen{$chr}){
		print "BN:\t$chr\n";
	}
	my $lenchr=$fastalen{$chr};
	my $sumunmap=0;
	my @loc=sort {$a<=>$b} keys %{$map{$chr}};
	my $fst=$loc[0];
	my $stx=1;
	my $enx=$fst-1;
	my $lenx=$enx-$stx+1;
	if($lenx>1){
		print OUT "$chr\t$stx\t$enx\t$lenx\t$lenchr\n";
		$sumunmap+=$lenx;
	}
	for (my $i=0;$i<@loc-1 ;$i++) {
		my $st1=$loc[$i];
		my $en1=$map{$chr}{$st1};
		my $st2=$loc[$i+1];
		my $en2=$map{$chr}{$st2};
		my $sty=$en1+1;
		my $eny=$st2-1;
		my $leny=$eny-$sty+1;
		$sumunmap+=$leny;
		print OUT "$chr\t$sty\t$eny\t$leny\t$lenchr\n";
	}
	my $stz=$loc[-1];
	my $enz=$map{$chr}{$stz};
	$stz=$enz+1;
	my $enz=$lenchr;
	my $lenz=$lenchr-$stz+1;
	if($lenz>0){
		$sumunmap+=$lenz;
		print OUT "$chr\t$stz\t$enz\t$lenz\t$lenchr\n";
	}
	my $rate=$sumunmap/$lenchr;
	print OUT "sum:\t$chr\t$rate\t0\t$sumunmap\t$lenchr\n";
}
foreach my $chr (sort keys %fastalen) {
	if(!exists $map{$chr}){
		print OUT "$chr\t1\t$fastalen{$chr}\t$fastalen{$chr}\t$fastalen{$chr}\n";
		print OUT "sum:\t$chr\t1\t0\t$fastalen{$chr}\t$fastalen{$chr}\n";
	}
}
close OUT;




 


 