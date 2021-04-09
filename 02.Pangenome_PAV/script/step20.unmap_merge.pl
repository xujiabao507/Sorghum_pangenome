use strict;
my $unmapfile=shift;
my $output=shift;
 
open OUT,">$output" or die $!;

open MATCH, $unmapfile =~/gz$/ ? "zcat $unmapfile|" : "<$unmapfile" or die $!;
my %map;
my %lens;
while (<MATCH>) {
	chomp;
	next if($_=~/sum/);
	my @a=split/\s+/;
	my $chr=$a[0];
	my $st=$a[1];
	my $en=$a[2];
	my $len=$a[3];
	next if($len<100);
	#my $lenchr=$a[5];
	$map{$chr}{$st}=$en;
	$lens{$chr}=$len;
}
close MATCH;

foreach my $chr (sort keys %map) {
	my %unmap;
#	next if($chr ne "Contig1036_pilon");
	my @loc=sort {$a<=>$b} keys %{$map{$chr}};
	if(scalar(@loc)<=1){
		my $st1=$loc[0];
		my $en1=$map{$chr}{$st1};
		my $lensum=$en1-$st1+1;
		print OUT "1\t$chr\t100\t$lensum\t$st1\t$en1\t$lensum\n";
		print "$chr\t1\t$lensum\t$st1\t$en1\t$lensum\n";
		next;
	}
	for (my $i=0;$i<@loc ;$i++) {
		my $st1=$loc[$i];
		my $en1=$map{$chr}{$st1};
		$unmap{$st1}{$st1}=0;
		for (my $j=$i+1;$j<@loc ;$j++) {
			my $st2=$loc[$j];
			my $en2=$map{$chr}{$st2};
			if(($st2-$en1)<100000){
				$unmap{$st1}{$st2}=0;
				$unmap{$st1}{$st1}=0;
			}else{
				last;
			}
		}
	}
	my %coverage;
	foreach my $loc (sort {$a<=>$b} keys %unmap) {
		my @loc2=sort {$a<=>$b} keys %{$unmap{$loc}};
		my $st=$loc2[0];
		my $en=$map{$chr}{$loc2[-1]};
		my $lenall=$en-$st+1;
		my $lensum=0;
		foreach my $loc2 (@loc2) {
			my $st2=$loc2;
			my $en2=$map{$chr}{$loc2};
			$lensum+=($en2-$st2+1);
		}
		my $rate=$lensum/$lenall;
		$rate=int($rate*10000+0.5)/100;
		print "-1\t$chr\t$rate\tLoc\t@loc2\tEnd\n";
		if($rate>0.6){
			$coverage{$loc2[0]}=$loc2[-1];
			print "0\t$chr\t$rate\tLoc\t$loc2[0]\t$map{$chr}{$loc2[0]}\t$loc2[-1]\t$map{$chr}{$loc2[-1]}\tlen:$lensum\t$lenall\t$en\n";
		}
	}
	my %merge;
	my @locall=sort {$a<=>$b} keys %coverage;
	for (my $i=0;$i<@locall ;$i++) {
		my $st1=$locall[$i];
		my $stmp=$coverage{$st1};
		my $en1=$map{$chr}{$stmp};
		$merge{$st1}=$en1;
		for (my $j=$i+1;$j<@locall ;$j++) {
			my $st2=$locall[$j];
			my $stmp2=$coverage{$st2};
			my $en2=$map{$chr}{$stmp2};
			print "0.45\t$chr\t$st1\t$en1\t$st2\t$en2\n";
			if($st2<=($en1+1)){
				$merge{$st1}=$en2;
				print "0.5\t$chr\t$st1\t$en1\t$st2\t$en2\n";
				$i=$j
			}else{
				last;
			}
		}
		print "0.75\t$chr\t$st1\t$merge{$st1}\n";
	}
	
	foreach my $stx (sort {$a<=>$b} keys %merge) {
		my $enx=$merge{$stx};
		my $lensum=0;
		foreach my $st2 (sort {$a<=>$b} keys %{$map{$chr}}) {
			my $en2=$map{$chr}{$st2};
			if($st2>=$stx and $en2<=$enx){
				$lensum+=($en2-$st2+1);
			}
		}
		my $lenall=$enx-$stx+1;
		my $rate=$lensum/$lenall;
		$rate=int($rate*10000+0.5)/100;
		print "2\t$chr\t$rate\t$lensum\t$lenall\t$stx\t$enx\n";
		if($lensum>500){
			print OUT "2\t$chr\t$rate\t$lensum\t$lenall\t$stx\t$enx\n";
		}
	}
#	last if($chr  eq "Contig1000_pilon");;
}
close OUT;




 


 