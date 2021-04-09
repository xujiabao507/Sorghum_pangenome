my $mapfile=shift;
my $gapout=shift;
open OUT,">$gapout" or die $!;
my %scale;
open MAP,$mapfile or die $!;
while (<MAP>) {
	chomp;
	next if($_!~/\w/);
	my @a=split/\s+/;
	my $chr1=$a[0];
#	next if($chr1 ne "Chr08");
	my $st1=$a[1];
	my $en1=$a[2];
	my $min=$st1;
	my $max=$en1;
	if($st1>$en1){
		$min=$en1;
		$max=$st1;
	}
	$max=$max+500;
	$scale{$chr1}{$min}=$max;
	#print "0\t$chr1\t$min\t$max\n";
}
close MAP;

foreach my $chr (sort keys %scale) {
	my %scalemerge;
	my @loc=sort {$a<=>$b} keys %{$scale{$chr}};
	foreach my $st (@loc) {
		my $en=$scale{$chr}{$st};
		#print "1\t$chr\t$st\t$en\n";
	}
	my $num=scalar(@loc);
	for (my $i=0;$i<@loc ;$i++) {
		my $st1=$loc[$i];
		my $en1=$scale{$chr}{$st1};
		$scalemerge{$st1}=$en1;
		for (my $j=$i+1;$j<@loc;$j++) {
			my $st2=$loc[$j];
			my $en2=$scale{$chr}{$st2};
			print "1\t$num\t$i\t$j\t$chr\t$st1\t$scalemerge{$st1}\t$st2\t$en2\n";
			my $dis=$st2-$scalemerge{$st1};
			if($dis<0){
				if($scalemerge{$st1}<$en2){
					$scalemerge{$st1}=$en2;
					print "2\t$num\t$i\t$j\t$chr\t$st1\t$scalemerge{$st1}\t$st2\t$en2\t$dis\n";
					print "2\t$num\t$i\t$j\t$chr\t$st1\t$scalemerge{$st1}\n";
				}else{
					print "3\t$num\t$i\t$j\t$chr\t$st1\t$scalemerge{$st1}\t$st2\t$en2\t$dis\n";
					print "3\t$num\t$i\t$j\t$chr\t$st1\t$scalemerge{$st1}\n";
				}
				if($j==(scalar(@loc)-1)){
					$i=$j;
				}
				#print "2\t$i\t$j\n";
			}else{
				$i=$j-1;
				print "4\t$num\t$i\t$j\t$chr\t$st1\t$en1\t$st2\t$en2\t$dis\n";
				#print "3\t$i\t$j\n";
				last;
			}
		}
	}
	my @loc2=sort {$a<=>$b} keys %scalemerge;
	for (my $i=0;$i<(scalar(@loc2)-2) ;$i++) {
		my $loc=$loc2[$i];
		my $dis=$loc2[$i+1]-$scalemerge{$loc};
		print OUT "$chr\t$scalemerge{$loc}\t$loc2[$i+1]\t$dis\n";
	}

	foreach my $st (sort {$a<=>$b} keys %scalemerge) {
		my $dis=$scalemerge{$st}-$st;
		#print "3\t$st\t$scalemerge{$st}\t$dis\n";
	}
}

close OUT;













