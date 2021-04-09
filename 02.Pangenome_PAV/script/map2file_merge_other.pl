use strict;
my $cov=shift;
my $map=shift;
my $mapout=shift;

my %covs;
open COV,$cov or die $!;
while (<COV>) {
	chomp;
	my @a=split/\s+/;
	my $chr1=$a[0];
	#if($chr1 ne "Chr01"){next};
	my $st1=$a[1];
	my $en1=$a[2];
	$covs{$chr1}{$st1}=$en1;
	#print "$chr1\t$st1\t$en1\n";
}
close COV; 

my %maps;
open OUT ,">$mapout" or die $!;
open MAP,"$map" or die $!;
while (<MAP>) {
	chomp;
	my @a=split/\s+/;
	my $chr1=$a[0];
	#if($chr1 ne "Chr01"){next};
	my $st1=$a[1];
	my $en1=$a[2];
	my $chr2=$a[3];
	my $st2=$a[4];
	my $cigar=$a[5];
	$maps{$chr1}{$st1}="$en1\t$chr2\t$st2\t$cigar";
	#print "$chr1\t$st1\t$en1\t$chr2\t$st2\t$cigar\n";
}
close MAP;


foreach my $chr1 (sort keys %maps) {
	my @loc=sort {$a<=>$b} keys %{$covs{$chr1}};
	foreach my $st1 (sort {$a<=>$b} keys %{$maps{$chr1}}) {
		my $en1=(split/\s+/,$maps{$chr1}{$st1})[0];
		my $flag=0;
		for (my $i=0;$i<@loc ;$i++) {
			my $segst=$loc[$i];
			my $segen=$covs{$chr1}{$segst};
			if($segst>$en1){last};
			if($segen<$st1){next};
			#print "1\t$st1\t$en1\t$segst\t$segen\n";
			if(($en1>=$segst and $en1<=$segen) or ($st1>=$segst and $st1<=$segen)){
				print OUT "$chr1\t$st1\t$maps{$chr1}{$st1}\t$segst\t$segen\n";
			#	print "2\t$st1\t$en1\t$segst\t$segen\n";
				$flag=1;
				last;
			}
			
			#if($segst<$st1){next};
		}
		if($flag==0){
			print OUT "$chr1\t$st1\t$maps{$chr1}{$st1}\tOther\n";
		}
	}
}
close OUT;

















