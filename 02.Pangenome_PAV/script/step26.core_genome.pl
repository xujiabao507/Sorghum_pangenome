use strict;

 
open OUT,">$ARGV[-1].CG.depth" or die $!;
my %lens;
open LEN,$ARGV[0] =~/gz$/ ? "zcat $ARGV[0]|" : "<$ARGV[0]" or die $!;
while (<LEN>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $len=$a[1];
	$lens{$chr}=$len;
}
close LEN;

my %map;
for (my $i=1;$i<(@ARGV-1) ;$i++) {
	my $mapfile=$ARGV[$i];
	my $base=(split/\//,$mapfile)[-1];
	my $sam=(split/\./,$base)[0];
	open MATCH, $mapfile =~/gz$/ ? "zcat $mapfile|" : "<$mapfile" or die $!;
	while (<MATCH>) {
		chomp;
		my @a=split/\s+/;
		my $chr=$a[0];
		my $st=$a[1];
		my $en=$a[2];
		$map{$chr}{$st}{$en}{$sam}=0;
	}
	close MATCH;
}

foreach my $chr (sort keys %map) {
	my $interval=1000000;
	for (my $i=1;$i*$interval<=$lens{$chr} ;$i++) {
		my $st=($i-1)*$interval+1;
		my $en=$i*$interval;
		my %depth;
		foreach my $st1(sort {$a<=>$b} keys %{$map{$chr}}) {
			foreach my $en1 (sort {$a<=>$b}keys %{$map{$chr}{$st1}}) {
				if($en1<$st){next};
				if($st1>$en){last};
				my $flag=0;
				if($en1>=$st and $en1<=$en or $st1>=$st and $st1<=$en1){
					$flag=1;
				}
				if($flag==1){
					foreach my $sam (sort keys %{$map{$chr}{$st1}{$en1}}) {
						for (my $loc=$st1;$loc<=$en1 ;$loc++) {
							$depth{$loc}++;
						}
					}
				}
			}
		}
		for (my $pos=$st;$pos<=$en ;$pos++) {
			if(exists $depth{$pos}){
				print OUT "$chr\t$pos\t$depth{$pos}\n";
			}else{
				print OUT "$chr\t$pos\t0\n";
			}
		}
	}
}

close OUT;









