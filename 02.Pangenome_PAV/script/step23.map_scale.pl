use strict;
my $fastafile=shift;
my $mapfile=shift;
my $output=shift;
 

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


open OUT,">$output" or die $!;
open MATCH, $mapfile =~/gz$/ ? "zcat $mapfile|" : "<$mapfile" or die $!;
my %scale;
my %covs;
while (<MATCH>) {
	chomp;
	my @a=split/\s+/;
	my $contig=$a[0];
	if($contig=~/chr_/){
		$contig=~s/chr_//g;
	}
	my $st=$a[1];
	if($st==0){$st=1};
	my $en=$a[2];
	my $rate=$a[3];
	my $cov=$a[4];
	if($rate>0.6){
		$scale{$contig}{$st}=$en;
		$covs{$contig}{$st}=$cov;
	}
}
close MATCH;

foreach my $contig (sort keys %scale) {
#	if($contig ne "Contig1030_pilon"){next};
	my @loc=sort {$a<=>$b} keys %{$scale{$contig}};
	my %locs;
	for (my $i=0;$i<@loc ;$i++) {
		my $st1=$loc[$i];
		my $en1=$scale{$contig}{$st1};
		my $enx=$en1;
		$locs{$st1}=$enx;
		for (my $j=$i+1;$j<@loc ;$j++) {
			my $st2=$loc[$j];
			my $en2=$scale{$contig}{$st2};
			if($st2<=($enx+1)){
				$enx=$en2;
				$locs{$st1}=$enx;
				$i=$j;
			}else{
				#$locs{$st1}=$enx;
				#print "0\t$contig\t$i\t$st1\t$enx\n";
				last;
			}
		}
		#$locs{$st1}=$enx;
		#print "1\t$contig\t$i\t$st1\t$enx\n";
	}
	#print "2\n";
	my $lensum=0;
	foreach my $st (sort {$a<=>$b} keys %locs) {
		my $len=$locs{$st}-$st+1;
		$lensum+=$len;
		#print OUT "$contig\t$st\t$locs{$st}\t$len\t$fastalen{$contig}\n";
		my $sumbp=0;
		my $covsum=0;
		for (my $loc=$st;$loc<=$locs{$st} ;$loc+=100) {
			if(exists $covs{$contig}{$loc}){
				$covsum+=$covs{$contig}{$loc};
				$sumbp+=($scale{$contig}{$loc}-$loc+1);
				#print "$contig\t$loc\t$scale{$contig}{$loc}\n";
			}
		}
		my $covrate=$covsum/$sumbp;
		#print "3\t$contig\t$st\t$locs{$st}\n";
		print OUT "$contig\t$st\t$locs{$st}\t$len\t$covrate\t$covsum\t$sumbp\t$fastalen{$contig}\n";
	}
	my $rate=$lensum/$fastalen{$contig};
	print OUT "sum:\t$contig\t$lensum\t$fastalen{$contig}\t$rate\n";
}
close OUT;




 


 
