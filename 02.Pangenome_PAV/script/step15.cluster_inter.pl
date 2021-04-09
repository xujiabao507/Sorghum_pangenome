use strict;

my $map=shift;
my $out=shift;
my $out_contig=shift;
open OUT,">$out" or die  $!;
open OUT2,">$out_contig" or die  $!;
print "$map start\n";
my %clu_mainchr;
my %tar;
my %maps;
open MAP,$map or die $!;
while (<MAP>) {
	chomp;
	if($_=~/NA/){next};
	my @a=split/\s+/;
	my $chr1=$a[0];
	#if($chr1 ne "Chr01"){last};
	my $st1=$a[1];
	my $ctnum=$a[2];
	if($ctnum eq "Gap"){next};
	my $chrt=$a[4];
	my $stt=$a[5];
	$clu_mainchr{$chr1}{$ctnum}{$chrt}++;
	$tar{$chr1}{$ctnum}{$st1}=0;
	$maps{$chr1}{$st1}{$chrt}=$stt;
}
close MAP;

my %cluster2chr;
my %cluster2sum;
foreach my $chr1 (sort keys %clu_mainchr) {
	foreach my $ctnum (sort {$a<=>$b} keys %{$clu_mainchr{$chr1}}) {
		my @chr=sort {$clu_mainchr{$chr1}{$ctnum}{$b}<=>$clu_mainchr{$chr1}{$ctnum}{$a}} keys %{$clu_mainchr{$chr1}{$ctnum}};
		my $sum=0;
		foreach my $chr (@chr) {
			$sum+=$clu_mainchr{$chr1}{$ctnum}{$chr};
		}
		my $rate=int(($clu_mainchr{$chr1}{$ctnum}{$chr[0]}/$sum)*100+0.5)/100;
		$cluster2chr{$chr1}{$ctnum}="$chr[0]\t$rate\t$clu_mainchr{$chr1}{$ctnum}{$chr[0]}\t$sum";
	}
}
my %scale;
foreach my $chr1 (sort keys %tar) {
	foreach my $ctnum (sort {$a<=>$b} keys %{$tar{$chr1}}) {
		my @loc=sort {$a<=>$b} keys %{$tar{$chr1}{$ctnum}};
		print "$chr1\t$ctnum\t@loc\n";
		my ($chr2,$rate,$num,$sum)=(split/\s+/,$cluster2chr{$chr1}{$ctnum})[0,1,2,3];
		my %clusterscale;
		foreach my $loc (@loc) {
			if(exists $maps{$chr1}{$loc}{$chr2}){
				$clusterscale{$loc}=$maps{$chr1}{$loc}{$chr2};
			}
		}
		my @loc2=sort {$a<=>$b} keys %clusterscale;
		my $len=abs($loc[-1]-$loc[0])+501;
		#################### 出现最多的 100k区域
		my %main_100k;
		my %main_100k2site;
		foreach my $loc (@loc2) {
			my $k100=int($clusterscale{$loc}/100000)*100000;
			$main_100k{$k100}++;
			$main_100k2site{$k100}{$clusterscale{$loc}}=0;
			#print "-1\t$loc\t$clusterscale{$loc}\n";
		}
		my @k100=sort {$main_100k{$b}<=>$main_100k{$a}} keys %main_100k;
		my $kmost=$k100[0];
		my $sum=0;
		my @loctar=sort {$a<=>$b} keys %{$main_100k2site{$kmost}};
		foreach my $loctar (@loctar) {
			$sum+=$loctar;
		}
		my $mean=int($sum/scalar(@loctar));
		my @loc3;
		foreach my $loc (@loc2) {
			if(abs($maps{$chr1}{$loc}{$chr2}-$mean)<$len){
				push @loc3,$loc;
			}
		}
		@loc2=();
		@loc2=@loc3;
		####################
		my ($locmin,$min);
		$locmin=$loc2[0];
		$min=$clusterscale{$locmin};
		my ($locmax,$max);
		$locmax=$loc2[-1];
		$max=$clusterscale{$locmax};
		my $len2=abs($locmax-$locmin)+500;
		my $len3=abs($max-$min)+500;
		print OUT2 "$chr1\t$ctnum\t$len\t$len2\t$len3\t$locmin\t$locmax\t$chr2\t$min\t$max\n";
		#last;
		if($min=~/\d/ and $max=~/\d/){
			if($min<$max){
				$scale{$chr2}{$min}=$max+500;
			}else{
				$scale{$chr2}{$max}=$min+500;
			}
		}
		#if($min =~ /\d/){
			#$scale{$chr1}{$ctnum}{0}="$chr2\t$min";
			#$scale{$chr1}{$ctnum}{1}="$chr2\t$max";
		#	$scale{$chr2}{$min}=$max+500;
		#}
	}
	#last;
}

foreach my $chr (sort keys %scale) {
	my %scalemerge;
	my @loc=sort {$a<=>$b} keys %{$scale{$chr}};
	for (my $i=0;$i<@loc ;$i++) {
		my $loc=$loc[$i];
		print OUT "$chr\t$loc\t$scale{$chr}{$loc}\n";
	}
#	for (my $i=0;$i<@loc ;$i++) {
#		my $st1=$loc[$i];
#		my $en1=$scale{$chr}{$st1};
#		$scalemerge{$st1}=$en1;
#		for (my $j=$i+1;$j<@loc;$j++) {
#			my $st2=$loc[$j];
#			my $en2=$scale{$chr}{$st2};
#			if($st2<=$en1){
#				if($en1<$en2){
#					$scalemerge{$st1}=$en2;
#					print "1\t$chr\t$st1\t$en1\t$st2\t$en2\n";
#				}else{
#					print "2\t$chr\t$st1\t$en1\t$st2\t$en2\n";
#					next;
#				}
#				#print "$st1\t$en1\t$st2\t$en2\n";
#			}else{
#				$i=$j-1;
#				#$scalemerge{$st2}=$en2;
#				last;
#			}
#		}
#	#	print OUT "$chr\t$scale{$chr}{$loc}\t$loc[$i+1]\n";
#	}
#	my @locs2=sort {$a<=>$b} keys %scalemerge;
#	for (my $i=0;$i<scalar(@locs2)-1;$i++) {
#		my $loc=$locs2[$i];
#		print OUT "$chr\t$scalemerge{$loc}\t$locs2[$i+1]\n";
#	}
}
 
close OUT;
close OUT2;






