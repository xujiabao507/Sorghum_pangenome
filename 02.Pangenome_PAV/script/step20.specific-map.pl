use strict;
my $contigfile=shift;
my $mapfile=shift;
my $output=shift;

open my $OUT,">$output" or die $!;

my %queryfa;
my %querylen;
open CONTIG, $contigfile =~/gz$/ ? "zcat $contigfile|" : "<$contigfile" or die $!;
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
	$querylen{$name}=$len;
	$queryfa{$name}=$seq;
}
close CONTIG;


my %map1;
open MAP, $mapfile =~/gz$/ ? "zcat $mapfile|" : "<$mapfile" or die $!;
while (<MAP>) {
	chomp;
	my @a=split/\s+/;
	my $query=$a[0];
	my $cst=$a[1];
	my $cen=$a[2];
	my $chr=$a[4];
	my $chst=$a[5];
	my $chen=$a[6];

	################
	my $mmin=$chst;
	my $mmax=$chen;
	if($mmin>$chen){
		$mmin=$chen;
		$mmax=$chst;
	}
	$map1{$chr}{$mmin}{$mmax}=0;
}
close MAP;

my %map;
my $sum1=0;
foreach my $chr (sort keys %map1) {
	foreach my $st (sort {$a<=>$b} keys %{$map1{$chr}}) {
		my @loc=sort {$b<=>$a} keys %{$map1{$chr}{$st}};
		$map{$chr}{$st}=$loc[0];
	#	$sum1+=($loc[0]-$st+1);
	#	my $sub=$loc[0]-$st+1;
	#	print "0\tsum\t$chr\t$sum1\t$sub\t$st\t@loc\n";
		#if(scalar(@loc)>1){
		#	print "0\t$chr\t$st\t@loc\n";
		#}
	}
}

foreach my $chr (sort keys %map) {
	my @st=sort {$a<=>$b} keys %{$map{$chr}};
	my %start;
	my $len=$querylen{$chr};
	for (my $i=0;$i<@st ;$i++) {
		if($st[$i]==0){next};
		$start{$st[$i]}=0;
		my $st1=$st[$i];
		my $en1=$map{$chr}{$st1};
		for (my $j=$i+1;$j<@st ;$j++) {
			my $st2=$st[$j];
			my $en2=$map{$chr}{$st2};
			if($st1<$st2 and $en1>$en2){
				$st[$j]=0;
			}
		}
	}
	@st=();
	@st=sort {$a<=>$b} keys %start;
	my $lensum=0;
	# my $lensum2=0;
	for (my $i=0;$i<@st-1 ;$i++) {
		my $st1=$st[$i];
		my $en1=$map{$chr}{$st1};
		my $st2=$st[$i+1];
		my $en2=$map{$chr}{$st2};
		my $enx=$en1;
		if($st2<$en1){
			$enx=$st2;
		};
		$lensum+=($enx-$st1+1);
	#	$lensum2+=($en1-$st1+1);
	#	my $dis=$st2-$en1-1;
	#	print "1\t$chr\t$st1\t$en1\t$st2\t$en2\t$lensum\t$lensum2\n";
	#	print "0\t$st1\t$en1\t$st2\t$en2\t$dis\n";
	}
	my $rate=$lensum/$len;
####################################
	my $st=$st[0];
	if(($st-1)>500){
		my $st=$st-1;
		print $OUT "$chr\t1\t$st\t$st\t$lensum\t$len\t$rate\n";
	}
	for (my $i=0;$i<@st-1 ;$i++) {
		my $st1=$st[$i];
		my $en1=$map{$chr}{$st1};
		my $st2=$st[$i+1];
		my $en2=$map{$chr}{$st2};
		my $lenpart=($st2-1)-($en1+1)+1;
		
		if($lenpart>500){
			my $loc1=$en1+1;
			my $loc2=$st2-1;
			print $OUT "$chr\t$loc1\t$loc2\t$lenpart\t$lensum\t$len\t$rate\n";
		}
	}
	my $en=$map{$chr}{$st[-1]};
	
	my $loc1=$en+1;
	my $lenpart=$len-($en+1)+1;
	if($lenpart>500){
		print $OUT "$chr\t$loc1\t$len\t$lenpart\t$lensum\t$len\t$rate\n";
	}
}

##################################################   没比对上的
foreach my $contig (sort keys %querylen) {
	if(!exists $map{$contig}){
		my $lensub=$querylen{$contig};
		print $OUT "$contig\t1\t$lensub\t$lensub\t0\t0\t0\n";
	}
}

close $OUT;






	 