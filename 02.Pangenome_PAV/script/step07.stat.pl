use strict;
my $obo=shift;
my $out=shift;

my %match;
my %lens;
	open my $OUT,">$out" or die $!;
	open OBO, $obo =~/gz$/ ? "zcat $obo|" : "<$obo" or die $!;
	while (<OBO>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		my @a=split/\s+/;
		my $chr=$a[19];
		my $contig=$a[20];
		my $percent=$a[15];
		my $contig_st=$a[3];
		my $contig_en=$a[4];
		my $ind=$a[9];
		$lens{$contig}=$a[12];
		$match{$contig}{$chr}+=$percent;
	}
	close OBO;
	my $sum_seq=0;
	my %perc;
	my %num;
	my %perc_contig;
	my %num_contig;
	my $sum=0;
	foreach my $contig (sort keys %match) {
		my $num=scalar(keys %{$match{$contig}});
		$sum_seq+=$lens{$contig};
		$num{$num}++;
		$num_contig{$num}{$contig}=0;
		my $max=0;
		$sum++;
		foreach my $chr(sort keys %{$match{$contig}}) {
			if($max<$match{$contig}{$chr}){
				$max=$match{$contig}{$chr};
			}
		}
		my $int=int($max/5);
		$perc{$int}++;
		$perc_contig{$int}{$contig}=0;
	}

	
	print $OUT "Total seq:\t$sum\t$sum_seq\n";
	my $sum11=0;
	my $sum12=0;
	foreach my $num (sort {$a<=>$b} keys %num) {
		my $rate1=$num{$num}/$sum;
		my $numseq_sum=0;
		foreach my $contig (sort keys %{$num_contig{$num}}) {
			$numseq_sum+=$lens{$contig};
		}
		my $rate2=$numseq_sum/$sum_seq;
		$sum11+=$rate1;
		$sum12+=$rate2;
		my $sum111=1-$sum11;
		my $sum122=1-$sum12;
		print $OUT "Number\t$num\t$rate1\t$sum111\t$num{$num}\t$sum\t$numseq_sum\t$rate2\t$sum122\n";
	}
	$sum11=0;
	$sum12=0;
	foreach my $perct (sort {$a<=>$b} keys %perc) {
		my $rate1=$perc{$perct}/$sum;
		my $numseq_sum=0;
		foreach my $contig (sort keys %{$perc_contig{$perct}}) {
			$numseq_sum+=$lens{$contig};
		}
		my $perctage=$perct*5;
		my $rate2=$numseq_sum/$sum_seq;
		$sum11+=$rate1;
		$sum12+=$rate2;
		my $sum111=1-$sum11;
		my $sum122=1-$sum12;
		print $OUT "Perc\t$perctage\t$rate1\t$sum111\t$perc{$perct}\t$sum\t$numseq_sum\t$rate2\t$sum122\n";
	}
	close $OUT;