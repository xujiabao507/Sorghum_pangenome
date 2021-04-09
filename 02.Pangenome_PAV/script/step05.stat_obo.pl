use strict;
my $obo=shift;
my $out=shift;

my %match;
my %cont;
my %scale;
my %lens;
	open my $OUT,">$out" or die $!;
	open OBO, $obo =~/gz$/ ? "zcat $obo|" : "<$obo" or die $!;
	<OBO>;<OBO>;<OBO>;<OBO>;<OBO>;
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
		next if($ind<85);
		$lens{$contig}=$a[12];
		$match{$contig}{$chr}+=$percent;
		$cont{$contig}{$chr}{$contig_st}=$_;
		$scale{$contig}{$chr}{$contig_st}=0;
		$scale{$contig}{$chr}{$contig_en}=0;
	}
	close OBO;
	my $sum=0;
	my %contigs;

	my %print;
	foreach my $contig (sort keys %match) {
			my $flag=0;
		foreach my $chr (sort {$match{$contig}{$b}<=>$match{$contig}{$a}} keys %{$match{$contig}}) {
			my @scale=sort {$a<=>$b} %{$scale{$contig}{$chr}};
			$contigs{$contig}=0;
			if($match{$contig}{$chr}>50){
				foreach my $st (sort {$a<=>$b} keys %{$cont{$contig}{$chr}}) {
					print $OUT "$cont{$contig}{$chr}{$st}\t$match{$contig}{$chr}\t$st\t$scale[0]\t$scale[-1]\n";
				}
				$print{$contig}=0;
			}else{
				if(!exists $print{$contig}){
					foreach my $st (sort {$a<=>$b} keys %{$cont{$contig}{$chr}}) {
						print $OUT "$cont{$contig}{$chr}{$st}\t$match{$contig}{$chr}\t$st\t$scale[0]\t$scale[-1]\n";
					}
					$flag=1;
					$print{$contig}=0;
				}else{
					next;
				}
			}
			if($flag==1){last};
		}
	}
	close $OUT;

	my $num_contig=0;
	foreach my $contig (sort keys %contigs) {
		$num_contig++;
		$sum+=$lens{$contig};
	}
	print "$num_contig\t$sum\n";