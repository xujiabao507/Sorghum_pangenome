use strict;
my $obo=shift;
my $delta=shift;
my $out=shift;

	my %match;
	open my $OUT,">$out" or die $!;
	open OBO, $obo =~/gz$/ ? "zcat $obo|" : "<$obo" or die $!;
	while (<OBO>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		my @a=split/\s+/;
		my $chr=$a[19];
		my $contig=$a[20];
		my $contig_st=$a[3];
		my $contig_en=$a[4];
		my $chr_st=$a[0];
		my $chr_en=$a[1];
		$match{"$contig\t$contig_st\t$contig_en"}{"$chr\t$chr_st\t$chr_en"}=0;
	}
	close OBO;
	
	open DELTA, $delta =~/gz$/ ? "zcat $delta|" : "<$delta" or die $!;
	
	my $firsthead1=<DELTA>;
	my $firsthead2=<DELTA>;
	print $OUT "$firsthead1$firsthead2";
	$/=">";
	while (my $line=<DELTA>) {
		$/="\n";
		$line=~s/>//g;
		my @line=split/\n+/,$line;
		my $head=$line[0];
		my ($tar,$query)=(split/\s+/,$head);
		my $strout;
		for (my $i=1;$i<@line ;$i++) {
			my @str=split/\s+/,$line[$i];
			if(scalar(@str)>1){
				my $tar_loc="$str[0]\t$str[1]";
				my $query_loc="$str[2]\t$str[3]";
				if(exists $match{"$query\t$query_loc"}{"$tar\t$tar_loc"}){
					#print "$tar\t$tar_loc\t$query\t$query_loc\n";
					$strout.="$line[$i]\n";
					for (my $j=$i+1;$j<@line ;$j++) {
						if($line[$j]!=0){
							$strout.="$line[$j]\n";
						}else{
							$strout.="$line[$j]\n";
							$i=$j;
							last;
						}
					}
				}
			}
		}
		if($strout=~/\w/){
			print $OUT ">$head\n$strout";
		}
		$/=">";
	}
	close DELTA;
	close $OUT;

	
	