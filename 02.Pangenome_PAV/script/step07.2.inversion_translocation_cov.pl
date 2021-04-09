use strict;
my $cov=shift;
my $sv_inter=shift;
my $sv_out=shift;
my $out=shift;
	open OUT,">$out" or die $!;
	print "$cov is loading\n";
	my %covs;
	open COV, $cov =~/gz$/ ? "zcat $cov|" : "<$cov" or die $!;
	while (<COV>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		my @a=split/\s+/;
		my $id=$a[0];
		my $contig=$a[1];
		my $st=$a[2];
		my $en=$a[3];
		$covs{$contig}{$st}{$en}=$id;
	}
	close COV;
	print "$sv_inter is loading\n";
	my %svmat;
	my %scale;
	my %type;
	open INTER, $sv_inter =~/gz$/ ? "zcat $sv_inter|" : "<$sv_inter" or die $!;
	while (<INTER>) {
		chomp;
		my @a=split/\s+/;
		my $type1=$a[0];
		my $type2=$a[1];
		my $ctg=$a[2];
		my $st=$a[3];
		my $en=$a[4];
		if($st>$en){
			my $tmp=$en;
			$en=$st;
			$st=$tmp;
		}
		if($type1  eq "TransLoc" or $type2  eq "Inversion" or $type2 eq "TransLoc"){
			$svmat{$ctg}{$st}=$en;
			$type{$ctg}{$st}="IN_SV\t$type1\t$type2";
			$type{$ctg}{$en}="IN_SV\t$type1\t$type2";
		}
		$scale{$ctg}{$st}=$en;
	}
	close INTER;
	print "$sv_out is loading\n";
	open OUTSV, $sv_out =~/gz$/ ? "zcat $sv_out|" : "<$sv_out" or die $!;
	while (<OUTSV>) {
		chomp;
		my @a=split/\s+/;
		my $type1=$a[0];
		my $type2=$a[1];
		my $ctg=$a[2];
		my $st=$a[3];
		my $en=$a[4];
		if($st>$en){
			my $tmp=$en;
			$en=$st;
			$st=$tmp;
		}
		if($type1  eq "TransLoc" or $type2  eq "Inversion" or $type2 eq "TransLoc"){
			$svmat{$ctg}{$st}=$en;
			$type{$ctg}{$st}="OutTL\t$type1\t$type2";
			$type{$ctg}{$en}="OutTL\t$type1\t$type2";
		}
		$scale{$ctg}{$st}=$en;
	}
	close OUTSV;

	print "Ouput\t".localtime()."\n";
	foreach my $ctg (sort keys %scale) {
		my (%up,%down,%edge);
		my @sts=sort {$a<=>$b} keys %{$scale{$ctg}};
		my $num=scalar(@sts)-1;
		for (my $i=0;$i<@sts ;$i++) {
			my $st=$sts[$i];
			my $en=$scale{$ctg}{$st};
			if(exists $svmat{$ctg}{$st}){
				my $h=$i-1;
				my $j=$i+1;
				my ($st1,$en1,$st2,$en2);
				if($h>=0 and $h<=$num){
					$st1=$sts[$h];
					$en1=$scale{$ctg}{$st1};
				}
				if($j>=0 and $h<=$num){
					$st2=$sts[$j];
					$en2=$scale{$ctg}{$st2};
				}
				if($i==0){
					$down{$en}=$st2;
					$edge{$st}{$en}=$st2;
				}
				if($i==$num){
					$up{$st}=$en1;
					$edge{$st}{$st}=$en1;
				}
				if($i>0 and $i<$num){
					$up{$st}=$en1;
					$down{$en}=$st2;
					$edge{$st}{$st}=$en1;
					$edge{$st}{$en}=$st2;
				}
			}
		}
		foreach my $st (sort {$a<=>$b} keys  %{$svmat{$ctg}}) {
			my $en=$svmat{$ctg}{$st};
			my ($up,$down);
			if(!exists $edge{$st}{$st} ){
				$up="NA";
			}else{
				$up=$edge{$st}{$st};
				if($up>$st){
					$up=$st-2000;
				}
			}
			if(!exists $edge{$st}{$en} ){
				$down="NA";
			}else{
				$down=$edge{$st}{$en};
				if($down<$en){
					$down=$en+2000;
				}
			}
			my ($up_read,$down_read,$up_read_high,$down_read_high);
			$up_read=$down_read=$up_read_high=$down_read_high=0;
			print OUT "Summary_Start:\t$ctg\t$st\t$en\t$up\t$down\t$type{$ctg}{$st}\n";
			###################  
			if(exists $edge{$st}{$st} ){
				my $en1=$edge{$st}{$st};
				my $dis=$st-$en1+1;
				my $start=$en1;
				my $end=$st;
				if($start>$end){
					$start=$end-2000;
				}
				foreach my $subst1 (sort {$a<=>$b} keys %{$covs{$ctg}}) {
					foreach my $suben1 (sort {$a<=>$b} keys %{$covs{$ctg}{$subst1}}){
						if($start>=$subst1 and $start<=$suben1 and $end>=$subst1 and $end<=$suben1){
							my $dis1=$start-$subst1+1;
							my $dis2=$suben1-$end+1;
							if($dis1>500 and $dis2>500){
								$up_read_high++;
							}
							print OUT "Up\tCov_High\t$start\t$end\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis1\t$dis2\n";
							$up_read++;
							next;
						}
						if($start>=$subst1 and $start<=$suben1){
							my $dis1=$start-$subst1+1;
							print OUT "Up\tCov_Up\t$start\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis1\n";
							$up_read++;
							next;
						}
						if($end>=$subst1 and $end<=$suben1){
							my $dis2=$suben1-$end+1;
							print OUT "Up\tCov_Down\t$end\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis2\n";
							$up_read++;
							next;
						}
						if($subst1>=$start and $suben1<=$end){
							my $dis1=$start-$subst1;
							my $dis2=$suben1-$end;
							print OUT "Up\tCov_Inter\t$start\t$end\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis1\t$dis2\n";
							$up_read++;
							next;
						}
					}
				}
			}
			if(exists $edge{$st}{$en} ){
				my $st2=$edge{$st}{$en};
				my $dis=$st2-$en+1;
				my $start=$en;
				my $end=$st2;
				if($start>$end){
					$end=$start+2000;
				}
				foreach my $subst1 (sort {$a<=>$b} keys %{$covs{$ctg}}) {
					foreach my $suben1 (sort {$a<=>$b} keys %{$covs{$ctg}{$subst1}}){
						if($start>=$subst1 and $start<=$suben1 and $end>=$subst1 and $end<=$suben1){
							my $dis1=$start-$subst1+1;
							my $dis2=$suben1-$end+1;
							if($dis1>500 and $dis2>500){
								$down_read_high++;
							}
							print OUT "Down\tCov_High\t$start\t$end\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis1\t$dis2\n";
							$down_read++;
							next;
						}
						if($start>=$subst1 and $start<=$suben1){
							my $dis1=$start-$subst1+1;
							print OUT "Down\tCov_Up\t$start\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis1\n";
							$down_read++;
							next;
						}
						if($end>=$subst1 and $end<=$suben1){
							my $dis2=$suben1-$end+1;
							print OUT "Down\tCov_Down\t$end\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis2\n";
							$down_read++;
							next;
						}
						if($subst1>=$start and $suben1<=$end){
							my $dis1=$start-$subst1;
							my $dis2=$suben1-$end;
							print OUT "Down\tCov_Inter\t$start\t$end\t$ctg\t$subst1\t$suben1\t$covs{$ctg}{$subst1}{$suben1}\t$dis1\t$dis2\n";
							$down_read++;
							next;
						}
					}
				}
			}
			if(($up_read+$down_read)>0){
				print OUT "Summary_End:\t$ctg\t$st\t$en\t$up\t$down\t$type{$ctg}{$st}\tRead_Sup:\t$up_read\t$up_read_high\t$down_read\t$down_read_high\n";
			}else{
				print OUT "Summary_End:\t$ctg\t$st\t$en\t$up\t$down\t$type{$ctg}{$st}\tRead_Zero:\t$up_read\t$down_read\n";
			}
		}
	}
	close OUT;
	