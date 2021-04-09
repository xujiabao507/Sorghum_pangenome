use strict;
my $blast_map=shift;
my $output=shift;


my $query_tmp;
my $num=0;
my %coverread;
my %contigs;
open my $OUT,">$output" or die $!;
open my $BLAST,$blast_map=~/gz$/?"zcat $blast_map|":"<$blast_map" or die $!;
while (<$BLAST>) {
	chomp;
	my @a=split/\s+/;
	my $query=$a[0];
	my $chr=$a[1];
	my $ind=$a[2];
	if($ind<90){next};
	my $query_st=$a[6];
	my $query_en=$a[7];
	my $tar_st=$a[8];
	my $tar_en=$a[9];
	my $score=$a[11];
	my $len_query=abs($query_st-$query_en)+1;
	next if($len_query<100);
	if($num==0){$query_tmp=$query};
	if(!exists $contigs{$query}){
		if($num>0){
			print "$num\t1.1\t$query_tmp\t1.2\t$query\n";
			my %readscov;
			my %locs;
			my @len=sort {$b<=>$a} keys %coverread;
			########## 初始化
			my $len=$len[0];
			foreach my $str (sort keys %{$coverread{$len}}) {
				#print $OUT  "0\t$str\n";
				my ($qst,$qen)=(split/\s+/,$str)[6,7];
				$readscov{$qst}{$str}=0;
				for (my $loc=$qst;$loc<=$qen ;$loc++) {
					$locs{$loc}{"$qst\t$qen"}=0;
				}
			}
			########## next length
			for (my $i=1;$i<@len ;$i++) {
				$len=$len[$i];
				foreach my $str (sort  keys %{$coverread{$len}}) { ### 该长度的every line
					my ($qst,$qen)=(split/\s+/,$str)[6,7];
					my %lencover;
					for (my $loc=$qst;$loc<=$qen ;$loc++) {
						if(exists $locs{$loc}){
							foreach my $loc_tar (sort keys %{$locs{$loc}}) {
								$lencover{$loc_tar}++;
							}
						}
					}
					my @scf=sort keys %lencover;
					if(@scf==0){
						#print $OUT  "1\t$str\t0\n";
						$readscov{$qst}{$str}=0;
						for (my $loc=$qst;$loc<=$qen ;$loc++) {
							$locs{$loc}{"$qst\t$qen"}=0;
						}
						next;
					}else{
						my $flag=0;
						foreach my $scf (sort {$lencover{$b}<=>$lencover{$a}} keys %lencover) {
							my ($sta,$ena)=(split/\s+/,$scf)[0,1];
							my $rate1=$lencover{$scf}/(abs($sta-$ena)+1);
							my $rate2=$lencover{$scf}/$len;
							if($rate1<0.2 && $rate2<0.2){
								$readscov{$qst}{$str}=0;
								#print $OUT "2\t$str\n";
								#print "2\t$qst\t$qen\t$lencover{$scf}\t$rate1\t$rate2\t$lencover{$scf}\t$scf\n";
								$flag=1;
							}
							last;
						}
						if($flag==0){next};
						for (my $loc=$qst;$loc<=$qen ;$loc++) {
							$locs{$loc}{"$qst\t$qen"}=0;
							#print $OUT "3\t$loc\t$qst\t$qen\n";
						}
					}
				}
				#my @loc=sort keys %locs;
				#my $locnum=scalar(@loc);
				#print "4\t$locnum\n";
			}
			foreach my $loc (sort {$a<=>$b} keys %readscov) {
				foreach my $str (sort keys %{$readscov{$loc}}) {
					print $OUT "$str\n";
				}
			}
		}
		%coverread=();
		$query_tmp=$query;
		my $len=abs($query_en-$query_st)+1;
		$contigs{$query}=0;
		$coverread{$len}{$_}=0;
		#last;
	}else{
		my $len=abs($query_en-$query_st)+1;
		$contigs{$query}=0;
		$coverread{$len}{$_}=0;
	}
	$num++;
	#if($query ne "Contig1122_pilon_total_2886851"){last};
}
		if($num>0){
			print "the last one\n";
			my %readscov;
			my %locs;
			my @len=sort {$b<=>$a} keys %coverread;
			########## 初始化
			my $len=$len[0];
			foreach my $str (sort keys %{$coverread{$len}}) {
				#print $OUT  "0\t$str\n";
				my ($qst,$qen)=(split/\s+/,$str)[6,7];
				$readscov{$qst}{$str}=0;
				for (my $loc=$qst;$loc<=$qen ;$loc++) {
					$locs{$loc}{"$qst\t$qen"}=0;
				}
			}
			########## next length
			for (my $i=1;$i<@len ;$i++) {
				$len=$len[$i];
				foreach my $str (sort  keys %{$coverread{$len}}) { ### 该长度的every line
					my ($qst,$qen)=(split/\s+/,$str)[6,7];
					my %lencover;
					for (my $loc=$qst;$loc<=$qen ;$loc++) {
						if(exists $locs{$loc}){
							foreach my $loc_tar (sort keys %{$locs{$loc}}) {
								$lencover{$loc_tar}++;
							}
						}
					}
					my @scf=sort keys %lencover;
					if(@scf==0){
						#print $OUT  "1\t$str\t0\n";
						$readscov{$qst}{$str}=0;
						for (my $loc=$qst;$loc<=$qen ;$loc++) {
							$locs{$loc}{"$qst\t$qen"}=0;
						}
						next;
					}else{
						my $flag=0;
						foreach my $scf (sort {$lencover{$b}<=>$lencover{$a}} keys %lencover) {
							my ($sta,$ena)=(split/\s+/,$scf)[0,1];
							my $rate1=$lencover{$scf}/(abs($sta-$ena)+1);
							my $rate2=$lencover{$scf}/$len;
							if($rate1<0.2 && $rate2<0.2){
								$readscov{$qst}{$str}=0;
								#print $OUT "2\t$str\n";
								#print "2\t$qst\t$qen\t$lencover{$scf}\t$rate1\t$rate2\t$lencover{$scf}\t$scf\n";
								$flag=1;
							}
							last;
						}
						if($flag==0){next};
						for (my $loc=$qst;$loc<=$qen ;$loc++) {
							$locs{$loc}{"$qst\t$qen"}=0;
							#print $OUT "3\t$loc\t$qst\t$qen\n";
						}
					}
				}
				#my @loc=sort keys %locs;
				#my $locnum=scalar(@loc);
				#print "4\t$locnum\n";
			}
			foreach my $loc (sort {$a<=>$b} keys %readscov) {
				foreach my $str (sort keys %{$readscov{$loc}}) {
					print $OUT "$str\n";
				}
			}
		}
close $BLAST;
close $OUT;
















