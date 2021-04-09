use strict;
my $obo=shift;
my $out=shift;

	my %match;
	my %match_rev;
	my %chr_reg_pre;
	my %contig_len;
	
	open my $OUT1,">$out.TransLocation" or die $!;
	open my $OUT2,">$out.Inversion" or die $!;
	open OBO, $obo =~/gz$/ ? "zcat $obo|" : "<$obo" or die $!;
	while (<OBO>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		my @a=split/\s+/;
		my $chr=$a[19];
		my $contig=$a[20];
		my $rate=$a[21];
		my $contiglen=$a[24];
		next if($rate<80);
		my $chrst=$a[0];
		my $chren=$a[1];
		my $contigst=$a[3];
		my $contigen=$a[4];
		$contig_len{$contig}=$contiglen;
		$match{$contig}{$contigst}="$chr\t$chrst";
		$match{$contig}{$contigen}="$chr\t$chren";
		$match_rev{$chr}{$chrst}{$contig}="$contigst";
		$match_rev{$chr}{$chren}{$contig}="$contigen";
		$chr_reg_pre{$contig}{$contigst}=$contigen;
		#$chr_reg{$contig}{$chr}{$chrst}=$chren;
	}
	close OBO;
	my %contig_reg;
	my %chr_reg;
	foreach my $contig (sort keys %chr_reg_pre) {
		my @st=sort {$a<=>$b} keys %{$chr_reg_pre{$contig}};
		for (my $i=0;$i<@st ;$i++) {
			my $flag=0;
			my $st1=$st[$i];
			my $en1=$chr_reg_pre{$contig}{$st1};
			my $min1=$st1;
			my $max1=$en1;
			if($st1>$en1){
				$min1=$en1;$max1=$st1;
			}
			for (my $j=0;$j<@st ;$j++) {
				next if($i==$j);
				my $st2=$st[$j];
				my $en2=$chr_reg_pre{$contig}{$st2};
				if($st1==$st2 and $en1==$en2){next};
				my $min2=$st1;
				my $max2=$en1;
				if($st2>$en2){
					$min2=$en2;$max2=$st2;
				}
				if($st1>=$st2 and $en1<=$en2){
					$flag=1;
					last;
				}
			}
			if($flag==0){
				$contig_reg{$contig}{$st1}=$en1;
				my ($chr,$chrst)=(split/\s/,$match{$contig}{$st1})[0,1];
				my $chren=(split/\s/,$match{$contig}{$en1})[1];
				$chr_reg{$contig}{$chr}{$chrst}=$chren;
				#if($contig eq "Contig1000_pilon"){
				#	print "$contig\t$st1\t$en1\n"
				#}
			}
		}
	}

	foreach my $contig (sort keys %contig_reg) {
		my %reg;
		my $chr_num=scalar(keys %{$chr_reg{$contig}});
		if($chr_num>1){next};
		my $chrtar=0;
		foreach my $chr (sort keys %{$chr_reg{$contig}}) {
			$chrtar=$chr;
			foreach my $st (sort {$a<=>$b} keys %{$chr_reg{$contig}{$chr}} ) {
				my $en=$chr_reg{$contig}{$chr}{$st};
				$reg{$st}=$en;
				#if($contig eq "Contig1000_pilon"){
				#	print "2\t$contig\t$st\t$en\n";
				#}
			}
		}
		my $len=$contig_len{$contig};
		my @st=sort {$a<=>$b} keys %reg;
		if(scalar(@st)==1){next};
		my %cluster;
		my %cluster_cov;
		my $sttmp=$st[0];
		my $entmp=$reg{$sttmp};
		my $cont_st=$match_rev{$chrtar}{$sttmp}{$contig};
		my $cont_en=$match_rev{$chrtar}{$entmp}{$contig};
		my $seglen=abs($cont_en-$cont_st)+1;
		#print "3.5\t$contig\t$chrtar\t$sttmp\t$entmp\t$len\n";
		$cluster{$sttmp}{$sttmp}=0;
		$cluster_cov{$sttmp}+=$seglen;
		for (my $i=1;$i<@st ;$i++) {
			my $st=$st[$i];
			my $en=$reg{$st};
			my $dis=$st-$entmp;
			$cont_st=$match_rev{$chrtar}{$st}{$contig};
			$cont_en=$match_rev{$chrtar}{$en}{$contig};
			$seglen=abs($cont_en-$cont_st)+1;
			if($dis<=$len){
				$cluster{$sttmp}{$st}=0;
				$cluster_cov{$sttmp}+=$seglen;
				$entmp=$en;
			#	print "4.1\t$contig\t$chrtar\t$sttmp\t$st\t$entmp\t$dis\n";
			}else{
				$sttmp=$st;
				$entmp=$en;
				$cluster{$sttmp}{$sttmp}=0;
				$cluster_cov{$sttmp}+=$seglen;
			#	print "4.2\t$contig\t$chrtar\t$sttmp\t$st\t$entmp\t$dis\n";
			}
		}
		my $sv_ok="Disc";
		my $cluster_num=0;
		foreach my $num (sort {$a<=>$b} keys %cluster) {
			my $rate=$cluster_cov{$num}/$len;
			$cluster_num++;
			if($rate>0.6){
				$sv_ok="Main";
			}
		}
		if($cluster_num>1){
			$sv_ok.="_SV_Yes"
		}else{
			$sv_ok.="_SV_No"
		}
		
		foreach my $num (sort {$a<=>$b} keys %cluster) {
			my $rate=$cluster_cov{$num}/$len;
			my $flag="MainsLoc";
			if($rate<0.6){
				$flag="TransLoc";
			}
			$rate=int($rate*100000+0.5)/1000;
			foreach my $st (sort {$a<=>$b} keys %{$cluster{$num}}) {
				my $en=$reg{$st};
				my $cont_st=$match_rev{$chrtar}{$st}{$contig};
				my $cont_en=$match_rev{$chrtar}{$en}{$contig};
				my $seg_len=abs($cont_en-$cont_st)+1;
				if($cont_en>$cont_st){
					print $OUT1 "$sv_ok\t$flag\t$contig\t$cont_st\t$cont_en\t$seg_len\tF\t$chrtar\tcluster_st:$num\t$st\t$en\tcluster_len:$cluster_cov{$num}\t$rate\t$len\n";
				}else{
					print $OUT1 "$sv_ok\t$flag\t$contig\t$cont_st\t$cont_en\t$seg_len\tR\t$chrtar\tcluster_st:$num\t$st\t$en\tcluster_len:$cluster_cov{$num}\t$rate\t$len\n";
				}
			}
		}
		foreach my $num (sort {$a<=>$b} keys %cluster) {
			my $rate=$cluster_cov{$num}/$len;
			my %direction;
			if($rate>0.2){
				my @sts=sort {$a<=>$b} keys %{$cluster{$num}};
				foreach my $st (@sts) {
					my $en=$reg{$st};
					my $cont_st=$match_rev{$chrtar}{$st}{$contig};
					my $cont_en=$match_rev{$chrtar}{$en}{$contig};
					my $forward;
					if($cont_st<$cont_en){
						$forward="F";
					}else{
						$forward="R";
					}
					$direction{$forward}+=abs($cont_en-$cont_st)+1;
				}
				my $direct="F";
				if($direction{"F"}<$direction{"R"}){
					$direct="R";
				}
				my %inversion;
				foreach my $st (@sts) {
					my $en=$reg{$st};
					my $cont_st=$match_rev{$chrtar}{$st}{$contig};
					my $cont_en=$match_rev{$chrtar}{$en}{$contig};
					my $seg_len=abs($cont_en-$cont_st)+1;
					my $forward;
					if($cont_st<$cont_en){
						$forward="F";
					}else{
						$forward="R";
					}
				#	if($contig eq "Contig1000_pilon"){
				#		print "3\t$contig\t$cont_st\t$cont_en\t\n";
				#	}
					if($forward ne $direct){
						$inversion{$cont_st}="Inversion\t$contig\t$cont_st\t$cont_en\t$seg_len\t$forward\t$chrtar\tcluster_st:$num\t$st\t$en\tcluster_len:$cluster_cov{$num}\t$rate\t$len";
						#print $OUT2 "Inversion\t$contig\t$cont_st\t$cont_en\t$seg_len\t$forward\t$chrtar\tcluster_st:$num\t$st\t$en\tcluster_len:$cluster_cov{$num}\t$rate\t$len\n";
					}else{
						$inversion{$cont_st}="Rankingss\t$contig\t$cont_st\t$cont_en\t$seg_len\t$forward\t$chrtar\tcluster_st:$num\t$st\t$en\tcluster_len:$cluster_cov{$num}\t$rate\t$len";
						#print $OUT2 "Maindirec\t$contig\t$cont_st\t$cont_en\t$seg_len\t$forward\t$chrtar\tcluster_st:$num\t$st\t$en\tcluster_len:$cluster_cov{$num}\t$rate\t$len\n";
					}
				}
				#################################
				my %hash1;
				my %hash1_rev;
				my $index=0;
				foreach my $st (@sts) {
					my $cont_st=$match_rev{$chrtar}{$st}{$contig};
					$hash1{$cont_st}=$index;
					$hash1_rev{$index}=$cont_st;
					$index++;
				}
				my %contighash;
				my %contighash_rev;
				$index=0;
				foreach my $cst1 (sort {$a<=>$b} keys %hash1) {
					$contighash{$cst1}=$index;
					$contighash_rev{$index}=$cst1;
				#	print "1\t$contig\t$index\t$cst1\n";
					$index++;
				}
				my %tls;
				if($direct eq "F"){
					my $cnum=0;
					my %clustersub;
					for (my $i=1;$i<(@sts-1);$i++) {
						my $index1=$contighash{$match_rev{$chrtar}{$sts[$i-1]}{$contig}};
						my $index2=$contighash{$match_rev{$chrtar}{$sts[$i]}{$contig}};
						my $index3=$contighash{$match_rev{$chrtar}{$sts[$i+1]}{$contig}};
					#	print "2\t$contig\t$sts[$i-1]\t$sts[$i]\t$sts[$i+1]\t$index1\t$index2\t$index3\t$direct\n";
						if($index1<$index2 and $index2<$index3){
							$clustersub{$cnum}{$index1}=0;
							$clustersub{$cnum}{$index2}=0;
							$clustersub{$cnum}{$index3}=0;
						}else{
							$cnum++;
						}
					}
					my %cnums;
					foreach $cnum (sort keys %clustersub) {
						my @index=sort {$a<=>$b} keys %{$clustersub{$cnum}};
						my $covlen=0;
						foreach my $index (@index) {
							my $cst=$contighash_rev{$index};
							my $cen=$contig_reg{$contig}{$cst};
							$covlen+=(abs($cen-$cst)+1);
						}
						$cnums{$covlen}=$cnum;
					}

					my @cov=sort {$b<=>$a} keys %cnums;
					my $covmax=$cov[0];
					my $cnummax=$cnums{$covmax};
					#print "@cov\n$covmax\t$cnummax\n";
					foreach my $index (sort {$a<=>$b} keys %contighash_rev) {
						my $cont_stsub=$contighash_rev{$index};
						if(exists $clustersub{$cnummax}{$index}){
							#print $OUT2 "MainsLoc\t$index\t$inversion{$cont_stsub}\n";
							$tls{$cont_stsub}="MainsLoc";
						}else{
							#print $OUT2 "TransLoc\t$index\t$inversion{$cont_stsub}\n";
							$tls{$cont_stsub}="TransLoc";
						}
					}
					if(scalar(@sts)==2){
						my $cst1=$match_rev{$chrtar}{$sts[0]}{$contig};
						my $cen1=$chr_reg_pre{$contig}{$cst1};
						my $cst2=$match_rev{$chrtar}{$sts[1]}{$contig};
						my $cen2=$chr_reg_pre{$contig}{$cst2};
						my $len1=abs($cen1-$cst1)+1;
						my $len2=abs($cen2-$cst2)+1;
						if($len1>$len2){
							$tls{$cst1}="MainsLoc";
							$tls{$cst2}="TransLoc";
						}else{
							$tls{$cst2}="MainsLoc";
							$tls{$cst1}="TransLoc";
						}
					}
				}else{
					#print "3\t$contig\t$direct\n";
					my $cnum=0;
					my %clustersub;
					for (my $i=1;$i<(@sts-1);$i++) {
						my $index1=$contighash{$match_rev{$chrtar}{$sts[$i-1]}{$contig}};
						my $index2=$contighash{$match_rev{$chrtar}{$sts[$i]}{$contig}};
						my $index3=$contighash{$match_rev{$chrtar}{$sts[$i+1]}{$contig}};
						if($index1>$index2 and $index2>$index3){
							$clustersub{$cnum}{$index1}=0;
							$clustersub{$cnum}{$index2}=0;
							$clustersub{$cnum}{$index3}=0;
						}else{
							$cnum++;
						}
					}
					my %cnums;
					foreach $cnum (sort keys %clustersub) {
						my @index=sort {$a<=>$b} keys %{$clustersub{$cnum}};
						my $covlen=0;
						foreach my $index (@index) {
							my $cst=$contighash_rev{$index};
							my $cen=$contig_reg{$contig}{$cst};
							$covlen+=(abs($cen-$cst)+1);
						}
						$cnums{$covlen}=$cnum;
					}
					my @cov=sort {$b<=>$a} keys %cnums;
					my $covmax=$cov[0];
					my $cnummax=$cnums{$covmax};
					foreach my $index (sort {$a<=>$b} keys %contighash_rev) {
						my $cont_stsub=$contighash_rev{$index};
						if(exists $clustersub{$cnummax}{$index}){
							#print $OUT2 "MainsLoc\t$index\t$inversion{$cont_stsub}\n";
							$tls{$cont_stsub}="MainsLoc";
						}else{
							$tls{$cont_stsub}="TransLoc";
							#print $OUT2 "TransLoc\t$index\t$inversion{$cont_stsub}\n";
						}
					}
					if(scalar(@sts)==2){
						my $cst1=$match_rev{$chrtar}{$sts[0]}{$contig};
						my $cen1=$chr_reg_pre{$contig}{$cst1};
						my $cst2=$match_rev{$chrtar}{$sts[1]}{$contig};
						my $cen2=$chr_reg_pre{$contig}{$cst2};
						my $len1=abs($cen1-$cst1)+1;
						my $len2=abs($cen2-$cst2)+1;
						if($len1>$len2){
							$tls{$cst1}="MainsLoc";
							$tls{$cst2}="TransLoc";
						}else{
							$tls{$cst2}="MainsLoc";
							$tls{$cst1}="TransLoc";
						}
					}
				}
				foreach my $st (@sts) {
					my $cont_st=$match_rev{$chrtar}{$st}{$contig};
					print $OUT2 "$tls{$cont_st}\t$inversion{$cont_st}\n";
				}
			}
		}
		#if($contig eq "Contig1001_pilon"){last};
	}
	
	
	close $OUT1;
	close $OUT2;



