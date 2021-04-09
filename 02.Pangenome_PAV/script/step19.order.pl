use strict;
##########  /hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S10.PAV/S273/step05.S2732ref/03.inter
my $map=shift;
my $orderfile=shift;
my $svout=shift;
open ORDER ,">$orderfile" or die  $!;
open SV,">$svout" or die  $!;
print "$map start\n";
#open CLU,">cluster.txt" or die $!;
my %sammap;
my %refmap;
my %sam2ref;
my %ref2sam;
my %lens;
open MAP,$map or die $!;
while (<MAP>) {
	chomp;
	next if($_!~/\d/);
	my @a=split/\s+/;
	my $chr=$a[0];
#	next if($chr ne "Chr01");
	my $len=$a[2];
	my $qst=$a[4];
	my $qen=$a[5];
	my $tst=$a[6];
	my $ten=$a[7];
	$sammap{$chr}{$qst}=$qen;
	$refmap{$chr}{$tst}=$ten;
	$sam2ref{$chr}{$qst}=$tst;
	$ref2sam{$chr}{$tst}=$qst;
	$lens{$chr}{$qst}=$len;
}
close MAP;
print "loading is ok\n";
########################## cross number  :: delete high number and remain cross number == 0 
my %qst2tst;
foreach my $chr (sort keys %sam2ref) {
	my %crossnum;
	my %crossloc;
	my %order;
	my @qst=sort {$a<=>$b} keys %{$sam2ref{$chr}};
	%qst2tst=();
	foreach my $qst (@qst) {
		my $tst=$sam2ref{$chr}{$qst};
		$qst2tst{$qst}=$tst;
	}
	print "start first\t".localtime()."\n";
	for (my $i=0;$i<@qst ;$i++) {
		my $qst1=$qst[$i];
		my $tst1=$qst2tst{$qst1};
		for (my $j=$i+1;$j<@qst ;$j++) {
			my $qst2=$qst[$j];
			my $tst2=$qst2tst{$qst2};
			if(($qst2<$qst1 and $tst2>$tst1) or ($qst2>$qst1 and $tst2<$tst1)){
				$crossloc{$qst1}{$qst2}=0;
				$crossloc{$qst2}{$qst1}=0;
				$crossnum{$qst1}++;
				$crossnum{$qst2}++;
			}
		}
	}
	print "end first\t".localtime()."\n";
	foreach my $qst (@qst) {
		if(!exists $crossloc{$qst}){
			$order{$qst}=0;
			print "primary order0\t$qst\t$qst2tst{$qst}\n";
		}
	}
	my %delete;
	my %relate;
	my $num=scalar(keys %crossloc);
	my $ind=0;
	while (1) {
		$ind++;
		print "Start\n";
		%relate=();
		my @qst=sort {$a<=>$b} keys %crossloc;
		my $qstnum1=scalar(keys %crossloc);
		print "$ind before delete $qstnum1\n";
		my @qst_crn=sort {$crossnum{$b}<=>$crossnum{$a}} keys %crossnum;
		my $q1=$qst_crn[0];
		foreach my $qst (@qst_crn) {
			if($crossnum{$qst}>10000){
				my $ind2=0;
				foreach my $qst2 (sort keys %{$crossloc{$qst}}) {
					$ind2++;
					delete $crossloc{$qst2}{$qst};
					$delete{$qst}=0;
					$relate{$qst2}=0;
				}
				delete $crossloc{$qst};
				delete $crossnum{$qst};
				next;
			}else{
				my $ind2=0;
				foreach my $qst2 (sort keys %{$crossloc{$qst}}) {
					$ind2++;
					delete $crossloc{$qst2}{$qst};
					my $qst2num=0;
					if(scalar( keys %{$crossloc{$qst2}})>0){
						$qst2num=1;
					}
					if($qst2num==0){
						delete $crossloc{$qst2};
						delete $crossnum{$qst2};
					}
					$delete{$qst}=0;
					$relate{$qst2}=0;
				}
		#		print "$ind\tdelete2\tcross:$crossnum{$qst}\t$qst\t$qst2tst{$qst}\n";
				delete $crossloc{$qst};
				delete $crossnum{$qst};
				last;
			}
		}
		my $qstnum2=scalar(keys %crossloc);
		print "$ind after delete $qstnum2 \n";
		my @val_tmp;
		foreach my $qst (sort {$a<=>$b} keys %relate) {
			if(!exists $crossloc{$qst} and !exists $delete{$qst}){
				$order{$qst}=0;
				print "New_loc\t$qst\t$qst2tst{$qst}\n";
			}
		}
		my $num=scalar(keys %crossloc);
		my $ordnum=scalar(keys %order);
		#my @order=sort {$a<=>$b} keys %order;
		#my $ordnum=scalar(@order);
		#my $ordind=0;
		#foreach my $ord (@order) {
		#	$ordind++;
		#	print "Order\tindex:$ind\t$ordind\tcrossnum:$num\tordnum:$ordnum\t$ord\t$qst2tst{$ord}\n";
		#}
		print "summary_index:\t$ind\t$num\t$ordnum\n";
		print "End\n\n";
		if($num==0){
			last;
		}
	}
	my @sort=sort {$a<=>$b} keys %order;
	foreach my $st (@sort) {
		my $tst=$qst2tst{$st};
		print ORDER "$chr\t$st\t$sammap{$chr}{$st}\t$tst\t$refmap{$chr}{$tst}\n";
	}
	my $indexnum=scalar(@sort);
	my @samloc=sort {$a<=>$b} keys %{$sam2ref{$chr}};
	foreach my $loc (@samloc) {
		if(!exists $order{$loc}){
			for (my $h=0;$h<$indexnum-1 ;$h++) {
				if($loc>$sort[$h] and $loc<$sort[$h+1]){
					my $qst=$loc;
					my $qen=$sammap{$chr}{$qst};
					my $tst=$qst2tst{$loc};
					my $ten=$refmap{$chr}{$tst};
					my $len1=$qen-$qst;
					my $len2=$ten-$tst;
					my $tag="Normal";
					if($ten<$tst){
						$tag="Inversion";
					}
					print SV "Translocation:\t$chr\t$sort[$h]\t$sort[$h+1]\t$loc\t$sammap{$chr}{$loc}\t$len1\t$tst\t$ten\t$len2\t$tag\n";
					last;
				}
			}
		}else{
			my $qst=$loc;
			my $qen=$sammap{$chr}{$qst};
			my $tst=$qst2tst{$qst};
			my $ten=$refmap{$chr}{$tst};
			my $len1=$qen-$qst;
			my $len2=$ten-$tst;
			my $tag="Normal";
			if($ten<$tst){
				$tag="Inversion";
			}
			print SV "Order:\t$chr\tNA\tNA\t$qst\t$qen\t$len1\t$tst\t$ten\t$len1\t$tag\n";
		}
	}
}





######################### 
#close CLU;
close ORDER;
close SV;


sub crossnum{
	my ($tagst,$arr)=@_;
	my $tagen=$qst2tst{$tagst};
	my @arr=@$arr;
	my $cronum=0;
	for (my $l=0;$l<@arr ;$l++) {
		my $stloc=$arr[$l];
		my $enloc=$qst2tst{$stloc};
		if($stloc eq $tagst){next};
		if(($stloc<$tagst and $enloc>$tagen) or ($stloc>$tagst and $enloc<$tagen)){
			$cronum++;
		}
	}
	return $cronum 
}




