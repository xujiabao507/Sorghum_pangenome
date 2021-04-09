use strict;
my $readmap =shift;
my $sv_out=shift;
my $svdetail_out=shift;

open SVOUT,">$sv_out" or die $!;
open SVDOUT,">$svdetail_out" or die $!;


my %mapread;
my %tags;
 
open MR, $readmap =~/gz$/ ? "zcat $readmap|" : "<$readmap" or die $!;
while(<MR>){
	chomp;
	my @a=split/\s+/;
	my $tag=$a[0];
	$tag=~s/\://g;
	my $chr=$a[1];
	my $st=$a[2];
	my $en=$a[3];
	my $tchr=$a[4];
	my $tloc=$a[5];
	$mapread{$chr}{$st}="$tchr\t$tloc";
	$tags{$chr}{$st}=$tag;
}
close MR;

foreach my $chr (sort keys %mapread) {
	my %svs;
	my @st=sort {$a<=>$b} keys %{$mapread{$chr}};
	my $svnum=0;
	for (my $i=0;$i<(scalar(@st)-1);$i++) {
		my $st1=$st[$i];
		my $st2=$st[$i+1];
		if($tags{$chr}{$st1}=~/Order/ and $tags{$chr}{$st2}=~/SV/){
			$svnum++;
		}
		if($tags{$chr}{$st1}=~/SV/){
			$svs{$svnum}{$st1}=0;
		}
	}
	
	
	foreach my $svnum (sort {$a<=>$b} keys %svs) {
		my %main;
		my %maps;
		my @svscale=sort {$a<=>$b} keys %{$svs{$svnum}};
		foreach my $st (@svscale) {
			my ($tchr,$tloc)=(split/\s+/,$mapread{$chr}{$st})[0,1];
			$main{$tchr}++;
			$maps{$tchr}{$st}=$tloc;
		#	print "0\t$chr\t$svscale[0]\t$svscale[-1]\t$svnum\t$st\t$mapread{$chr}{$st}\n";
		}
		my $mainchr;
		my %rate;
		foreach my $tchr (sort {$main{$b}<=>$main{$a}} keys %main) {
			my $rate=$main{$tchr}/scalar(@svscale);
			if($rate>0.5){
				my %mainreg;
				my %mainread;
				my @mainloc=sort {$a<=>$b} keys %{$maps{$tchr}};
				my $len=$mainloc[-1]-$mainloc[0]+500;
				#print SVDOUT "$svscale[0]\t$svscale[-1]\t@mainloc\n";
				foreach my $st (@mainloc) {
					my $tst=$maps{$tchr}{$st};
					my $int=int($tst/10000)*10000;
					$mainreg{$int}++;
					$mainread{$int}{$st}=0;
					print SVDOUT "1\t$chr\t$svnum\t$st\t$tchr\t$rate\t$maps{$tchr}{$st}\t$int\n";
				}
				my @k100=sort {$mainreg{$b}<=>$mainreg{$a}} keys %mainreg;
				my $k100=$k100[0];
				my $sum=0;
				my $num=0;
				foreach my $st (sort keys %{$mainread{$k100}}) {
					$sum+=$maps{$tchr}{$st};
					$num++;
				}
				my $mean=$sum/$num;
				my %result;
				print SVDOUT "x\t$svscale[0]\t$svscale[-1]\t$mean\t@mainloc\n";
				foreach my $st (@mainloc) {
					my $tst=$maps{$tchr}{$st};
					if(abs($tst-$mean)>2*$len){next};
					$result{$st}=$tst;
					print SVDOUT "2\t$chr\t$svnum\t$st\t$mean\t$tchr\t$rate\t$maps{$tchr}{$st}\n";
				}
				my @resloc=sort {$a<=>$b} keys %result;
				if(scalar(@resloc)<1){next};
				my $len2=$resloc[-1]-$resloc[0]+500;
				my $st1=$resloc[0];
				my $st2=$resloc[-1];
				my $tst1=$maps{$tchr}{$st1};
				my $tst2=$maps{$tchr}{$st2};
				my $len3=abs($tst2-$tst1)+500;
				if($tst1>$tst2){
					print SVOUT "TL\tINVERS\t$chr\t$svscale[0]\t$svscale[-1]\t$st1\t$st2\t$tchr\t$tst1\t$tst2\t$len2\t$len3\t$rate\n";
					print SVDOUT "3\tTL\tINVERS\t$chr\t$svscale[0]\t$svscale[-1]\t$st1\t$st2\t$tchr\t$tst1\t$tst2\t$len2\t$len3\t$rate\n";
				}else{
					print SVOUT "TL\tNORMAL\t$chr\t$svscale[0]\t$svscale[-1]\t$st1\t$st2\t$tchr\t$tst1\t$tst2\t$len2\t$len3\t$rate\n";
					print SVDOUT "3\tTL\tINVERS\t$chr\t$svscale[0]\t$svscale[-1]\t$st1\t$st2\t$tchr\t$tst1\t$tst2\t$len2\t$len3\t$rate\n";
				}
			}
			last;
		}
	}
	
}
close SVOUT;
close SVDOUT;







