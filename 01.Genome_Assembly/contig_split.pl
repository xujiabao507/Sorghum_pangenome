use strict;
use File::Basename;
my $fa=shift;
my $interval=5000;

open FA,$fa=~/gz$/?"zcat $fa|":"<$fa" or die $!;
$/=">";
<FA>;
$/="\n";
while(my $chr=<FA>){
	chomp $chr;
	$chr=(split/\s+/,$chr)[0];
	$/=">";
	my $seq=<FA>;
	$seq=uc($seq);
	$seq=~s/[>\s+]//g;
	$/="\n";
	my $length=length($seq);
	for (my $loc=1;$loc<$length ;$loc+=$interval) {
		my $loc2=$loc+2*$interval-1;
		if($loc2>$length){$loc2=$length};
		my $len=$loc2-$loc+1;
		my $subseq=substr($seq,$loc-1,$len);
		my $qual="J" x $len;
		#print ">$chr\_$loc\_$loc2\n$subseq\n";
		print "\@$chr\_$loc\_$loc2\n$subseq\n+\n$qual\n";
	}
}
close FA;
 





