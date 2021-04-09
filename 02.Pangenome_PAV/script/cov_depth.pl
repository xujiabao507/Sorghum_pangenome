use strict;


my $depthfile=shift;
my $out=shift;
my %chrs;
my $chrold="Chr01";
my %depth;
open OUT,">$out" or die $!;
open DF,$depthfile or die $!;
while (<DF>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $pos=$a[1];
	my $dep=$a[2];
	if(!exists $chrs{$chr}){
		#print "1\t$chr\t$pos\n";
		my @loc=sort {$a<=>$b} keys %depth;
		my %scale;
		my $st=0;
		my $en=0;
		foreach my $loc (@loc) {
			my $loc1=$loc-1;
			my $loc2=$loc;
			my $loc3=$loc+1;
			if(!exists $depth{$loc1} and exists $depth{$loc2} and exists $depth{$loc3}){
				$st=$loc2;
				$scale{$st}=0;
			}
			if(exists $depth{$loc1} and exists $depth{$loc2} and !exists $depth{$loc3}){
				$en=$loc2;
				$scale{$st}=$en;
			}
		}
		my $sum=0;
		foreach my $st (sort {$a<=>$b} keys %scale) {
			my $len=($scale{$st}-$st)+1;
			$sum+=$len;
			print OUT "$chrold\t$st\t$scale{$st}\t$len\t$sum\n";
		}
		$chrs{$chr}=0;
		%depth=();
		#next;
	}else{
		$chrold=$chr;
		# print "0\t$chr\t$pos\n";
	}
	if($dep==0){
		$depth{$pos}=0;
	}
	#next if($pos>100000);
}
close DF;


		my @loc=sort {$a<=>$b} keys %depth;
		my %scale2;
		my $st=0;
		my $en=0;
	#	print "$loc[0]\t$loc[-1]\n";
		foreach my $loc (@loc) {
			my $loc1=$loc-1;
			my $loc2=$loc;
			my $loc3=$loc+1;
			if(!exists $depth{$loc1} and exists $depth{$loc2} and exists $depth{$loc3}){
				$st=$loc2;
				$scale2{$st}=0;
	#			print "2\t$st\t";
			}
			if(exists $depth{$loc1} and exists $depth{$loc2} and !exists $depth{$loc3}){
				$en=$loc2;
				$scale2{$st}=$en;
	#			print "$en\n";
			}
		}
		my $sum2=0;
		my $interst=1;
		foreach my $st (sort {$a<=>$b} keys %scale2) {
			my $len=($scale2{$st}-$st)+1;
			my $len2=$st-$interst+1;
			$sum2+=$len;
			print OUT "$chrold\t$st\t$scale2{$st}\t$len\t$sum2\n";
			$interst=$scale2{$st}+1;
		}
close OUT;

