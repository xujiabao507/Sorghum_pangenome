use strict;
my $core_depth=shift;
my $depthnum=shift;
my $output=shift;
open OUT,">$output" or die $!;
my %scale;
open CD,$core_depth=~/gz$/ ? "zcat $core_depth|" : "<$core_depth" or die $!;
my $num=0;
my $depth_tmp=0;
my $pos_tmp=0;
my $pos_start=0;
while (<CD>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $pos=$a[1];
	my $depth=$a[2];
	if($depth eq $depthnum){
		$num++;
		if($num eq 1){
			$depth_tmp=$depth;
			$pos_tmp=$pos;
			$pos_start=$pos;
			$scale{$chr}{$pos_start}=$pos;
			next;
		}
		if(($pos-$pos_tmp)==1){
			$scale{$chr}{$pos_start}=$pos;
			$pos_tmp=$pos;
		}else{
			$pos_start=$pos;
			$scale{$chr}{$pos_start}=$pos;
			$pos_tmp=$pos;
		}
	}
}
close CD;
foreach my $chr (sort keys %scale) {
	foreach my $st (sort {$a<=>$b} keys %{$scale{$chr}}) {
		my $en=$scale{$chr}{$st};
		my $len=$en-$st+1;
		print OUT "$chr\t$st\t$en\t$len\n";
	}
}
close OUT;









