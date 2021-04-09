use strict;
my $chrlist=shift;
my $map=shift;
my $mapout=shift;

my %chrlen;
open CL,$chrlist or die $!;
while (<CL>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $len=$a[1];
	$chrlen{$chr}=$len;
}
close CL;


 
my %pavs;
open MAP,$map or die $!;
while (<MAP>) {
	chomp;
	my @a=split/\s+/,$_;
	#my @b=split/\_/,$a[0];
	my $chr1=$a[0];
	my $st1=$a[1];
	my $en1=$a[2];
	$pavs{$chr1}{$st1}=$en1;
}
close MAP; 


open OUT,">$mapout" or die $!;
foreach my $chr (sort keys %pavs) {
	my %unmap;
	my $st=0;
	my @map=sort {$a<=>$b} keys %{$pavs{$chr}};
	#print "@map\n";
	if($st<$map[0]){
		$unmap{$st}=($map[0]-1);
	#	print "0\t$chr\t$st\t$map[0]-1\n";
	}
	for (my $i=0;$i<(scalar(@map)-1) ;$i++) {
		my $st1=$map[$i];
		my $en1=$pavs{$chr}{$st1};
		my $st2=$map[$i+1];
		$unmap{$en1+1}=($st2-1);
	#	print  "1\t$chr\t$en1+1\t$st2-1\n";
	}
	my $en=$pavs{$chr}{$map[-1]};
	if($en<$chrlen{$chr}){
		$unmap{$en+1}=$chrlen{$chr};
	#	print  "2\t$chr\t$en+1\t$chrlen{$chr}\n";
	}
	foreach my $st (sort {$a<=>$b} keys %unmap) {
		my $len=$unmap{$st}-$st+1;
		print OUT "$chr\t$st\t$unmap{$st}\t$len\n";
	}
}
close OUT;



 