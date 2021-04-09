use strict;
my $varfile=shift;
my $output=shift;

my %ins;
my %del;
my %locs;
open IN,"$varfile" or die $!;
while (<IN>) {
	chomp;
	my @a=split/\s+/;
	my $loc1=$a[0];
	my $ref=$a[1];
	my $var=$a[2];
	my $loc2=$a[3];
	my $tar=$a[8];
	my $contig=$a[9];
	if($ref eq "."){
		$ins{$tar}{$loc1}.=$var;
	}
	if($var eq "."){
		$del{$contig}{$loc2}.=$ref;
		$locs{$contig}{$loc2}{$tar}{$loc1}=0;
	}
}
close IN;

my %del2;
my %del2str;
foreach my $contig (sort keys %locs) {
	foreach my $loc (sort keys %{$locs{$contig}}) {
		foreach my $tar (sort keys %{$locs{$contig}{$loc}}) {
			my @locs=sort {$a<=>$b} keys %{$locs{$contig}{$loc}{$tar}};
			$del2{$tar}{$locs[0]}=$locs[-1];
			#print "$tar\t$locs[0]\t$locs[-1]\n";
			$del2str{$tar}{$locs[0]}=$del{$contig}{$loc};
		}
	}
}

my %lendisdel;
my $sumdel=0;
foreach my $chr (sort keys %del2) {
	foreach my $loc (sort {$a<=>$b} keys %{$del2{$chr}}) {
		my $len=$del2{$chr}{$loc}-$loc+1;
	#	print "$len\n";
		$lendisdel{$len}++;
		$sumdel++;
	}
}
my %lendisins;
my $sumins=0;
foreach my $chr (sort keys %ins) {
	foreach my $loc (sort {$a<=>$b} keys %{$ins{$chr}}) {
		my $len=length($ins{$chr}{$loc});
		$lendisins{$len}++;
		$sumins++;
	}
}

foreach my $num (sort {$a<=>$b} keys %lendisdel) {
	my $rate=$lendisdel{$num}/$sumdel;
	print "Deletion\t$num\t$lendisdel{$num}\t$sumdel\t$rate\n";
}

foreach my $num (sort {$a<=>$b} keys %lendisins) {
	my $rate=$lendisins{$num}/$sumins;
	print "Insertion\t$num\t$lendisins{$num}\t$sumins\t$rate\n";
}