use strict;
my $bwamem=shift;

my %maps;
my %contigmaps;
my %maplens;
my %nums=0;
my $line=0;
open SAM,$bwamem=~/gz$/?"zcat $bwamem|":"<$bwamem" or die $!;
while (<SAM>) {
	chomp;
	my @a=split/\s+/,$_;
	my $id=$a[0];
	my @ids=split/\_/,$id;
	my $contig="$ids[0]";
	my $chr=$a[4];
	my $cst=$ids[-2];
	my $cen=$ids[-1];
	$contigmaps{$contig}{$cen}=0;
	my $rate=$a[1];
	my $chrst=$a[5];
	my $len=(split/\//,$a[2])[0];
	next if($chr!~/Chr/);
	$maps{$contig}{$chr}{$chrst}="$cst\t$cen\t$a[2]";
	$maplens{$contig}{$chr}+=$len;
	$nums{$contig}{$chr}++;
	$line++;
}
close SAM;
foreach my $contig (sort keys %maplens) {
	my $flag=0;
	my @loci=sort {$b<=>$a} keys %{$contigmaps{$contig}};
	my $len=$loci[0];
	my $num=scalar(@loci);
	foreach my $chr (sort {$maplens{$contig}{$b}<=>$maplens{$contig}{$a}} keys %{$maplens{$contig}}) {
		my $num_chr=0;
		my $mean=int($maplens{$contig}{$chr}/$nums{$contig}{$chr});
		foreach my $chrst (sort {$a<=>$b} keys %{$maps{$contig}{$chr}}) {
			$num_chr++;
		}
		my $rate=$num_chr/$num;
		#next if($rate<0.5 or $mean<6000);
		foreach my $chrst (sort {$a<=>$b} keys %{$maps{$contig}{$chr}}) {
			print "$contig\t$len\t$chr\t$chrst\t$maps{$contig}{$chr}{$chrst}\t$maplens{$contig}{$chr}\t$nums{$contig}{$chr}\t$mean\t$num\t$num_chr\t$rate\n";
		}
		last;
	}
}

