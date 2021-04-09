use strict;
my $depthfile=shift;
my $output=shift;

print "$depthfile is start\n";
my %genedepth;
my %genecover;
my %chr2gene;
open my $DPH,$depthfile=~/gz$/?"zcat $depthfile|":"<$depthfile" or die $!;
while (<$DPH>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $depth=$a[2];
	my $gene=$a[3];
	$genedepth{$gene}{$depth}++;
	$genecover{$gene}++;
	$chr2gene{$chr}{$gene}=0;
}
close $DPH;
print "$depthfile is over\n";


open my $OUT,">$output" or die $!;
foreach my $chr (sort keys %chr2gene) {
	foreach my $gene (sort keys %{$chr2gene{$chr}}) {
		my $sum=0;
		foreach my $dep (sort {$a<=>$b} keys %{$genedepth{$gene}}) {
			$sum+=($dep*$genedepth{$gene}{$dep});
		}
		my $mean=$sum/$genecover{$gene};
		my $var=0;
		foreach my $dep (sort {$a<=>$b} keys %{$genedepth{$gene}}) {
			$var+=(($dep-$mean)**2)*$genedepth{$gene}{$dep};
		}
		$var=$var/$genecover{$gene};
		my $sdev=$var**0.5;
		print $OUT "$chr\t$gene\t$mean\t$sdev\t$genecover{$gene}\t$var\t$sum\n";
	}
}
close $OUT;






 