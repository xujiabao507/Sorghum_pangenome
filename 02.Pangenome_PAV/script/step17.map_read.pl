use strict;

my $map=shift;
my $out=shift;
open OUT,">$out" or die  $!;
print "$map start\n";
my %readmap;
open MAP,"zcat $map|" or die $!;
while (<MAP>) {
	chomp;
	my @a=split/\s+/;
	my $read=$a[0];
	my @b=split/\_/,$read;
	my $chr_read=$b[0];
	my $st_read=$b[1];
	my $en_read=$b[2];
	my $cigar=$a[1];
	my $chr_tar=$a[2];
	my $st_tar=$a[3];
	$readmap{$chr_read}{$st_read}="$en_read\t$cigar\t$chr_tar\t$st_tar";
#	print OUT "$chr_read\t$st_read\t$en_read\t$cigar\t$chr_tar\t$st_tar\n";
}
close MAP;
foreach my $chr (sort keys %readmap) {
	foreach my $st (sort {$a<=>$b} keys %{$readmap{$chr}}) {
		print OUT "$chr\t$st\t$readmap{$chr}{$st}\n";
	}
}
close OUT;

