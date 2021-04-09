use strict;

my $map=shift;
my $read=shift;
my $out=shift;

print "$map start\n";
my %map;
open MAP,"zcat $map|" or die $!;
while (<MAP>) {
	chomp;
	my @a=split/\s+/;
	my $read=$a[0];
	my ($chr,$st)=(split/\_/,$read)[0,1,2];
	my $tarchr=$a[2];
	my $tarloc=$a[3];
	$map{$chr}{$st}="$tarchr\t$tarloc";
#	print "$chr\t$st\t$tarchr\t$tarloc\n";
}
close MAP;
open OUT,">$out" or die  $!;
print "$read start\n";
open READ,"$read" or die $!;
while (<READ>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $st=$a[1];
	if(exists $map{$chr}{$st}){
		print OUT "$_\t$map{$chr}{$st}\n";
	}else{
		print OUT "$_\tNA\tNA\n";
	}
}
close READ;
close OUT;







