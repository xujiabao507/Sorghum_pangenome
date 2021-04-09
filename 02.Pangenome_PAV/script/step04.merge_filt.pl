use strict;
my $indir=shift;
my $genome=shift;
my $query=shift;
my $output=shift;

my $sam=(split/\//,$indir)[-1];
my $ref="/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S02.DB/$sam.fa";

my @filt=glob "$indir/*filt";
open my $OUT,">$output" or die $!;
print $OUT "$genome $query\n";
print $OUT "NUCMER\n";
foreach my $filt (@filt) {
	print "$filt\n";
	open FILT, $filt =~/gz$/ ? "zcat $filt|" : "<$filt" or die $!;
	<FILT>;<FILT>;
	while (<FILT>) {
		print $OUT "$_";
	}
	close FILT;
}
close $OUT;
