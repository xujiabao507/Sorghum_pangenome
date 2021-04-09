use strict;
open OUT,">pav.xls" or die $!;

my $sam_specific_dir="/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S10.PAV/";
my @samspe=glob "$sam_specific_dir/*/step11.*.unmap_reads_merge.xls";

my %pavstat;
foreach my $spe (@samspe) {
	print "$spe\n";
	my $sam=(split/\//,$spe)[-2];
	open SPE,$spe or die $!;
	while (<SPE>) {
		chomp;
		next if($_!~/\w/);
		my @a=split/\s+/;
		my $len=$a[3];
		if($len>100){
			$pavstat{$sam}+=$len;
		}
	}
	close SPE;
}
foreach my $sam (sort keys %pavstat) {
	print OUT "Sam-specifit\t$sam\t$pavstat{$sam}\n";
}

my @pavs=glob "$sam_specific_dir/*/step10.*.pav_cov.xls";

my %pavstat2;
foreach my $pav (@pavs) {
	my $sam=(split/\//,$pav)[-2];
	print "$pav\n";
	open PAV,$pav or die $!;
	while (<PAV>) {
		chomp;
		next if($_!~/\w/);
		my @a=split/\s+/;
		my $dep=$a[4];
		if($dep==0){
			$pavstat2{$sam}++
		}
	}
	close PAV;
}

foreach my $sam (sort keys %pavstat2) {
	print OUT "ref-specifit\t$sam\t$pavstat2{$sam}\n";
}

close OUT;
