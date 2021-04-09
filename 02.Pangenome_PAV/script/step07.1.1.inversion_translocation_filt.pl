use strict;
my $obo=shift;
 
my %chrs;
my %cont;
	open my $OUT,">$obo.filt" or die $!;
	open OBO, $obo =~/gz$/ ? "zcat $obo|" : "<$obo" or die $!;
	while (<OBO>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		my @a=split/\s+/;
		my $type1=$a[0];
		my $type2=$a[1];
		my $chr=$a[2];
		$cont{$chr}.="$_\n";
		if($type1 eq "TransLoc" or $type2 eq "Inversion" or $type2 eq "TransLoc"){
			$chrs{$chr}=0;
		}
	}
	foreach my $chr (sort keys %cont) {
		if(exists $chrs{$chr}){
			print $OUT "$cont{$chr}";
		}
	}
	close OBO;
	close $OUT;



