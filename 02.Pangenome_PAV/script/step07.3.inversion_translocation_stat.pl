use strict;

my $sv_tl=shift;
my $stat=shift;

open ITL, $sv_tl =~/gz$/ ? "zcat $sv_tl|" : "<$sv_tl" or die $!;
my (%outsv,%insv,%outnum,%innum,%inv,%inv_num);
while (<ITL>) {
	chomp;
	my @a=split/\s+/;
	my $tag=$a[0];
	next if($tag!~/Summary_End/);
	my $scf=$a[1];
	my $pos=$a[2];
	my $type=$a[6];
	my $type1=$a[7];
	my $type2=$a[8];
	my $read1=$a[11];
	my $read2=$a[13];
	my $readall=$read1+$read2;
	
	if($type eq "OutTL" and $type2 eq "TransLoc"){
		$outsv{0}++;
		$outnum{0}{$scf}=0;
		if($readall>0){
			$outsv{1}++;
			$outnum{1}{$scf}=0;
		}
	}
	if($type eq "IN_SV" and $type1 eq "TransLoc"){
		$insv{0}++;
		$innum{0}{$scf}=0;
		if($readall>0){
			$insv{1}++;
			$innum{1}{$scf}=0;
		}
	}
	if($type eq "IN_SV" and $type1 ne "TransLoc" and $type2  eq "Inversion"){
		$inv{0}++;
		$inv_num{0}{$scf}=0;
		if($readall>0){
			$inv{1}++;
			$inv_num{1}{$scf}=0;
		}
	}
}
close ITL;

open STAT,">$stat" or die $!;
my $scf_out_all=scalar(keys %{$outnum{0}});
my $scf_out_high=scalar(keys %{$outnum{1}});
print STAT "Type\tScf_num_All\tScf_num_readsup\tVar_num\tVar_num_Readsup\n";
print STAT "OutSV_Translocation\t$scf_out_all\t$scf_out_high\t$outsv{0}\t$outsv{1}\n";
my $scf_in_all=scalar(keys %{$innum{0}});
my $scf_in_high=scalar(keys %{$innum{1}});
print STAT "InSV_Translocation\t$scf_in_all\t$scf_in_high\t$insv{0}\t$insv{1}\n";
my $scf_in_inv_all=scalar(keys %{$inv_num{0}});
my $scf_in_inv_high=scalar(keys %{$inv_num{1}});
print STAT "InSV_Inversion\t$scf_in_inv_all\t$scf_in_inv_high\t$inv{0}\t$inv{1}\n";
close STAT;

















