use strict;
my $chrin=shift;
my $pavdir=shift;

if(!-e "$pavdir/PAV_Merge"){
	`mkdir -p "$pavdir/PAV_Merge"`;
}
open PAVOUT,">$pavdir/PAV_Merge/Pav.$chrin.merge.list" or die $!;
my @pav=glob "$pavdir/*/*.ref.unmap.merge";
my %pavs;
my %pav;
my %scales;
my @sam;
foreach my $pav (@pav) {
	my $base=(split/\//,$pav)[-1];
	my $sam=(split/\./,$base)[0];
	print "$pav\n";
	push @sam,$sam;
	open PAV,$pav or die $!;
	while (<PAV>) {
		chomp;
		my @a=split/\s+/;
		my $chr=$a[1];
		next if($chr ne $chrin);
		my $st=$a[5];
		my $en=$a[6];
#		print "$chr\t$st\t$";
		$pavs{$st}=$en;
		$scales{$st}=0;
		$scales{$en}=0;
		for (my $loc=$st;$loc<=$en ;$loc++) {
			$pav{$loc}{$sam}=0;
		}
	}
	close PAV;
}
my @loc=sort keys %pav;
my $num=scalar(@loc);
if($num==0){
	exit;
}

my @scale=sort {$a<=>$b} keys %scales;
my $samstr=join "\t",@sam;
print PAVOUT"CHR\tP1\tP2\tLen\t$samstr\n";
print "@scale\n";
for (my $i=0;$i<@scale-1 ;$i++) {
	my $loc1=$scale[$i];
	my $loc2=$scale[$i+1];
	my $loc3=$loc1+1;
	my @val;
	my $flag=0;
	foreach my $sam (@sam) {
		if(exists $pav{$loc3}{$sam}){
			push @val,"Y";
			$flag=1;
		}else{
			push @val,"N";
		}
	}
	my $vals=join "\t",@val;
	if($flag==1){
		my $len=$loc2-$loc1+1;
		print PAVOUT "$chrin\t$loc1\t$loc2\t$len\t$vals\n";
	}
}
close PAVOUT;























