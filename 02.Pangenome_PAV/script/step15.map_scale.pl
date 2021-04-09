use strict;
my $mummer_map=shift;
my $blast_map=shift;
my $output=shift;

open my $OUT,">$output" or die $!;

my %nucmer;
my %lensca;
open my $MUMMER,$mummer_map=~/gz$/?"zcat $mummer_map|":"<$mummer_map" or die $!;
while (<$MUMMER>) {
	chomp;
	my @a=split/\s+/;
	my $ref_loc1=$a[0];
	my $ref_loc2=$a[1];
	my $len=$a[6];
	my $scalen=$a[12];
	if($len<5000){next};
	my $ref=$a[19];
	my $sca=$a[20];
	$lensca{$sca}=$scalen;
	$nucmer{$sca}{$ref_loc1}=$ref;
	$nucmer{$sca}{$ref_loc2}=$ref;
	#print "0\t$sca\t$scalen\t$ref\t$ref_loc1\t$ref_loc2\n";
}
close $MUMMER;
my %nucmer_scale;
foreach my $sca (sort keys %nucmer) {
	my @loc=sort {$a<=>$b} keys %{$nucmer{$sca}};
	my $loc=$loc[0];
	my $st=$loc[0]-$lensca{$sca};
	my $en=$loc[-1]+$lensca{$sca};
	my $ref=$nucmer{$sca}{$loc};
	$nucmer_scale{$sca}="$ref\t$st\t$en";
	#print "1\t$ref\t$sca\t$st\t$en\t$lensca{$sca}\n";
}

my %blast;
open my $BLAST,$blast_map=~/gz$/?"zcat $blast_map|":"<$blast_map" or die $!;
while (<$BLAST>) {
	chomp;
	my @a=split/\s+/;
	my $query=$a[0];
	my @query;
	if ($query=~/part/) {
		@query=split/_part/,$query;
	}else{
		@query=split/_total/,$query;
	}
	my $contig=$query[0];
	if(!exists $nucmer_scale{$contig}){next};
	my ($chr,$st,$en)=(split/\s+/,$nucmer_scale{$contig})[0,1,2];
	my $chr_tar=$a[1];
	if($chr ne $chr_tar){next};
	my $ind=$a[2];
	if($ind<90){next};
	my $len_query=$a[3];
	if($len_query<100){next};
	my $loc_query=$a[6];
	my $loc1=$a[8];
	my $loc2=$a[9];
	my $min=$loc1;
	my $max=$loc2;
	if($min>$loc2){
		$min=$loc2;
		$max=$loc1;
	}
	if($st<$min and $en>$max){
		$blast{$query}{$loc_query}{$_}=0;
	}
}
close $BLAST;
foreach my $query (sort keys %blast) {
	foreach my $loc (sort {$a<=>$b} keys %{$blast{$query}}) {
		foreach my $str (sort keys %{$blast{$query}{$loc}}) {
			print $OUT "$str\n";
		}
	}
}
close $OUT;
















