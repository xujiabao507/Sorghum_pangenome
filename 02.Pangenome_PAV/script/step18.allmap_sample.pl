use strict;
my $contigfile=shift;
my $obo=shift;
my $blastmap=shift;
my $blastmap_stat=shift;
my $output=shift;

open my $OUT,">$output" or die $!;
my %querylen;
open CONTIG, $contigfile =~/gz$/ ? "zcat $contigfile|" : "<$contigfile" or die $!;
$/=">";
<CONTIG>;
$/="\n";
while(<CONTIG>){
	chomp;
	my $name=(split/\s+/,$_)[0];
	$/=">";
	my $seq=<CONTIG>;
	$/="\n";
	$seq=~s/\s+//g;
	$seq=~s/\>//g;
	my $len=length($seq);
	$querylen{$name}=$len;
}
close CONTIG;
	
	my %smallcontig;
	my %match;
	my %sam_match;
	my %nucmer;
	open OBO, $obo =~/gz$/ ? "zcat $obo|" : "<$obo" or die $!;
	while (<OBO>) {
		chomp;
		$_=~s/(^\s+|\s+$)//g;
		my @a=split/\s+/;
		my $chr=$a[19];
		my $contig=$a[20];
		my $maprate=$a[21];
		my $contig_st=$a[3];
		my $contig_en=$a[4];
		my $chr_st=$a[0];
		my $chr_en=$a[1];
#		if($contig_st>$contig_en){
#			my $tmp=$contig_st;
#			$contig_st=$contig_en;
#			$contig_en=$tmp;
#		}
		my $len_contig=$a[12];
		my $len=abs($chr_en-$chr_st)+1;
		next if($len<5000);
		if($maprate<20){
			print "$_\n";
			next;
		}
		$sam_match{$contig}{$contig_st}=$len;
		$match{$contig}{$chr}{$chr_st}=0;
		$match{$contig}{$chr}{$chr_en}=0;
		$nucmer{$contig}{$contig_st}{$contig_en}="$chr\t$chr_st\t$chr_en\t$len\tnucmer";
	}
	close OBO;
	

my %blast;
open BLAST, $blastmap_stat =~/gz$/ ? "zcat $blastmap_stat|" : "<$blastmap_stat" or die $!;
while (<BLAST>) {
	chomp;
	my @a=split/\s+/;
	my $query=$a[0];
	my $len=$a[1];
	if($query=~/part/){
		my ($contig,$tail)=(split/_part_/,$query)[0,1];
		my ($st,$en)=(split/\_/,$tail)[0,1];
		$sam_match{$contig}{$st}=$len;
	}else{
		my ($contig,$tail)=(split/_total_/,$query)[0,1];
		$sam_match{$contig}{1}=$len;
	}
}
close BLAST;

open BLAST, $blastmap =~/gz$/ ? "zcat $blastmap|" : "<$blastmap" or die $!;
while (<BLAST>) {
	chomp;
	my @a=split/\s+/;
	my $query=$a[0];
	my $chr=$a[1];
	my $contig_st=$a[6];
	my $contig_en=$a[7];
		#if($contig_st>$contig_en){
		#	my $tmp=$contig_st;
		#	$contig_st=$contig_en;
		#	$contig_en=$tmp;
		#}
	my $st=$a[8];
	my $en=$a[9];
	my $len=abs($en-$st)+1;
	if($query=~/part/){
		my ($contig,$tail)=(split/_part_/,$query)[0,1];
		my ($stsub,$ensub)=(split/\_/,$tail)[0,1];
		my $contig_st=$stsub+$contig_st-1;
		my $contig_en=$stsub+$contig_en-1;
		$match{$contig}{$chr}{$st}=0;
		$match{$contig}{$chr}{$en}=0;
		$nucmer{$contig}{$contig_st}{$contig_en}="$chr\t$st\t$en\t$len\tblast";
	}else{
		my ($contig,$tail)=(split/_total_/,$query)[0,1];
		$match{$contig}{$chr}{$st}=0;
		$match{$contig}{$chr}{$en}=0;
		$nucmer{$contig}{$contig_st}{$contig_en}="$chr\t$st\t$en\t$len\tblast";
	}
}
close BLAST;

my %rates;
foreach my $contig (sort keys %querylen) {
	my $lenall=$querylen{$contig};
	my $lensum=0;
	if(exists $sam_match{$contig}){
		foreach my $loc (sort keys %{$sam_match{$contig}}) {
			$lensum+=$sam_match{$contig}{$loc};
		}
	}
	my $rate=$lensum/$lenall;
	$rates{$contig}="$rate\t$lensum\t$lenall";
#	if(exists $match{$contig}){
#		foreach my $chr (sort keys %{$match{$contig}}) {
#			my @loc=sort {$a<=>$b} keys %{$match{$contig}{$chr}};
#			my $len=$loc[-1]-$loc[0]+1;
#			print $OUT "1\t$contig\t$lensum\t$lenall\t$rate\t$chr\t$loc[0]\t$loc[-1]\t$len\n";
#		}
#	}else{
#		print $OUT "0\t$contig\t$lensum\t$lenall\t$rate\tNA\tNA\tNA\tNA\n";
#	}
}

foreach my $contig (sort keys %nucmer) {
	foreach my $st (sort {$a<=>$b} keys %{$nucmer{$contig}}) {
		foreach my $en (sort {$a<=>$b} keys %{$nucmer{$contig}{$st}} ) {
			my $len1=abs($en-$st)+1;
			print $OUT "$contig\t$st\t$en\t$len1\t$nucmer{$contig}{$st}{$en}\t$rates{$contig}\n";
		}
	}
}
close $OUT;





 