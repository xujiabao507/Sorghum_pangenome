use strict;
my $fasta=shift;
my $sv=shift;
my $out=shift;


my %faseq;
open FASTA, $fasta =~/gz$/ ? "zcat $fasta|" : "<$fasta" or die $!;
$/=">";
<FASTA>;
$/="\n";
while(<FASTA>){
	chomp;
	my $name=(split/\s+/,$_)[0];
	$/=">";
	my $seq=<FASTA>;
	$/="\n";
	$seq=~s/\s+//g;
	$seq=~s/\>//g;
	$faseq{$name}=$seq;
}
close FASTA;

my %svs;
open SV,$sv or die $!;
while (<SV>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[7];
	my $st=$a[8];
	my $en=$a[9];
	my $min=$st;
	my $max=$en;
	if($st<=$en){
		$max=$max+500;
	}
	if($st>$en){
		$min=$en;
		$max=$st;
		$min=$min-500;
	}
	$svs{$chr}{$min}=$max;
}
close SV;

open OUT,">$out" or die $!;
foreach my $chr (sort keys %svs) {
	my $seq=$faseq{$chr};
	foreach my $st (sort {$a<=>$b} keys %{$svs{$chr}}) {
		my $en=$svs{$chr}{$st};
		my $len=$en-$st;
		my $subseq=substr($seq,$st-1,$len);
		my $qual="J" x $len;
	#	print OUT ">ref_$chr\_$st\_$en\n$subseq\n";
		print OUT "\@ref_$chr\_$st\_$en\n$subseq\n+\n$qual\n";
	}
}
close OUT;














