use strict;
my $contig=shift;
my $output=shift;
my $outputdir=shift;

open CONTIG, $contig =~/gz$/ ? "zcat $contig|" : "<$contig" or die $!;
my $sum=0;
my $num=0;
my $str;
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
	###########################
	$sum+=length($seq);
	if($sum<10000000){
		$str.=">$name\n$seq\n";
	}else{
		$num++;
		my $tag="0" x (3-length($num))."$num";
		if(!-e "$outputdir/$output"){
			`mkdir -p "$outputdir/$output"`;
		}
		print "$outputdir/$output/$output.$tag.fasta\n";
		open OUT,">$outputdir/$output/$output.$tag.fasta";
		print OUT "$str";
		close OUT;
		$str="";
		$str.=">$name\n$seq\n";
		$sum=length($seq);
	}
}
$num++;
my $tag="0" x (3-length($num))."$num";
print "$outputdir/$output/$output.$tag.fasta\n";
open OUT,">$outputdir/$output/$output.$tag.fasta";
print OUT "$str";
close OUT;



close CONTIG;








