use strict;
my $contig=shift;
my $blast=shift;
my $output=shift;

open my $OUT,">$output" or die $!;
my %contig;
open CONTIG, $contig =~/gz$/ ? "zcat $contig|" : "<$contig" or die $!;
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
	$contig{$name}=$len;
}
close CONTIG;


my %blast;
open BLAST, $blast =~/gz$/ ? "zcat $blast|" : "<$blast" or die $!;
while (<BLAST>) {
	chomp;
	my @a=split/\s+/;
	my $query=$a[0];
	my $st=$a[6];
	my $en=$a[7];
	$blast{$query}{$st}=$en;
}
close BLAST;

foreach my $query (sort keys %contig) {
	my %query_len;
	if(exists $blast{$query}){
		foreach my $st (sort {$a<=>$b} keys %{$blast{$query}}) {
			my $en=$blast{$query}{$st};
			for (my $loc=$st;$loc<=$en ;$loc++) {
				$query_len{$loc}++;
			}
			
		}
	}
	my $lensub=scalar(keys %query_len);
	my $lenall=$contig{$query};
	my $rate=$lensub/$lenall;
	print $OUT "$query\t$lensub\t$lenall\t$rate\n";
}
close $OUT;












