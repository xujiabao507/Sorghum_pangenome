use strict;
my $fastafile=shift;
my $output=shift;

open my $OUT,"|gzip >$output" or die $!;

my %fasta;
open CONTIG, $fastafile =~/gz$/ ? "zcat $fastafile|" : "<$fastafile" or die $!;
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
	my @seq=split//,$seq;
	if($name=~/^\d/){$name="chr_$name"};
	for (my $loc=1;$loc<$len ;$loc=$loc+100) {
		my $st=$loc-1;
		my $en=$st+500-1;
		my @sub=@seq[$st..$en];
		my $seqsub= join "",@sub;
		my $qual="J" x length($seqsub);
		#print "$name\t$st\t$en\n";
		print $OUT "\@$name\_$st\_$en\n$seqsub\n+\n$qual\n";
	}
}
close CONTIG;
close  $OUT;



