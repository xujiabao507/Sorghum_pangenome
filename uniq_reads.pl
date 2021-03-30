use strict;
use File::Basename;
my $bam=shift;
my $fq1=shift;
my $fq2=shift;
my $out=shift;

print "$bam start\n";
my %reads;
open BAM,"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view $bam|" or die $!;
while(<BAM>){
	chomp;
	my @a=split/\s+/;
	my $read=$a[0];
	$reads{$read}++;
}
close BAM;

print "$out start\n";
open FQO1,"|gzip >$out\_1.fq.gz" or die $!;
open FQO2,"|gzip >$out\_2.fq.gz" or die $!;
open FQ1,$fq1=~/gz/?"zcat $fq1|":"<$fq1" or die $!;
open FQ2,$fq2=~/gz/?"zcat $fq2|":"<$fq2" or die $!;
while (my $name1=<FQ1>) {
	my $seq1=<FQ1>;
	my $str1=<FQ1>;
	my $qual1=<FQ1>;
	my $name2=<FQ2>;
	my $seq2=<FQ2>;
	my $str2=<FQ2>;
	my $qual2=<FQ2>;
	$name1=~/\@(.*)\s+(\d)/;
	my $query=$1;
	if(exists $reads{$query}){
		if($reads{$query} > 1){
			print FQO1 "$name1$seq1$str1$qual1";
			print FQO2 "$name2$seq2$str2$qual2";
		}
	}
}
close FQ1;
close FQ2;
close FQO1;
close FQO2;

