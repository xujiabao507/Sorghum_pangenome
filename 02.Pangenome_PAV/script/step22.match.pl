use strict;
my $fastafile=shift;
my $bam=shift;
my $output=shift;

print "$fastafile start\n";
open my $OUTMAP,"|gzip >$output" or die $!;
my %lens;
open CONTIG, $fastafile =~/gz$/ ? "zcat $fastafile|" : "<$fastafile" or die $!;
while(my $name=<CONTIG>){
	my $seq=<CONTIG>;
	<CONTIG>;<CONTIG>;
	$name=~s/\@//g;
	chomp $name;
	chomp $seq;
	my $len=length($seq);
	$lens{$name}=$len;
	#print "$name\t$len\n";
}
close CONTIG;

print "$bam start\n";
#open my $MATCH,"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view $bam|" or die $!;
open MATCH, $bam =~/gz$/ ? "zcat $bam|" : "<$bam" or die $!;
while (<MATCH>) {
	chomp;
	my @a=split/\s+/;
	my $readstr=$a[0];
	my @read=split/\_/,$readstr;
	my $num=scalar(@read);
	my $read=join "_",@read[0..$num-3];
	my $st=$read[-2];
	my $en=$read[-1];

	my $cigar=$a[5];
	my $cigar_tag=$a[5];
	my $match=0;
	while (1) {
		if($cigar_tag=~/(\d{1,3})M/){
			$match+=$1;
			$cigar_tag=$';
		}else{
			last;
		}
	}
	my $rate=$match/$lens{$readstr};
	print $OUTMAP "$read\t$st\t$en\t$rate\t$match\t$lens{$read}\t$cigar\n";
}
close MATCH;
close $OUTMAP;


