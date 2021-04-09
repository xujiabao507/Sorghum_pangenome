use strict;
my $fastafile=shift;
my $bam=shift;
my $output=shift;

print "$fastafile start\n";
open my $OUTMAP,"|gzip >$output.unmap.gz" or die $!;
open my $OUTMAP2,"|/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -S - -b >$output.map.bam" or die $!;
my %lens;
open CONTIG, $fastafile =~/gz$/ ? "zcat $fastafile|" : "<$fastafile" or die $!;
while(<CONTIG>){
	chomp;
	my @a=split/\s+/;
	my $name=$a[0];
	my $len=$a[1];
	$lens{$name}=$len;
}
close CONTIG;

print "$bam start\n";
open MATCH,"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h $bam|" or die $!;
while (<MATCH>) {
	chomp;
#	print "$_\n";
	my @a=split/\s+/;
	my $read=$a[0];
	my $cigar=$a[5];
	my $cigar_tag=$a[5];
	my $match=0;
	if($_=~/^@/){
		print $OUTMAP2 "$_\n";
		next;
	}	
	while (1) {
		if($cigar_tag=~/(\d{1,9})M/){
			$match+=$1;
			$cigar_tag=$';
		}else{
			last;
		}
	}
	my $rate=$match/$lens{$read};
	if($rate<0.60){
		print $OUTMAP "$read\t$rate\t$match\t$lens{$read}\t$cigar\t$_\n";
	}else{
		print $OUTMAP2 "$_\n";
	}
#	print  "$read\t$rate\t$match\t$lens{$read}\t$cigar\t$_\n";
}
close MATCH;
close $OUTMAP;
close $OUTMAP2;

#system("/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -S $output.map.sam -b -o $output.map.bam")
