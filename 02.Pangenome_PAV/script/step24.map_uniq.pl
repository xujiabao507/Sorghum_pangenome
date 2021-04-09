use strict;
my $fastafile =shift;
my $bam=shift;
my $out=shift;

open OUT,">$out" or die $!;

print "$fastafile start \n";
my %length;
open CONTIG, $fastafile =~/gz$/ ? "zcat $fastafile|" : "<$fastafile" or die $!;
my $sum=0;
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
	$length{$name}=$len;
}
close CONTIG;

open BAM,$bam=~/bam$/?"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h $bam|":"<$bam" or die $!;
my %readmap;
while (<BAM>) {
	if($_=~/^@/){
		print BAMOUT $_;
		next;
	}
	chomp;
	my @a=split/\s+/,$_,11;
	my $read=$a[0];
	my $chr=$a[2];
	next if($chr eq "*");
	my $st=$a[3];
	my $qual=$a[4];
	my $cigar=$a[5];
	my $match=&matchcigar($cigar);
	#if(exists $length{$read} and $length{$read}<10){
	#	print "$read\n";
	#}
#	next;
	my $rate=$match/$length{$read};
	if($qual>10){
		print OUT "$read\t$chr\t$st\t$qual\t$match\t$length{$read}\t$rate\t$cigar\n";
	}
}
close BAM;
close OUT;

 


sub matchcigar{
	my ($cigar_sub)=@_;
	my $matchsub=0;
	while (1) {
		if($cigar_sub=~/(\d{1,3})M/){
			$matchsub+=$1;
			$cigar_sub=$';
		}else{
			last;
		}
	}
	return $matchsub;
}















