use strict;
my $bamfile=shift;
my $outfile=shift;
my $rnabam=shift;
$rnabam||=0;

open OUT,"|/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h - -b -o $outfile" or die $!;
my %map1;
open SAM,"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h $bamfile|" or die $!;
while (<SAM>) {
	if($_=~/^\@/){
		print OUT "$_";
	#	print "$_";
		next;
	}
	my @a=split/\s+/,$_,12;
	my $read=$a[0];
	my $chr=$a[2];
	my $pos=$a[3];
	my $qual=$a[4];
	my $cigar=$a[5];
	my $mate=$a[6];
	my $seq=$a[9];
	my $feats=$a[11];
	my @feats=split/\:/,$feats;
	my $length=length($seq);
	next if($qual<20);
	next if($mate eq "=");
	next if($feats=~/XA:Z/ or $feats=~/SA:Z/);
	my $match=0;
	while (1) {
		if($cigar=~/(\d+)(M)/){
			$match+=$1;
			$cigar=$';
		}else{
			last;
		}
	}
	next if($match/$length<0.9);
	my $rate=$match/$length;
	################### for RNA bellow
	if($rnabam==1){
		my $hit=0;
		foreach my $feat (@feats) {
			my @feat=split/\:/,$feat;
			if($feat eq "NH"){
				$hit=$feat[2];
				last;
			}
		}
		next if($hit>1);
	}
	################### for RNA above
	print OUT "$_";
}
close SAM;
close OUT; 


