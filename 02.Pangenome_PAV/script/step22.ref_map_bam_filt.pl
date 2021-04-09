use strict;

my $bam=shift;
my $bamout=shift;
 

open BAMOUT,"|/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -b - >$bamout" or die $!;

print "$bam start\n";
open my $MATCH,"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h $bam|" or die $!;
my %readmap;
while (<$MATCH>) {
	if($_=~/^@/){
		print BAMOUT $_;
		next;
	}
	chomp;
	my @a=split/\s+/,$_,11;
	my $read=$a[0];
	my %cigar;
	my $chr=$a[2];
	my $st=$a[3];
	my $qual=$a[4];
#	next if($qual<20);
	my $cigar=$a[5];
	my $seq=$a[9];
	$cigar{$cigar}{"$chr\t$st"}=0;
	my @str=split/\s+/,$a[10];
	foreach my $str (@str) {
		my @b=split/\:/,$str;
		if($str=~/^SA/){
			my @c=split/\;/,$b[2];
			foreach my $cc (@c) {
				my @d=split/\,/,$cc;
				my $chr1=$d[0];
				my $st1=$d[1];
				my $cigar1=$d[3];
				$cigar{$cigar}{"$chr1\t$st1"}=0;
			}
		}
		if($str=~/^XA/){
			my @c=split/\;/,$b[2];
			foreach my $cc (@c) {
				my @d=split/\,/,$cc;
				my $chr1=$d[0];
				my $st1=$d[1];
				my $cigar1=$d[2];
				$st1=reverse $st1;
				chop $st1;
				$st1=reverse $st1;
				$cigar{$cigar}{"$chr1\t$st1"}=0;
			}
		}
	}
	foreach my $cigar (sort keys %cigar) {
		my $match=&matchcigar($cigar);
		my $rate=$match/150;
		if($rate>0.6){
			print BAMOUT "$_\n";
			last;
		}
	}
}
close $MATCH;

print "output start\n";
close BAMOUT;
 

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
