use strict;
my $fastafile=shift;
my $bamfile=shift;
my $outfile=shift;
my $rnabam=shift;
$rnabam||=0;
my %lens;
print "$fastafile start\n";
open FA,$fastafile=~/gz$/?"zcat $fastafile|":"<$fastafile" or die $!;
#$/=">";
#<FA>;
#$/="\n";
while(my $name=<FA>){
#	chomp $chr;
#	$chr=(split/\s+/,$chr)[0];
#	$/=">";
#	my $seq=<FA>;
#	$seq=uc($seq);
#	$seq=~s/[>\s+]//g;
#	$/="\n";
#	my $len=length($seq);
	my $seq=<FA>;<FA>;<FA>;
	chomp $name;
	chomp $seq;
	$name=~s/\@//g;
	my $len=length($seq);
	$lens{$name}=$len;
#	print "$name\t$len\n";
}
close FA;

print "$bamfile start\n";
my %maps;
open SAM,"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h $bamfile|" or die $!;
while (<SAM>) {
	if($_=~/^\@/){
	#	print OUT "$_";
	#	print "$_";
		next;
	}
	chomp;
	my @a=split/\s+/,$_,12;
	my @b=split/\s+/,$_,2;
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
	#next if($qual<20);
	#next if($mate eq "=");
	#next if($feats=~/XA:Z/ or $feats=~/SA:Z/);
	my $match=0;
	while (1) {
		if($cigar=~/(\d+)(M)/){
			$match+=$1;
			$cigar=$';
		}else{
			last;
		}
	}
	#next if($match/$length<0.6);
	my $rate=$match/$lens{$read};
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
	$maps{$b[0]}{$rate}="$b[0]\t$rate\t$match/$lens{$read}\t$b[1]";
	#print "$b[0]\t$rate\t$match/$lens{$read}\t$b[1]\n";
}
close SAM;
print "output:\t$outfile\n";
open OUT,"|gzip >$outfile" or die $!;
foreach my $read (sort keys %maps) {
	foreach my $rate  (sort {$b<=>$a} keys %{$maps{$read}}) {
		print OUT "$maps{$read}{$rate}\n";
		last;
	}
}
close OUT;

