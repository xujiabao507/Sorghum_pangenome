use strict;
my $fastafile =shift;
my $bam=shift;
my $bamout=shift;
my $mapout=shift;
my $unmapout=shift;

open BAMOUT,"|/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -b - >$bamout" or die $!;
open MAP,"|gzip >$mapout" or die $!;
open UNMAP,">$unmapout" or die $!;

print "$fastafile start\n";
my %lens;
my $num=0;
open CONTIG, $fastafile =~/gz$/ ? "zcat $fastafile|" : "<$fastafile" or die $!;
while(my $name=<CONTIG>){
	my $seq=<CONTIG>;
	my $tmp1=<CONTIG>;
	my $tmp2=<CONTIG>;
	$num+=4;
	$name=~s/\@//g;
	chomp $name;
	$name=~s/\@//g;
	chomp $seq;
	my $len=length($seq);
	$lens{$name}=$len;
	#if($name=~/\+/){
	#	print "$num\t$name\t$seq\t$_\n";
	#}
	chomp $tmp1;
	chomp $tmp2;
#	last;
#	print "$num\t$name\t$seq\t$tmp1\t$tmp2\n";
}
close CONTIG;

print "$bam start\n";
open my $MATCH,"/zfssz5/BC_PUB/Software/03.Soft_ALL/samtools-1.7/samtools view -h $bam|" or die $!;
my %readmap;
my %readrate;
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
	my $cigar=$a[5];
	my $cigarfirst=$cigar;
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
				$cigar{$cigar1}{"$chr1\t$st1"}=0;
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
				$cigar{$cigar1}{"$chr1\t$st1"}=0;
			}
		}
	}
	my ($chrread)=(split/\_/,$read)[0];
	#my %map;
	my %mapother;
	my $flag=0;
	my $match1=&matchcigar($cigarfirst);
	my $rate1=$match1/$lens{$read};
	$readrate{$read}="$rate1\t$cigarfirst\t$chr\t$st";
	if($rate1>0.6){
	#	$map{$chr}{$rate1}="$cigarfirst\t$chr\t$st\t$_";
		$mapother{$rate1}{$chr}="$cigarfirst\t$chr\t$st\t$_";
		$readmap{$read}=$rate1;
		print BAMOUT "$_\n";
	}else{
		foreach my $cigar (sort keys %cigar) {
			my $match=&matchcigar($cigar);
			my $rate=$match/$lens{$read};
			if($rate>0.6){
				$readmap{$read}=$rate;
				foreach my $loc (sort keys %{$cigar{$cigar}}) {
					my ($chrx,$pos)=(split/\s+/,$loc)[0,1];
	#				$map{$chrx}{$rate}="$cigar\t$chrx\t$pos\t$_";
					$mapother{$rate}{$chrx}="$cigar\t$chrx\t$pos\t$_";
				}
				$flag=1;
				print BAMOUT "$_\n";
				last;
			}
		}
	}
	foreach my $rate (sort {$b<=>$a} keys %mapother) {
		foreach my $chrother (sort keys %{$mapother{$rate}}) {
			print MAP "$read\t$mapother{$rate}{$chrother}\n";
		}
	}
}
close $MATCH;

print "output start\n";
close BAMOUT;



foreach my $read (sort keys %lens) {
	if(!exists $readmap{$read}){
		print UNMAP "$read\t$readrate{$read}\n";
	}
}
close MAP;
close UNMAP;

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
