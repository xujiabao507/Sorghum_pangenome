use strict;

my $bwamem=shift;
my $fasta=shift;
my $geneline=shift;
my $out=shift;

my $base=(split/\//,$bwamem)[-1];
$base=(split/\./,$base)[0];

my %maps;
print "$bwamem is output\n";
open SAM,$bwamem=~/gz$/?"zcat $bwamem|":"<$bwamem" or die $!;
while (<SAM>) {
	chomp;
	my @a=split/\s+/,$_;
	my $scf=$a[0];
	my $scfpos=$a[4];
	my $chr=$a[2];
	my $pos=$a[3];
	my $meanlen=$a[9];
	my $rate=$a[12];
	if($rate<0.6){next};
	if($meanlen<1000){next};
#	if($meanlen<5000){next};
	my $int=int($pos/200000);
	$maps{$scf}{$scfpos}="$chr\t$pos";
#	print "$chr\t$pos\t$scf\t$scfpos\t$rate\t$meanlen\n";
}
close SAM;
my %midnum;
my %strand;
print "best is output\n";
foreach my $scf (sort keys %maps) {
	my @win=sort {$a<=>$b} keys %{$maps{$scf}};
	my $num=0;
	my $winnum=scalar(@win);
	my $midpos=0;
	my $chr;
	foreach my $scfpos (@win) {
		my ($chrsub,$pos)=(split/\s+/,$maps{$scf}{$scfpos})[0,1];
		$chr=$chrsub;
		$num++;
		if($num>=$winnum/2){
			$midpos=$pos;
			last;
		}
	}
	$midnum{$chr}{$midpos}=$scf;
	#############
	my $for=0;
	my $rev=0;
	for (my $i=0;$i<@win-1 ;$i++) {
		my $loc1=$win[$i];
		my $loc2=$win[$i+1];
		my ($chr1,$pos1)=(split/\s+/,$maps{$scf}{$loc1})[0,1];
		my ($chr2,$pos2)=(split/\s+/,$maps{$scf}{$loc2})[0,1];
		if($pos1<$pos2){$for++};
		if($pos1>$pos2){$rev++};
	}
	my $direction=1;
	if($rev>$for){$direction=-1};
	$strand{$scf}=$direction;
	#print "$chr\t$scf\t$direction\t$midpos\n";
}

my %typenum;
my %scales;
my %scfgenes;
my %scfgenesall;
open GL,"$geneline" or die $!;
while (<GL>) {
	chomp;
	my @a=split/\s+/,$_;
	my $type=$a[0];
	my $gene=$a[1];
	my $scf=$a[2];
	my $chr=$a[4];
	if($type==1){
		$typenum{$scf}++;
		my $chrst=$a[11];
		my $direct=$a[6];
		$scales{$scf}{$chr}{$chrst}=0;
		$strand{$scf}=$direct;
		$scfgenes{$scf}{$chr}{$gene}=0;
	}
	$scfgenesall{$scf}{$chr}{$gene}=0;
}
close GL;
my %genemap;
foreach my $scf (sort keys %scales) {
	foreach my $chr (sort keys %{$scales{$scf}}) {
		my @locs=sort {$a<=>$b} keys %{$scales{$scf}{$chr}};
		my @gene=sort keys %{$scfgenes{$scf}{$chr}};
		my @geneall=sort keys %{$scfgenesall{$scf}{$chr}};
		my $pairgene=scalar(@gene);
		my $pairgeneall=scalar(@geneall);
		my $rate=$pairgene/$pairgeneall;
		next if($pairgene<5);
		next if($rate<0.6);
		my $int=int(scalar(@locs)/2);
		my $mid=$locs[$int];
		$midnum{$chr}{$mid}=$scf;
		$genemap{$scf}="$pairgene\t$rate";
	}
}



my %seqs;
open FA,$fasta=~/gz$/?"zcat $fasta|":"<$fasta" or die $!;
$/=">";
<FA>;
$/="\n";
while(my $chr=<FA>){
	chomp $chr;
	$chr=(split/\s+/,$chr)[0];
	$/=">";
	my $seq=<FA>;
	$seq=uc($seq);
	$seq=~s/[>\s+]//g;
	$/="\n";
	$seqs{$chr}=$seq;
}
close FA;

open LEN,">$out.length.agp" or die $!;
open FA,">$out" or die $!;
my $sumbp=0;
my %contigs;
foreach my $chr (sort keys %midnum) {
	my $chrseq;
	my $st=1;
	my $en=0;
	foreach my $midpos (sort {$a<=>$b} keys %{$midnum{$chr}}) {
		my $contig=$midnum{$chr}{$midpos};
		$contigs{$contig}=0;
		my $seq=$seqs{$contig};
		my $len=length($seq);
		my $en=$st+$len-1;
		if(exists $genemap{$contig}){
			print LEN "$chr\t$contig\t$len\t$st\t$en\t$strand{$contig}\tGeneMap\t$genemap{$contig}\n";
		}else{
			print LEN "$chr\t$contig\t$len\t$st\t$en\t$strand{$contig}\tGenomeMap\t0\t0\n";
		}
		if($strand{$contig} eq "-1"){
			$seq=reverse($seq);
			$seq=~tr/[ACGT]/[TGCA]/;
		}
		my $str="N" x 100;
		$chrseq.="$seq$str";
		$st=$en+1+100-1+1;
		
	}
	$sumbp+=length($chrseq);
	print FA ">$chr\n$chrseq\n";
}
my $sumchr=$sumbp;
foreach my $contig (sort keys %seqs) {
	if(!exists $contigs{$contig}){
		$sumbp+=length($seqs{$contig});
		print FA ">$contig\n$seqs{$contig}\n";
	}
}
close FA;
my $gzv=$sumchr/$sumbp;
print LEN "the genome length $base is:\t$sumchr\t$sumbp\t$gzv\n";
print  "the genome length for $base is:\t$sumchr\t$sumbp\t$gzv\n";
close LEN;



















