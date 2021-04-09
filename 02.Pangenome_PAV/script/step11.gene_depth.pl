use strict;
my $gff=shift;
my $depthfile=shift;
my $output=shift;

my %genest;
my %genename;
my %genescale;
open my $GFF,$gff=~/gz$/?"zcat $gff|":"<$gff" or die $!;
while (<$GFF>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $type=$a[2];
	my $st=$a[3];
	my $en=$a[4];
	my $strand=$a[6];
	my $feat=$a[8];
	next if($type ne "mRNA");
	$chr=~s/chromosome_//g;
	my @feat=split/\;/,$a[8];
	my $gene;
	foreach my $sub (@feat) {
		my @str=split/\=/,$sub;
		if($str[0] eq "ID"){
			$gene=$str[1];
			last;
		}
	}
	$genest{$chr}{$st}=$en;
	$genename{$chr}{$st}=$gene;
	$genescale{$chr}{$gene}="$st\t$en";
}
close $GFF;
print "$gff is over\n";


my $part=100;
my %chrbins;
my %scale;
foreach my $chr (sort keys %genest) {
	my $genenum=scalar(keys %{$genest{$chr}});
	my $num=int($genenum/$part)+1;
	my $index=0;
	my $bin=1;
	my @loc;
	foreach my $st (sort {$a<=>$b} keys %{$genest{$chr}}) {
		$index++;
		if($index%$num==0){
			push @loc,$st;
			foreach my $loc (@loc) {
				$chrbins{$chr}{$bin}{$loc}=$genest{$chr}{$loc};
			}
			$scale{$chr}{$bin}="$loc[0]\t$loc[-1]";
			$bin++;
			@loc=();
		}else{
			push @loc,$st;
		}
	}
	foreach my $loc (@loc) {
		$chrbins{$chr}{$bin}{$loc}=$genest{$chr}{$loc};
	}
	$scale{$chr}{$bin}="$loc[0]\t$loc[-1]";
}




my $cnum=1;
my %genelocdb;
my %cluster;
my %genedepth;
my %tagchr;
my $genenumber=0;
open my $DF,$depthfile=~/gz$/?"zcat $depthfile|":"<$depthfile" or die $!;
while (<$DF>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $pos=$a[1];
	my $depth=$a[2];
	if($depth==0){next};
	if(!exists $tagchr{$chr}){
		$cnum=1;
		$tagchr{$chr}=0;
	}
	if(!exists $scale{$chr}{$cnum}){
		next;
	}
	my ($wst,$wen)=(split/\s+/,$scale{$chr}{$cnum})[0,1];
	$wen=$genest{$chr}{$wen};
	#############  db ³ÉÁ¢
	if($pos>$wen){
		$cnum++;
		($wst,$wen)=(split/\s+/,$scale{$chr}{$cnum})[0,1];
		$wen=$genest{$chr}{$wen};
	}
	if(!exists $cluster{$chr}{$cnum}){
		%genelocdb=();
		my @st=sort {$a<=>$b} keys %{$chrbins{$chr}{$cnum}};
		$genenumber=scalar(@st);
		print "$chr\t$cnum\t$genenumber\t$wst\t$wen\t".localtime()."\n";
		foreach my $st (@st) {
			my $en=$genest{$chr}{$st};
			my $gene=$genename{$chr}{$st};
			for (my $loc=$st;$loc<=$en ;$loc++) {
				$genelocdb{$chr}{$loc}{$gene}=$en;
			}
		}
		$cluster{$chr}{$cnum}=0;
	}
	if(exists $genelocdb{$chr}{$pos}){
		foreach my $genet (sort keys %{$genelocdb{$chr}{$pos}}) {
			$genedepth{$chr}{$genet}{1}++;
			$genedepth{$chr}{$genet}{2}+=$depth;
		#	my $mean=$genedepth{$chr}{$genet}{2}/$genedepth{$chr}{$genet}{1};
		#	my ($st,$en)=(split/\s+/,$genescale{$chr}{$genet})[0,1];
		#	my $len=$en-$st+1;
		#	print "$chr\t$cnum\t$genenumber\t$wst\t$wen\t$genet\t$pos\t$genescale{$chr}{$genet}\t$len\t$mean\t$depth\t$genedepth{$chr}{$genet}{2}\t$genedepth{$chr}{$genet}{1}\n";
		}
	}
}
close $DF;




open my $OUT,">$output" or die $!;
foreach my $chr (sort keys %genedepth) {
	foreach my $gene (sort keys %{$genedepth{$chr}}) {
		my $meandepth=$genedepth{$chr}{$gene}{2}/$genedepth{$chr}{$gene}{1};
		my ($st,$en)=(split/\s+/,$genescale{$chr}{$gene})[0,1];
		my $len=$en-$st+1;
		my $coverage=$genedepth{$chr}{$gene}{1}/$len;
		print $OUT "$chr\t$gene\t$meandepth\t$coverage\t$genedepth{$chr}{$gene}{1}\t$genedepth{$chr}{$gene}{2}\n";
	}
}
close $OUT;
