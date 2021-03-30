use strict;

my $bwamem=shift;
my $fasta=shift;
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
	if($meanlen<5000){next};
	$maps{$scf}{$scfpos}="$chr\t$pos";
}
close SAM;
print "$bwamem is over\n";

open OUT,">$out" or die $!;
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
	if(exists $maps{$chr}){next};
	print OUT ">$chr\n$seq\n";
	#$seqs{$chr}=$seq;
}
close FA;
close OUT;


