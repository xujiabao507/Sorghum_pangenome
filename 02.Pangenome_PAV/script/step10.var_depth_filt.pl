use strict;
my $varfile=shift;
my $output=shift;

my %snps;
my ($snp,$all);
my %depth;
my %snpdep;
my $IN;
open $IN,$varfile=~/gz$/?"zcat $varfile|":"<$varfile" or die $!;
while (<$IN>) {
        chomp;
        next if($_=~/^#/);
        my @a=split/\s+/;
        my $chr=$a[0];
        my $pos=$a[1];
        my $ref=$a[3];
        my $alt=$a[4];
		my $depthsum;
		my @feat=split/\;/,$a[7];
		foreach my $feat (@feat) {
			my @str=split/\=/,$feat;
			if($str[0] eq "DP"){
				$depthsum=$str[1];
				last;
			}
		}
		my @str=split/\:/,$a[9];
		my $sum_allele=0;
		my @da=split/\,/,$str[2];
		my $max=0;
		for (my $i=0;$i<@da ;$i++) {
			$sum_allele+=$da[$i];
			if($i>0){
				if($max<$da[$i]){
					$max=$da[$i];
				}
			}
		}
		if($sum_allele/$depthsum<0.8){next};
		if($max/$sum_allele>0.1 and $max>10){
			#my $head=join "\t",@a[0..6];
			#print "$head\tDP=$depthsum\t$a[9]\n";
			print "$_\n";
		}
}
close $IN;
 
