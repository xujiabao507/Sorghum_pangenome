use strict;


my %genelen;
	my $fa="step04.cds.merge.largest.fa";
	### format 
	#>ORTHOMCL10052|16|AusTRCF317961|AusTRCF317961_CCG014398.1
	#ATGGCGGCGGAGGCAGGGAAGACGGCGGCGACGAAGAGCGGCGGCGGGCAGATCATGACGGTGGTGCGGGGCCTGGACGTGGCGCGGTACATGGGGCGGTGGTACGAGATCGCGTCGTTCCCGTCCTTCTTCCAGCCCCGGGACGGGCGGGACACGCGCGCCACCTACCGGCTGCTGGAGGACGGCGCCACGGTGCACGTCCTCAACGAGACGTGGAGCAAGGGGAAGCGCGACTACATCGAGGGCACCGCCTACAAGGCCGACCCGAACAGCGACGAGGCCAAGCTCAAGGTCAAGTTCTACCTCCCGCCCTTCCTCCCCGTCATCCCCGTCGTCGGCGACTACTGGGTGCTCTACGTCGACGACGACTACCAGTACGCGCTCGTCGGCGAGCCGCGCCGCAAGAACCTCTGGATTCTGTGCAGGAAGACGAGCATCGACGAGGAGGTGTACAACCAGCTGGTGGAGCGGGCAAAGGAGGAAGGCTACGACGTGAGCAAGCTGCACAGGACGCCGCAGGACGACCCGCCGCCGGAGAGCGACGCCGCGCCCACTGACACCAAGGGAGTATGGTGGTTCAAGTCGCTCTTTGGGAAATGA
	open FA, $fa =~/gz$/ ? "zcat $fa|" : "<$fa" or die $!;
	$/=">";
	<FA>;
	$/="\n";
	while(my $namestr=<FA>){
		chomp;
		my $gene=(split/\s+/,$namestr)[0];
		$/=">";
		my $seq=<FA>;
		$/="\n";
		$seq=~s/\s+//g;
		$seq=~s/\>//g;
		my $len=length($seq);
		my $fam=(split/\|/,$gene)[0];
		$genelen{$fam}{$len}=0;
	}
	close FA;
 
my $file="step12.cds.Fam.SNP.xls";  ########### SNP Frequency . the 8th and 9th coloum is the Frequency for allele in 6th and 7th coloum
#ORTHOMCL1007   16      2793    1       20      A       G       11      1
#ORTHOMCL1007   16      2793    2       89      A       T       15      1
#ORTHOMCL1007   16      2793    3       91      C       T       15      1
#ORTHOMCL1007   16      2793    4       120     A       G       15      1
#ORTHOMCL1007   16      2793    5       136     C       T       15      1
#ORTHOMCL1007   16      2793    6       152     A       G       6       10
#ORTHOMCL1007   16      2793    7       158     A       G       1       15
#ORTHOMCL1007   16      2793    8       165     A       G       1       15
#ORTHOMCL1007   16      2793    9       249     A       G       14      2
#ORTHOMCL1007   16      2793    10      254     C       G       14      2
#ORTHOMCL1007   16      2793    11      257     A       C       2       14
#ORTHOMCL1007   16      2793    12      260     G       T       14      2
my %div;
open FL,$file or die $!;
while (<FL>) {
	chomp;
	my @a=split/\s+/;
	my $fam=$a[0];
	my $samnum=$a[1];
	my $len=$a[2];
	my $num1=$a[7]*2;
	my $num2=$a[8]*2;
	my $sum=$num1+$num2;
	my $div=2*$num1*$num2/($sum*($sum-1));   ### ¦Ð(theta waterson) diversity 
	$div{$fam}{1}+=$div;
	$div{$fam}{2}=$len;
	$div{$fam}{3}++;
	$div{$fam}{4}=$samnum;
}
close FL;
my %fam2;
foreach my $fam (sort keys %div) {
	my $mean=$div{$fam}{1}/$div{$fam}{2};
	print "$fam\t$mean\t$div{$fam}{1}\t$div{$fam}{3}\t$div{$fam}{2}\n";
}
 
foreach my $fam (sort keys %genelen) {
	my @len=sort {$b<=>$a} keys %{$genelen{$fam}};
	if(!exists $div{$fam}){
		print "$fam\t0\t0\t0\t$len[0]\n";
	}
}



