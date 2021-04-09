use strict;
my $varset=shift;

my %insertion;
my %deletion;
my %refnt;
my %str;
my %indelp;
my %snps;
my %querynt;
open my $SNPS,">$varset.snp.homo" or die $!;
open my $INDELS,">$varset.indel.homo" or die $!;
open OBO, $varset =~/gz$/ ? "zcat $varset|" : "<$varset" or die $!;
<OBO>;<OBO>;<OBO>;<OBO>;
while (<OBO>) {
	chomp;
	$_=~s/(^\s+|\s+$)//g;
	my @a=split/\s+/;
	my $refpos=$a[0];
	my $refnt=$a[1];
	my $querynt=$a[2];
	my $querypos=$a[3];
	my $dis=$a[4];
	my $strand=$a[7];
	my $ref=$a[8];
	my $query=$a[9];
	$str{$query}{$querypos}=$strand;
	if($query eq "Contig1001_pilon" and $querypos == 9847){
		print "$query\t$strand\n";
	}
	if($querynt ne "." and $refnt ne "."){
		$snps{$query}{$querypos}="$query\t$querypos\t$querynt\t$ref\t$refpos\t$refnt";
		#print "$dis\t$query\t$querypos\t$querynt\t$ref\t$refpos\t$refnt\n";
		next;
	}
	if($querynt eq "."){  ### Ïà¶Ôref¶øÑÔ
		$refnt{$ref}{$refpos}=$refnt;
		$deletion{$query}{$querypos}{$ref}{$refpos}=0;
		$indelp{$query}{$querypos}=0;
		next;
	}
	if($refnt eq "."){
		$insertion{$ref}{$refpos}{$query}{$querypos}=0;
		$querynt{$query}{$querypos}=$querynt;
		$indelp{$query}{$querypos}=0;
		next;
	}
}
close OBO;
 
my %indels;
foreach my $query (sort keys %deletion) {
	foreach my $qpos (sort {$a<=>$b} keys %{$deletion{$query}}) {
		foreach my $chr (sort keys %{$deletion{$query}{$qpos}}) {
			my $seq;
			my @pos=sort {$a<=>$b} keys %{$deletion{$query}{$qpos}{$chr}};
			foreach my $pos (@pos) {
				$seq.=$refnt{$chr}{$pos};
			}
			my $strand=$str{$query}{$qpos};
			if($strand eq "-1"){
				$seq=reverse($seq);
			}
			my $seq_scf= "N" x length($seq);
			$indels{$query}{$qpos}="$query\t$qpos\t$strand\t$chr\t$pos[0]\t$pos[-1]\tdeletion\t$seq_scf\t$seq\n";
		}
	}
}
foreach my $chr (sort keys %insertion) {
	foreach my $pos (sort {$a<=>$b} keys %{$insertion{$chr}}) {
		foreach my $scf(sort keys %{$insertion{$chr}{$pos}}) {
			my $seq;
			my @scfpos=sort {$a<=>$b} keys %{$insertion{$chr}{$pos}{$scf}};
			foreach my $scfpos (@scfpos) {
				$seq.=$querynt{$scf}{$scfpos};
			}
			my $strand=$str{$scf}{$scfpos[0]};
			#if($str{$scf} eq "-1"){
			#	$seq=reverse($seq);
			#}
			my $refseq="N" x length($seq);
			$indels{$scf}{$scfpos[0]}="$scf\t$scfpos[0]\t$strand\t$chr\t$pos\tNA\tinsertion\t$seq\t$refseq\n";
		}
	}
}
foreach my $scf (sort keys %snps) {
	foreach my $pos(sort {$a<=>$b} keys  %{$snps{$scf}}) {
		if(!exists $indelp{$scf}{$pos}){
			print $SNPS "$snps{$scf}{$pos}\n";
		}
	}
}
foreach my $scf (sort keys %indels) {
	foreach my $pos(sort {$a<=>$b} keys  %{$indels{$scf}}) {
		print $INDELS "$indels{$scf}{$pos}";
	}
}
close $INDELS;
close $SNPS;