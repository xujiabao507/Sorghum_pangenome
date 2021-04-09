use strict;
my $map=shift;
my $fasta=shift;
my $output=shift;

open my $OUT,">$output" or die $!;
my %mapseq_total;
my %mapseq_part;
my %mapseq_partd;
open my $MAP,$map=~/gz$/?"zcat $map|":"<$map" or die $!;
while (<$MAP>) {
	chomp;
	my @a=split/\s+/;
	my $loc1=$a[3];
	my $loc2=$a[4];
	my $min=$loc1;
	my $max=$loc2;
	if($min>$loc2){
		$min=$loc2;
		$max=$loc1;
	}
	my $len=$max-$min+1;
	if($len<5000){
		next;
	}
	my $direction=$a[18];
	my $sca=$a[20];
	my $mrate=$a[21];
	if($mrate<20){next};
#	if($sca ne "Contig2418_pilon"){next};
	if($mrate>60){
		$mapseq_total{$sca}=0;
	}else{
		$mapseq_part{$sca}{$min}=$max;
		$mapseq_partd{$sca}{$min}=$direction;
	#	print "0\t$sca\t$min\t$max\t$mrate\n";
	}
}
close $MAP;



open FA, $fasta =~/gz$/ ? "zcat $fasta|" : "<$fasta" or die $!;
$/=">";
<FA>;
$/="\n";
while(<FA>){
	chomp;
	my $name=(split/\s+/,$_)[0];
	$/=">";
	my $seq=<FA>;
	$/="\n";
	$seq=~s/\s+//g;
	$seq=~s/\>//g;
	my $lenseq=length($seq);
	#if($name ne "Contig2418_pilon"){next};
	if(!exists $mapseq_total{$name} and !exists $mapseq_part{$name}){
		print $OUT ">$name\_total\_$lenseq\n$seq\n";
		#print "$name\n";
		next;
	}
	if(exists $mapseq_part{$name}){
		my %unmap;
		my @loc=sort {$a<=>$b} keys %{$mapseq_part{$name}};
		my $st=1;
		my $en=0;
		for (my $i=0;$i<@loc ;$i++) {
			$en=$loc[$i]-1;
			if($en>$st){
				my $len=abs($en-$st)+1;
				if($len>1000){
					$unmap{$st}=$en;
				}
		#		print "1\t$i\t$name\t$st\t$en\t$loc[$i]\t$len\n";
			}
			$st=$mapseq_part{$name}{$loc[$i]}+1;
		#	print "1.5\t$i\t$name\t$st\t$en\n";
		}
		$en=length($seq);
		#print "1\tend\t$name\t$st\t$en\n";
		if($en>$st){
			$unmap{$st}=$en;
		}
		foreach my $st1 (sort {$a<=>$b} keys %unmap) {
			my $en1=$unmap{$st1};
			my $len=$en1-$st1+1;
			my $substr=substr($seq,$st1-1,$len);
			print $OUT ">$name\_part\_$st1\_$en1\_$len\n$substr\n";
		#	print   "2\t$name\t$st1\t$en1\t$len\n";
		}
	}
}
close FA;
close $OUT;
















