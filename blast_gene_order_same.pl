use strict;

my $blast=shift;
my $out=shift;

my %scf2chr;
my %genenum;
my %chrgene;
my %scfgeneindex;
my %chrgeneindex;
my %genelocs;
open BLAST,"$blast" or die $!;
while (<BLAST>) {
	chomp;
	my @a=split/\s+/,$_;
	my $indexscf=$a[0];
	my $gene=$a[1];
	my $scf=$a[2];
	my $scfst=$a[3];
	my $scfen=$a[4];
	my $chr=$a[5];
	my $chrst=$a[6];
	my $chren=$a[7];
	$genelocs{$gene}{$scf}="$scf\t$scfst\t$scfen\t$chr\t$chrst\t$chren";
	$scf2chr{$scf}{$chr}++;
	$genenum{$scf}++;
	$scfgeneindex{$scf}{$gene}=$indexscf;
	$chrgene{$scf}{$chr}{$chrst}=$gene;
}
close BLAST;

my %mostsc;
foreach my $scf (sort keys %scf2chr) {
	foreach my $chr (sort {$scf2chr{$scf}{$b}<=>$scf2chr{$scf}{$a}} keys %{$scf2chr{$scf}}) {
		my $rate=$scf2chr{$scf}{$chr}/$genenum{$scf};
		if($rate>0.6){
			$mostsc{$scf}{$chr}++;
		}else{
			last;
		}
	}
}
my %chrgeneindex;
foreach my $scf (sort keys %mostsc) {
	foreach my $chr (sort keys %{$mostsc{$scf}}) {
		my %genes;
		my $index=0;
		foreach my $chrst (sort {$a<=>$b} keys %{$chrgene{$scf}{$chr}}) {
			my $gene=$chrgene{$scf}{$chr}{$chrst};
			if(!exists $genes{$gene}){
				$index++;
				$genes{$gene}=0;
				$chrgeneindex{$scf}{$chr}{$index}=$gene;
			}
		}
	}
}
open OUT1 ,">$out.1" or die $!;
open OUT2 ,">$out.2" or die $!;
foreach my $scf (sort keys %mostsc) {
	#if($scf ne "scaffold103|size389880"){next};
	foreach my $chr (sort keys %{$mostsc{$scf}}) {
		my %orders;
		foreach my $index (sort {$a<=>$b} keys %{$chrgeneindex{$scf}{$chr}}) {
			my $gene=$chrgeneindex{$scf}{$chr}{$index};
			my $scfindex=$scfgeneindex{$scf}{$gene};
			print OUT1 "$gene\t$scf\t$scfindex\t$chr\t$index\n";
			$orders{$index}=$scfindex;
	#		print "1\t$scf\t$index\t$scfindex\n";
		}
		my $direct=1;
		my @index=sort {$a<=>$b} keys %orders;
		my $for=0;
		my $rev=0;
		for (my $i=0;$i<@index-1 ;$i++) {
	#		print "$scf\t$index[$i]\t$orders{$index[$i]}\n";
			if($orders{$index[$i]}>$orders{$index[$i+1]}){
				$rev++;
			}else{
				$for++;
			}
		}
		if($rev>$for){$direct=-1};
		my %indexsub;
	#	print "$scf\t$direct\t$rev\t$for\n";
		if($direct==1){
			my $st=0;
			for (my $i=0;$i<@index-3 ;$i++) {
				if($orders{$index[$i]}<$orders{$index[$i+1]} and $orders{$index[$i+1]}<$orders{$index[$i+2]} ){
					$st=$i;
					last;
				}
			}
			my $locs=$orders{$index[$st]};
			$indexsub{$index[$st]}=$locs;
			for (my $i=$st+1;$i<@index;$i++) {
				if($locs<$orders{$index[$i]}){
					$locs=$orders{$index[$i]};
					$indexsub{$index[$i]}=$orders{$index[$i]};
				}
			}
		}else{
			my $st=0;
			for (my $i=0;$i<@index-3 ;$i++) {
				if($orders{$index[$i]}>$orders{$index[$i+1]} and $orders{$index[$i+1]}>$orders{$index[$i+2]} ){
					$st=$i;
					last;
				}
			}
			#print "x0\t$st\t$index[$st]\t$orders{$index[$st]}\n";
			my $locs=$orders{$index[$st]};
			$indexsub{$index[$st]}=$locs;
			for (my $i=$st+1;$i<@index;$i++) {
			#	print "x1\t$st\t$i\tY1\t$index[$i]\t$locs\t$orders{$index[$i]}\n";
				if($locs>$orders{$index[$i]}){
					$locs=$orders{$index[$i]};
					$indexsub{$index[$i]}=$orders{$index[$i]};
					#print "x2\t$st\t$i\tY2\t$locs\t$orders{$index[$i]}\n";
				}
				
			}
		}
		foreach my $index (sort {$a<=>$b} keys %{$chrgeneindex{$scf}{$chr}}) {
			my $gene=$chrgeneindex{$scf}{$chr}{$index};
			my $scfindex=$scfgeneindex{$scf}{$gene};
			if(exists $indexsub{$index}){
				print OUT2 "1\t$gene\t$scf\t$scfindex\t$chr\t$index\t$direct\t$genelocs{$gene}{$scf}\n";
			}else{
				print OUT2 "0\t$gene\t$scf\t$scfindex\t$chr\t$index\t$direct\t$genelocs{$gene}{$scf}\n";
			}
		}
	}
}
close OUT1;
close OUT2;






