use strict;
my $delmat=shift;
my $out=shift;
my %dis;
&dis();

my $num1=0;
my %dismat;
my %dismatnum;
open DEL,"zcat $delmat|" or die $!;
print "$delmat\n";
my $head=<DEL>;
my @head=split/\s+/,$head,2;
my @sam=split/\s+/,$head[1];
while (<DEL>) {
	chomp;
	my @a=split/\s+/,$_,2;
	my @val=split/\s+/,$a[1];
	for (my $i=0;$i<@val ;$i++) {
		for (my $j=$i+1;$j<@val ;$j++) {
			my $sam1=$sam[$i];
			my $sam2=$sam[$j];
			my $gt1=$val[$i];
			my $gt2=$val[$j];
			$dismat{$sam1}{$sam2}+=$dis{$gt1}{$gt2};
			$dismat{$sam2}{$sam1}+=$dis{$gt2}{$gt1};
			$dismatnum{$sam1}{$sam2}++;
			$dismatnum{$sam2}{$sam1}++;
		}
	}
	$num1++;
	if($num1%1000==0){
		print "$num1\t$_\n";
	}
}
close DEL;


my @samall=sort keys %dismat;
my $samnum=scalar(@samall);
print "OutputTree:\t$out\n";
open OUT,">$out" or die $!;
print OUT "$samnum\n";
foreach my $sam1 (@samall) {
	my @val;
	foreach my $sam2 (@samall) {
		if($sam1 eq $sam2){
			push @val,"0.00";
		}else{
			my $dis=$dismat{$sam1}{$sam2}/$dismatnum{$sam1}{$sam2};
			push @val,$dis;
		}
	}
	my $valstr=join "\t",@val;
	my $samtag=$sam1;
	#if(length($sam1)>0){
	#	$samtag=substr($sam1,0,7);
	#}
	#my $rest=" " x (10-length($samtag));
	print OUT "$samtag\t$valstr\n";
}
close OUT;

sub dis{
	$dis{0}{0}=0;
	$dis{0}{1}=1;
	$dis{1}{0}=1;
	$dis{1}{1}=0;
}












