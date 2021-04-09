use strict;

my $map_read_file=shift;
my $map_read_cov=shift;

 
open my $COV,"|gzip >$map_read_cov" or die $!;
open MATCH,$map_read_file=~/gz$/?"zcat $map_read_file|":"<$map_read_file" or die $!;
my %readcov;
my %read2ct;
my %readsold;
while (<MATCH>) {
	chomp;
	my @a=split/\s+/;
	my $contig=$a[0];
	my $pos=$a[1];
	my @reads=split/\,/,$a[7];
	my %readsnew;
	foreach my $read (@reads) {
		$readcov{$read}{$pos}=0;
		$read2ct{$read}=$contig;
		$readsnew{$read}=0;
	}
	foreach my $read (sort keys %readsold) {
		if(!exists $readsnew{$read}){
			my @scale=sort {$a<=>$b} keys %{$readcov{$read}};
			my $st=$scale[0];
			my $en=$scale[-1];
			my $len=$en-$st+1;
			my $contig=$read2ct{$read};
			print $COV "$read\t$contig\t$st\t$en\t$len\n";
			print   "$read\t$contig\t$st\t$en\t$len\n";
			delete $readcov{$read};
		}
	}
	%readsold=();
	%readsold=%readsnew;
}
close MATCH;
close $COV;
 