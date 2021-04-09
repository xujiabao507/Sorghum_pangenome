use strict;
my $bam=shift;
my $samtools=shift;

	open OBO,"$samtools view -h $bam|" or die $!;
	while (<OBO>) {
		if($_=~/^@/){
			print $_;
			next;
		}else{
			my @a=split/\s+/;
			my $qual=$a[4];
			if($qual>=20){
				print $_;
				next;
			}
		}
	}
	close OBO;