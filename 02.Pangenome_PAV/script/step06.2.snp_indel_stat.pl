use strict;
my $dir=shift;
my $sam=shift;
my $stat=shift;

my $snp_hom="$dir/S04.MumFilt_dir/$sam/$sam.filt.snps.snp.homo";
my $snp_het="$dir/S07.GATK/$sam/$sam.snp_result.vcf.gz";
print "$snp_hom\n$snp_het\n";
my %snp;
open HOM, $snp_hom =~/gz$/ ? "zcat $snp_hom|" : "<$snp_hom" or die $!;
while (<HOM>) {
	chomp;
	my @a=split/\s+/;
	$snp{hom}++;
}
close HOM;

open HET, $snp_het =~/gz$/ ? "zcat $snp_het|" : "<$snp_het" or die $!;
<HET>;<HET>;
while (<HET>) {
	chomp;
	my @a=split/\s+/;
	$snp{het}++;
}
close HET;

my $indel_hom="$dir/S04.MumFilt_dir/$sam/$sam.filt.snps.indel.homo";
my $indel_het="$dir/S07.GATK/$sam/$sam.indel_result.vcf.gz";
print "$indel_hom\n$indel_het\n";


my %indel;
open HOM, $indel_hom =~/gz$/ ? "zcat $indel_hom|" : "<$indel_hom" or die $!;
while (<HOM>) {
	chomp;
	my @a=split/\s+/;
	$indel{hom}++;
}
close HOM;

open HET, $indel_het =~/gz$/ ? "zcat $indel_het|" : "<$indel_het" or die $!;
<HET>;<HET>;
while (<HET>) {
	chomp;
	my @a=split/\s+/;
	$indel{het}++;
}
close HET;

open OUT,">$stat" or die $!;
print OUT "Sam\tVar_Type\tSNP_hom\tSNP_het\n";
print OUT "$sam\tSNP\t$snp{hom}\t$snp{het}\n";
print OUT "$sam\tINDEL\t$indel{hom}\t$indel{het}\n";
close OUT;









