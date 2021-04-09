my $fastafile=shift;
open CONTIG, $fastafile =~/gz$/ ? "zcat $fastafile|" : "<$fastafile" or die $!;
$/=">";
<CONTIG>;
$/="\n";
while(my $name=<CONTIG>){
        $name=(split/\s+/,$name)[0];
        $/=">";
        my $seq=<CONTIG>;
        $/="\n";
	#$name=~s/\@//g;
        chomp $name;
        chomp $seq;
        #print  "$name\n";
        my $len=length($seq);
        $lens{$name}=$len;
	print "$name\t$len\n";
   }
                                                close CONTIG;
