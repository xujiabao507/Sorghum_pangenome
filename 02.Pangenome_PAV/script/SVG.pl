use strict;  
use warnings;  
use lib "/zfssz5/BC_PUB/Software/03.Soft_ALL/PerlInfo/lib/perl5";
use SVG;

my $contig=shift;


my %contig;
my %cont2chr;
my %contigs;
my %chrs;
my %contscale;
my %chrscale;
open my $SEL,"/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S09.PAV/S353/blast/SBgenome.map_result" or die $!;
while (<$SEL>) {
	chomp;
	my @a=split/\s+/,$_;
	my $contig_sub=$a[0];
	my $cst=$a[1];
	my $cen=$a[2];
	my $cmid=($cst+$cen)/2;
	my $clen=$a[3];
	my $chr_in=$a[4];
	my $chrst=$a[5];
	my $chren=$a[6];
	my $chrmid=($chrst+$chren)/2;
	my $chrlen=$a[7];
	if($contig ne $contig_sub){next};
	$contig{$contig_sub}{$cmid}=$chrmid;
	$cont2chr{$contig_sub}=$chr_in;
	$contigs{$contig_sub}{$cmid}="$cst\t$cen";
	$chrs{$chr_in}{$chrmid}="$chrst\t$chren";
	$contscale{$cst}=0;
	$contscale{$cen}=0;
	$chrscale{$chrst}=0;
	$chrscale{$chren}=0;
}
close $SEL;


my %chrlen;
open my $SL,"/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S02.DB/SBgenome.fai" or die $!;
while (<$SL>) {
	chomp;
	my @a=split/\s+/;
	my $chr_in=$a[0];
	my $len=$a[1];
	$chrlen{$chr_in}=$len;
}
close $SL;

my %contiglen;
open   $SL,"/hwfssz4/BC_COM_P5/F17FTSNCKF1543/SORcsaD/S22.Pangenome/S02.DB/S353.genome.fa.fai" or die $!;
while (<$SL>) {
	chomp;
	my @a=split/\s+/;
	my $chr_in=$a[0];
	my $len=$a[1];
	$contiglen{$chr_in}=$len;
}
close $SL;



my $chr=$cont2chr{$contig};
my $chrlen=$chrlen{$chr};
my $clen=$contiglen{$contig};


my $width = 900;  
my $height = 900;#这儿是画布变量，真是画布可以加大，这儿只是为了方便下面是使用  
my $width1=800;
my $height1=800;

my $svg = SVG->new('width', $width, 'height', $height);  
my $x = $svg->group(id=>'group_x', style=>{stroke=>'black','stroke-width',1} );  
my $y = $svg->group(id=>'group_y', style=>{stroke=>'blue', 'stroke-width',0.5} );  
my $z = $svg->group(id=>'group_z', style=>{stroke=>'red', 'stroke-width',0.5} );  
my $p = $svg->group(id=>'group_p', style=>{stroke=>'purple', 'stroke-width',0.1} );  


my $base_x=50;
my $base_y=50;
my $h1=$height-2*$base_y;
my $h2=$height/2-2*$base_y;
my $meany=$h1/$chrlen;
my $meany2=$h2/$clen;
 
 
$x->line(x1=>$base_x, y1=>$base_y, x2=>$base_x, y2=>$base_y+$chrlen*$meany  );
for (my $i=0;$i*10000000<$chrlen;$i++) {
	my $loc=$i*10000000;
	my $x1=$base_x-10;
	my $x2=$base_x-5;
	my $y=$base_y+$loc*$meany;
	my $tag=int($loc/10000000);
	$x->line(x1=>$x1, y1=>$y, x2=>$x2, y2=>$y);
	$svg->text(x=>$x1-10, y=>$y, -cdata=>"$tag");
}

#
#	my @mid=sort {$a<=>$b} keys %{$chrs{$chr}};
#	my $chrx1=$base_x+($num-1)*$xinterval;
#	my $chry1=$base_y+$mmin*$meany;
#
#
#	$x->line(x1=>$chrx1, y1=>$chry1, x2=>$chrx1, y2=>$chry2 );#Chr左轴  
#	$svg->text(x=>$chrx1-20, y=>$chry1-15, -cdata=>"LG$num");
#	foreach my $dis (@dis) {
#		my $lx1=$chrx1-3;
#		my $lx2=$chrx2+3;
#		my $ly1=$base_y+$dis*$meany;
#		my $ly2=$base_y+$dis*$meany;
#		$y->line(x1=>$lx1, y1=>$ly1, x2=>$lx2, y2=>$ly2 );#Chr左轴
#	}
#	my $dismax=int($dis[-1]);
#	my $maxy=$base_y+$dismax*$meany+20;
#	$svg->text(x=>$chrx1-20, y=>$maxy+10, -cdata=>"$dismax CM");
#
#	my @bp=sort {$a<=>$b} keys %{$chr2scf{$num}};
#	my $mminbp=$bp[0];
#	my $mmaxbp=$bp[0];
#	my $chrx3=$base_x+($num-1)*$xinterval+100;
#	my $chrx4=$base_x+($num-1)*$xinterval+5+100;
#	my $chry3=$base_y+$mminbp*$meany2;
#	my $chry4=$base_y+$mmaxbp*$meany2;
#	$x->line(x1=>$chrx3, y1=>$chry3, x2=>$chrx3, y2=>$chry4 );#Chr左轴  
#	$x->line(x1=>$chrx4, y1=>$chry3, x2=>$chrx4, y2=>$chry4 );#Chr右轴
#	$svg->text(x=>$chrx3-20, y=>$chry1-15, -cdata=>"Chr$num");
#	foreach my $st (@bp) {
#		my $lx1=$chrx3-5;
#		my $lx2=$chrx4+5;
#		my $en=$chr2scf{$num}{$st};
#		my $ly1=$base_y+$st*$meany2;
#		my $ly2=$base_y+$en*$meany2;
#		my $wdth=10;
#		my $heth=($en-$st+1)*$meany2;
#		my $scf=$loc2scf{$st};
#		#if(!exists $strand{$scf}){
#		if(($en-$st+1)>2000000){
#			$svg->rectangle(x=>$lx1,y=>$ly1,width=>$wdth,height=>$heth,rx=>5,ry=>5,fill=>'red') #,stroke => 'red'
#		}else{
#			$svg->rectangle(x=>$lx1,y=>$ly1,width=>$wdth,height=>$heth,rx=>5,ry=>5,fill=>'orange') #,stroke => 'blue'
#		}
#	}
#	my $bpmax=int($bp[-1]/100000)/10;
#	my $bpmaxy=$base_y+$bpmax*$meany2*1000000+20;
#	$svg->text(x=>$chrx3-20, y=>$bpmaxy+10, -cdata=>"$bpmax Mb");
#
#	foreach my $dis (@dis) {
#		my $stx=$chrx2+4;
#		my $sty=$base_y+$dis*$meany;
#		my $enx=0;
#		my $eny=0;
#		my $flag=0;
#		my $scf2=0;
#		foreach my $scf (sort keys %{$dis2scf{$num}{$dis}}) {
#			if(exists $scfloc{$scf}){
#				my $loc2= $scfloc{$scf};
#				$enx=$chrx3-6;
#				$eny=$base_y+$loc2*$meany2;
#				$flag=1;
#				$scf2=$scf;
#			}
#		}
#		if($flag==0){
#			#print "$num\t$scf2\t$dis\n";
#			next;
#		}
#		$p->line(x1=>$stx, y1=>$sty, x2=>$enx, y2=>$eny );
#	}
#}

print $svg->xmlify;  
