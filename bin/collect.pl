#!/usr/bin/perl
use 5.010;

open infile,"<mapping/single1.1.sam" or die $!;


$readsmapped1 = 0;
$totalreads1 = 0;

while(<infile>){
	next if $_ eq "";
	next if /^@/;
	$totalreads1++;
	($seqname,$flag,$chrm,$pos,$lparam5,$lparam6,$lparam7,$lparam8,$lparam9,$lparam10) = split /\t/;
	
	$brevs = 0;
	$readsmapped1++ if !($flag & 0x4);
	$brevs = 1 if ($flag & 0x10);
	if( !($flag & 0x4) ){
		$mparr1[$seqname] = "$chrm\t$pos\t$brevs";
	}
}

close infile;

say "#1 totalreads1 = $totalreads1";
say "#1 readsmapped1 = $readsmapped1";

open infile,"<mapping/single1.2.sam" or die $!;


$readsmapped2 = 0;
$totalreads2 = 0;

while(<infile>){
	next if $_ eq "";
	next if /^@/;
	$totalreads2++;
	($seqname,$flag,$chrm,$pos,$lparam5,$lparam6,$lparam7,$lparam8,$lparam9,$lparam10) = split /\t/;
	
	$brevs = 0;
	$readsmapped2++ if !($flag & 0x4);
	if( !($flag & 0x4) ){
		$brevs = 1 if ($flag & 0x10);
		$mparr2[$seqname] = "$chrm\t$pos\t$brevs";
	}
}

close infile;

say "#2 totalreads2 = $totalreads2";
say "#2 readsmapped2 = $readsmapped2";

say "#construct array done!";

open outfile,">mapping/collectedSets" or die $!;

foreach(0..$#mparr1){
	$line1 = $mparr1[$_];
	$line2 = $mparr2[$_];
	next if($line1 eq "" || $line2 eq "");
	
	($chrm1, $pos1, $brevs1) = split /\s/,$line1;
	($chrm2, $pos2, $brevs2) = split /\s/,$line2;
	#seq [no] 0:left flank, 1:right flank, 2:full contig
	if($chrm1 =~ m/:/){
		$chrm1 =~ s/:/\t/;
	}else{
		$chrm1 = "$chrm1\t2";
	}
	if($chrm2 =~ m/:/){
		$chrm2 =~ s/:/\t/;
	}else{
		$chrm2 = "$chrm2\t2";
	}
	next if($chrm1 eq $chrm2);
	say outfile "$chrm1\t$pos1\t$brevs1\n$chrm2\t$pos2\t$brevs2";
}

exit(0);
