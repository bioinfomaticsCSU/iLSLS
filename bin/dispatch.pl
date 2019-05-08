#!/usr/bin/perl
use 5.010;

open LIB,"lib.info" or die $!;
$str = <LIB>;
($filenum, $kmer) = split /\s/,$str;
while($filenum--){
	$str = <LIB>;
	($filename, $min, $max,) = split /\s/,$str;
}
$mean = ($min + $max)/2;
#say $mean;
open infile,"<mapping/collectedSets" or die $!;

while(<infile>){
	$line1 = $_;
	$line2 = <infile>;
	next if($line1 eq "" || $line2 eq "");
	($readid, $chrm1, $pos1, $brevs1) = split /\s/,$line1;
	($readid, $chrm2, $pos2, $brevs2) = split /\s/,$line2;
	next if($chrm1 eq $chrm2);
	
	#seq [no] 0:left flank, 1:right flank, 2:full contig
	if($chrm1 =~ m/:/){
		$chrm1 =~ s/:/ /;
	}else{
		$chrm1 = "$chrm1 2";
	}
	if($chrm2 =~ m/:/){
		$chrm2 =~ s/:/ /;
	}else{
		$chrm2 = "$chrm2 2";
	}
	
	$hsh1{$chrm1} = "$hsh1{$chrm1}$chrm1\t$pos1\t$brevs1\n$chrm2\t$pos2\t$brevs2\n";
	$hsh2{$chrm2} = "$hsh2{$chrm2}$chrm1\t$pos1\t$brevs1\n$chrm2\t$pos2\t$brevs2\n";
}

open DSP1, ">mapping/dsp1" or die $!;
open DSP2, ">mapping/dsp2" or die $!;
foreach(keys %hsh1){
	say DSP1;
	print DSP1 $hsh1{$_};
	say DSP1 "-1 -1";
}
foreach(keys %hsh2){
	say DSP2;
	print DSP2 $hsh2{$_};
	say DSP2 "-1 -1";
}

exit(0);
