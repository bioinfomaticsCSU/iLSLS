#!/usr/bin/perl
use 5.010;

open INFILES,"mapping/inputfiles.info" or die $!;

$libline = <INFILES>;
($thread, $kmer, $filenum) = split /\s/,$libline;
say "$filenum\t$kmer";


$filenum = 1;
foreach(@ARGV){
	open in,"$_" or die $!;
	$histline = <in>;
	$histline = <in>;
	$histline = <in>;
	$histline = <in>;
	(undef,undef,undef,undef,$mean) = split /\s/,$histline;
	$mean = abs int $mean;
	$histline = <in>;
	(undef,undef,undef,$std) = split /\s/,$histline;
	$std = abs int $std;
	
	$libline = <INFILES>;
	($e1, $e2, $e3, $e4,$e5) = split /\s/,$libline;
	
	open SE, "se_readins_$filenum.txt" or open SE, "se_$filenum.txt" or die $!;
	$mean -= $e4;
	$read = <SE>;
	chomp $read;
	$e4 = length $read;
	close SE;
	$mean += $e4;
	
	if($mean>1000){
		$minspan = int ($mean - 2 * $std);
		$maxspan = int ($mean + 2 * $std);
	}else{
		$minspan = $mean - int (1.5 * $std);
		$maxspan = $mean + int (1.5 * $std);
	}
	
	push @records, "$mean\t$minspan\t$maxspan\t$std\t$e4\t$e5";
	$filenum++;
};

@records = sort {(split $a)[0] <=> (split $b)[0]} @records;

$filenum = 1;
foreach(@records){
	($e1, $e2, $e3, $e4, $e5, $e6) = split /\s/,$_;
	say "lib$filenum.fasta\t$e2\t$e3\t$e4\t$e5\t$e6";
	$filenum++;
	
}



