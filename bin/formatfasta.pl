#!/usr/bin/perl
use 5.010;

$cnt = -1;
$contig = "";
while(<>){
	if(/>/){
		if($cnt < 0){
			$cnt++;
			next;
		}else{
			say ">$cnt";
			say $contig;
			$contig = "";
			$cnt++;
		}
	}else{
		chomp;
		$contig = $contig."$_";
	}
}
say ">$cnt";
say $contig;

open RNM,">mapping/reads_num.info" or die $!;
say RNM $cnt;