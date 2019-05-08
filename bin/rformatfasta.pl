#!/usr/bin/perl
use 5.010;

$cnt = -1;
$contig = "";
while(<>){
	$_ = <>;
	$cnt++;
	say ">$cnt";
	while(length $_> 60){
		say substr($_,0,60);
		$_ = substr($_,60);
	}
}

open RNM,">mapping/reads_num.info" or die $!;
say RNM $cnt;
