#!/usr/bin/perl
use 5.010;

open POS,"ctgpos.info" or die $!;

while(<POS>){
	($ctg) = split ;
	$hshpre{$ctg} = $pre;
	$hshlast{$pre} = $ctg;
	$pre = $ctg;
}

while(<>){
	$line = <>;
	chomp $line;
	for(split /\s/,$line){
		s/\+//;
		s/\-//;
		$expected1 = $hshpre{$_};
		$expected2 = $hshlast{$_};
		if($expected1 != $observed && $expected2 != $observed){
			say "$observed\t$_";
		}
		$observed = $_;
	}
}