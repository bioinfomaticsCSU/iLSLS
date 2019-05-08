#!/usr/bin/perl
use 5.010;

open CTG,"combine.fasta" or die $!;

open POS,"ctgpos.info" or die $!;

while(<CTG>){
	$ctg = <CTG>;
	chomp $ctg;
	push @ctgs, $ctg;
}

while(<POS>){
	if(/=======/){
		while(<POS>){
			
		}
	}
}