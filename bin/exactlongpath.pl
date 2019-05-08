#!/usr/bin/perl
use File::Spec;
my $path_curf = File::Spec->rel2abs(__FILE__);
my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
#print "C Dir = ", $dirs,"\n";

open LIB,"lib.info" or die $!;

$line = <LIB>;
$maxspan = 0;
while(<LIB>){
	($libname, $min, $max, $std, $len, $ovlp) = split;
	$mean = ($min + $max)/2;
	if($mean > 1000){
		$maxspan = ($mean + 3*$std) > $maxspan?($mean + 3*$std):$maxspan;
	};
}

$maxspan *= 2;

open CTG,"$ARGV[0]" or die $!;
open BASE,">mapping/bases.fa" or die $!;
$cnt = 0;
while(<CTG>){
	next if />/;
	chomp;
	if(length > 2*$maxspan){
		print BASE ">$cnt:0\n";
		$subctg = substr($_, 0 , $maxspan);
		print BASE "$subctg\n";
		print BASE ">$cnt:1\n";
		$cnt++;
		$subctg = substr($_, (-1)*$maxspan);
		print BASE "$subctg\n";
	}else{
		print BASE ">$cnt\n";
		$cnt++;
		print BASE "$_\n";
	}
}
close BASE;
close CTG;



