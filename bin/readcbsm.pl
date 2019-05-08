#!/usr/bin/perl
use 5.010;

open BSE,"mapping/bases.fa" or die $!;
while(<BSE>){
	chomp;
	s/>//;
	$name = $_;
	$seq = <BSE>;
	chomp $seq;
	$hsh{$name} = length $seq;
}
close BSE;

open infile,"<mapping/bases.sam" or die $!;

while(<infile>){
	next if $_ eq "";
	next if /^@/;
	$totalreads1++;
	($seqname,$flag,$chrm,$pos,$lparam5,$lparam6,$lparam7,$lparam8,$lparam9,$lparam10) = split /\t/;
	
	$brevs = 0;
	$readsmapped1++ if !($flag & 0x4);
	$brevs = 1 if ($flag & 0x10);
	$endpos = $hsh{$seqname} + $pos - 1;
	if( !($flag & 0x4) ){
		#say "$seqname\t$brevs\t$chrm\t$pos\t$endpos";
		$hsh2{$seqname} = "$seqname\t$brevs\t$chrm\t$pos\t$endpos"
	}
}
close infile;

foreach(keys %hsh2){
	if(/:/){
		s/:\d//;
		$line0 = $hsh2{"$_:0"};
		$line1 = $hsh2{"$_:1"};
		($seq, $brevs, $ref, $pos1, $pos2) = split /\s/,$line0;
		($seq, $brevs, $ref, $pos3, $pos4) = split /\s/,$line1;
		if($brevs){
			$startpos = $pos3;
			$endpos = $pos2;
		}else{
			$startpos = $pos1;
			$endpos = $pos4;
		}
		
		$hsh3{$_} = "$_\t$brevs\t$ref\t$startpos\t$endpos";
	}else{
		$hsh3{$_} = $hsh2{$_};
	}
}

foreach(sort {$a<=>$b} keys %hsh3){
	say "$hsh3{$_}";
	push @arr,"$hsh3{$_}";
}

@arr = sort {(split /\s/,$a)[3] <=> (split /\s/,$b)[3]} @arr;

say "======================================";
$pre = -1;
foreach(@arr){
	($a, $b, $c, $d, $e) = split;
	if($pre == -1){
		$gap = 0;
		$pre = $e;
	}else{
		$gap = $d - $pre;
		$pre = $e;
	}
	print;
	say "\t|| $gap";
}
