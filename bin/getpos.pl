#!/usr/bin/perl
use 5.010;

$cnt = -1;
$contig = "";
$prename = "";
while(<>){
	if(/>/){
		chomp;
		s/>//;
		$seqname = $_;
		if($cnt < 0){
			$prename = $seqname;
			$cnt++;
			next;
		}else{
			$hsh{$prename} = $contig;
			$len = length $contig;
			#say "$prename\t$len";
			$hsh2{"$len:$prename"} = 1;
			$hshlen{$prename} = $len;
			$contig = "";
			$cnt++;
			$prename = $seqname;
		}
	}else{
		chomp;
		$contig = $contig."$_";
	}
}
$hsh{$prename} = $contig;
$len = length $contig;
$hshlen{$prename} = $len;
$len = length $contig;
$hsh2{"$len:$prename"} = 1;


$no = 0;
foreach(sort { (split /:/,$b)[0]<=>(split /:/,$a)[0]} keys %hsh2){
	s/[^\s:]+://;
	$hsh3{$_} = $no;
	$no++;
}

$prename = "";
$cnt = 0;
foreach(sort {length $hsh{$b} <=> length $hsh{$a} } keys %hsh){
	$hshno{$_} = $cnt;
	$cnt++;
	$prename = $seqname;
}

$prepos = 1;
foreach(sort {(split /:/,$a)[0] cmp (split /:/,$b)[0] or (split /:/,$a)[1]<=>(split /:/,$b)[1]} keys %hsh){
	print $hshno{$_};
	print "\t";
	print $hshlen{$_};
	print "\t";
	(undef,undef,$pos) = split /:/,$_;
	$gap = $pos - $prepos;
	(undef,undef,undef,$prepos) = split /:/,$_;
	$gap = 0 if($gap < -10000);
	print $gap;
	print "\t";
	say ;
}
