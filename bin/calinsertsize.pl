#!/usr/bin/perl
use 5.010;

if(scalar @ARGV != 1){
	say STDERR "Usage: histogram <ins.info.file>";
	exit(0);
}

while(<>){
	($elemins, $flag) = split;
	#say $ins;
	if(/all_proper/){
		push @ins,$elemins;
	}elsif(/all_not_proper/){
		push @ins,$elemins;
	}elsif(/not_all_proper/){
		push @ins,$elemins;
	}elsif(/#/ or /^\s*$/){
		next;
	}else{
		say STDERR "invalid record: $_";
		last;
	}
}

@ins = sort {$a <=> $b} @ins;
say "Totally, there are ".(scalar @ins)." pairs with insert sizes";

$szelems = int( 0.99 * (scalar @ins) );
$szelems = ($szelems < 0.99 * (scalar @ins)) ? $szelems + 1: $szelems;
say "99% of size = ".$szelems;

$insSpan = $ins[-1] - $ins[0];
$lspanindex = 0;
$rspanindex = $#ins;

for( $szelems - 1..$#ins ){
	if( ($ins[$_] - $ins[$_ - $szelems + 1]) < $insSpan ){
		$insSpan = $ins[$_] - $ins[$_ - $szelems + 1];
		$lspanindex = $_ - $szelems + 1;
		$rspanindex = $_;
	}
}

say "insert span, the min zone contain 99% elems: min = $ins[$lspanindex]\tmax = $ins[$rspanindex]";

$top2pop = $#ins - $rspanindex;
$bottom2shift = $lspanindex;

@ins95 = @ins;

while($top2pop--){
	pop @ins95;
}

while($bottom2shift--){
	shift @ins95;
}

$sum = 0;
$sum += $_ for @ins95;
$avg = $sum/(scalar @ins95);
say "average insert size = $avg";

$qsum = 0;
$qsum += ($_ - $avg)*($_ - $avg) for @ins95;
$diff = $qsum/$#ins95;
$std = sqrt $diff;
say "standard deviation = $std";

$cnt = 1;
for(1..$#ins){
	if($ins[$_] == $ins[$_ - 1]){
		$cnt++;
	}else{
		print "$ins[$_ - 1]\t$cnt\n";
		$cnt = 1;
	}
}
say "";
