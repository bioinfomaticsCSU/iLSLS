#!/usr/bin/perl
use 5.010;

$cnt = 0;
while (<>) {
	if (/^@/) {
	  print "\@$cnt\n";
	  $_ = <>; print; $_ = <>; print;$_ = <>;print;
	  $cnt++;
	}
}

open RNM,">mapping/reads_num.info" or die $!;
print RNM $cnt;
