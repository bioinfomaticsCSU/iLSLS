#!/usr/bin/perl
use 5.010;

if(scalar @ARGV < 3){
	say "Usage:\n\tcountins <infile> <minins> <maxins> -mp\n";
	exit(0);
}

open infile,"<${ARGV[0]}" or die $!;

$min = $ARGV[1];
$max = $ARGV[2];
$readlen = 0;

$mp = 0;

if(scalar @ARGV == 4){
	$mp = 1;
	say STDERR "# library is mate pair";
}
$readsmapped = 0;
$totalreads = 0;
$pairmapped = 0;
$samestrand = 0;
$difchrm = 0;
$insfalse = 0;
$left = 1;
while(<infile>){
	next if $_ eq "";
	next if /^@/;
	$totalreads++;
	if( $readlen == 0 ){
		($preseqname,$flag1,$prechrm,$lparam4,$lparam5,$lparam6,$lparam7,$lparam8,$lparam9,$lparam10) = split /\t/;
		$readlen = length $lparam10;
	}
	if($left){
		($preseqname,$flag1,$prechrm,$lparam4,$lparam5,$lparam6,$lparam7,$lparam8,$lparam9) = split /\t/;
	}else{
		($seqname,$flag2,$chrm,$rparam4,$rparam5,$rparam6,$rparam7,$rparam8,$rparam9) = split /\t/;
	}
	if(!$left){
		if($preseqname ne $seqname){
			say "error: records from one pair should be neigbors!";
			last;
		};
		
		$readsmapped++ if !($flag1 & 0x4);
		$readsmapped++ if !($flag2 & 0x4);
		if( !($flag1 & 0x4) && !($flag2 & 0x4) ){ #两个都 map 上
			if( ($flag1 & 0x2) && !($flag2 & 0x2) ){
					say STDERR "strange 1"; # 左边说全合适地 map 上，但右边却说没有
					last;
			};
			
			if( !($flag1 & 0x2) && ($flag2 & 0x2) ){
					say STDERR "strange 2"; # 左边说未合适地 map 上，但右边却说有
					last;
			};
			
			if( ($flag1 & 0x2) && ($flag2 & 0x2) ){ #全说合适地 map 上
				$pairmapped++ ;
				if( !!($flag1 & 0x10) ^ $mp ){ # left read is in the front of right read
					$ins = $lparam4 - $rparam4;
				}else{
					$ins = $rparam4 - $lparam4;
				}
				$ins += $readlen;
				say "$ins\tall_proper";
				$insfalse ++ if( $ins < $min || $ins > $max);
			}else{
				if( ($flag1 & 0x10) == ($flag2 & 0x10) ){
					$samestrand++;
					say STDERR "$preseqname is samestrand"; #cannot calculate insert size
				}elsif($prechrm ne $chrm) {
					say STDERR "$preseqname is difchrm" ; #cannot calculate insert size
				}else{
					if( !!($flag1 & 0x10) ^ $mp ){ # left read is in the front of right read
						$ins = $lparam4 - $rparam4;
					}else{
						$ins = $rparam4 - $lparam4;
					}
					if( !($flag1 & 0x2) && !($flag2 & 0x2) ){
						$ins += $readlen;
						say "$ins\tall_not_proper"; #still output insert size to stdout
						say STDERR "$preseqname is all_not_proper";
					}else{
						$ins += $readlen;
						say "$ins\tnot_all_proper"; #still output insert size to stdout
						say STDERR "$preseqname is not_all_proper";
					};
					if( $ins <= $min || $ins >= $max){
						say STDERR "$preseqname is insfalse"; #cannot calculate insert size
					};
				};
			}
		};
	}
	
	$left = !$left;
}


say "# totalreads = $totalreads";
say "# readsmapped = $readsmapped";
say "";
say "# total pairs = ".($totalreads/2);
say "# pairmapped = $pairmapped";
say "# samestrand = $samestrand";
say "# insfalse = $insfalse";

exit(0);

# if($seqname =~ m/1663059/){
	# say "$preseqname,$flag1,$lparam4,$lparam5,$lparam6,$lparam7,$lparam8,$lparam9";
	# say "$seqname,$flag2,$rparam4,$rparam5,$rparam6,$rparam7,$rparam8,$rparam9";
	# say ( !!($flag1 & 0x10) ^ $mp );
	# say ($flag1 & 0x10);
	# say $mp;
	# say "$lparam4\t$rparam4";
# }