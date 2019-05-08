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

# $maxspan *= 2;

# open CTG,"$ARGV[0]" or die $!;
# open BASE,">mapping/bases2.fa" or die $!;
# $cnt = 0;
# $contig = "";

# while(<CTG>){
	# next if />/;
	# chomp;
	# if(length > 2*$maxspan){
	# #if(0){
		# print BASE ">$cnt:0\n";
		# $subctg = substr($_, 0 , $maxspan);
		# print BASE "$subctg\n";
		# print BASE ">$cnt:1\n";
		# $cnt++;
		# $subctg = substr($_, (-1)*$maxspan);
		# print BASE "$subctg\n";
	# }else{
		# print BASE ">$cnt\n";
		# $cnt++;
		# print BASE "$_\n";
	# }
# }
# close BASE;
# close CTG;

open LOG,">mapping/domp.sh" or die $!;

$command = "python ${dirs}bowtie2-2.2.1/bowtie2-build combine.fasta index > mapping/bowtie2-build-bases.log 2>&1 & \nwait\n";
`$command`;
say LOG $command;

open INP,"mapping/inputfiles.info" or die $!;

$line = <INP>;
($threads, undef, undef) = split /\s/,$line;
$cnt = 1;
while(<INP>){
	($type, $left, $right, $len, $ovlp) = split;
	if($type eq "mp"){
		$command = "perl ${dirs}bowtie2-2.2.1/bowtie2 -k 10 --no-head -x index --threads $threads -U $left -S mapping/single.2.$cnt.1.sam > mapping/singlesummary$cnt.1 2>&1 \nwait\n";
		say LOG $command;
		
		$command = "perl ${dirs}bowtie2-2.2.1/bowtie2 -k 10 --no-head -x index --threads $threads -U $right -S mapping/single.2.$cnt.2.sam > mapping/singlesummary$cnt.2 2>&1 \nwait\n";
		say LOG $command;
		
		$cnt++;
	}
}

`sh mapping/domp.sh`;


