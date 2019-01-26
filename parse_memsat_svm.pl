#!/usr/bin/env perl
use strict;
use warnings;

open(F, $ARGV[0]) or die $!; 
my ($s, $t, $c, $n, $f);

while (<F>) {
	if (/^Signal peptide:\s+(.+)/) { 
		$s =$1; $s=$s =~ /^Not detected/ ? "False" : "True"; 
	}
	elsif (/^Topology:\s+(((\d+-\d+),?)+)/) {
		$t=$1;
	} 
	elsif (/^Helix count:\s+(\d+)/) {
		$c=$1;
	} 
	elsif (/^N-terminal:\s+(.+)$/) {
		$n=$1;
	}
} 
($f) = ($ARGV[0]=~ /([^\/]+)\/output\/[^\/]+.memsat_svm/); 

print join("\t", $f,"memsat_svm", $s, $c, $n, $t),"\n" if $t;  
