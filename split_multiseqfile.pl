#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  split_multiseqfile.pl
#
#        USAGE:  ./split_multiseqfile.pl 
#
#  DESCRIPTION: split a multiple sequence file in a given format into several files (file extension = *.format)
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (Daniel Lang), <daniel.lang@biologie.uni-freiburg.de>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  03/10/06 15:02:03 CET
#     REVISION:  ---
#===============================================================================

use strict;
use Bio::SeqIO;
use Getopt::Std;

my $USAGE = "USAGE: $0 -h | -f format -i seq_field -a [uc|lc] [-c chunk_size] -o outpath FILE\n\n";

my %opts;

getopts('a:i:o:f:c:h', \%opts);

$opts{o}	= `pwd` unless $opts{o};
$opts{f}	= 'fasta' unless $opts{f};
$opts{i}	= 'display_id' unless $opts{i};

die $USAGE if ($opts{'h'});

die "$USAGE $!\n" unless -e $ARGV[0] && -f $ARGV[0];
die "$USAGE $!\n" unless -e $opts{o} && -d $opts{o};
$opts{o}=~ s/\/$//;

die $USAGE. "-a must be one of: 'uc' or 'lc'!" if $opts{a} && $opts{a} !~ /uc|lc/i;
$opts{a} = lc($opts{a}) if $opts{a};

my $in = Bio::SeqIO->new(-file=>$ARGV[0], -format=>$opts{f});
my ($chunk,$out,$chunk_nr);
while (my $seq = $in->next_seq) {
	$chunk++;
	my $identifier;
	if (eval("\$identifier = \$seq->$opts{i}")) {
		if ($opts{a}) {
			eval("\$identifier = $opts{a}(\$identifier)");
			$seq->display_id($identifier);
		}
		if (exists $opts{c} and $opts{c} =~ /^\d+$/ and $chunk == 1) {
			$chunk_nr++;
			$out = Bio::SeqIO->new(-file=>">$opts{o}/$chunk_nr.$opts{f}");
		}
		elsif (! exists $opts{c} or $opts{c} !~ /^\d+$/) {
			$out = Bio::SeqIO->new(-file=>">$opts{o}/$identifier.$opts{f}");
		}
		if (exists $opts{c} and $chunk == $opts{c}) {
			$chunk = 0;
		}
		$out->write_seq($seq);
	}
	else {
		die $@;
	}
}
if (exists $opts{c} and $opts{c} =~ /^\d+$/) {
	print "$chunk_nr chunks written to $opts{o}\n";
}
