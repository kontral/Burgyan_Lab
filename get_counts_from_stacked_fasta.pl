#!/usr/bin/perl -w
use strict;
my $input = shift;
	if (!$input){
	print STDERR "$0: <input file - stacked fasta format>\n";
	print STDERR "No input file given - exiting\n";
	exit();
	}
my $total_seqs=0;
open (INPUT, "$input") or die "Can't open input file $input\n";
while (<INPUT>){
	if ($_=~m/^>\S+\((\S+)\)/){
	my $count = $1;
	$total_seqs+=$count;
	}
}
print "$total_seqs Sequences in $input\n";

# by Simon Moxon
