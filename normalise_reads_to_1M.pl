#!/usr/bin/perl -w
use strict;
my $input = shift;
	if (!$input){
		print STDERR "$0: <input file - stacked fasta format>\n";		
		print STDERR "no file given - exiting\n";
		exit();
	}
my $total = 0;
open (INPUT, "$input") or die "Can't open input file $input\n";
while (<INPUT>){
    if ($_=~/^>(\w+)\((\S+)\)/){
        my $count = $2;
        chomp ($count);
        $total+=$count;
	}
}
close INPUT;

open (INPUT, "$input") or die "Can't open input file $input\n";
while (<INPUT>){
	if ($_=~/^>(\w+)\((\S+)\)/){
		my $seq = $1;
		my $count = $2;
		chomp ($seq, $count);
		my $norm = ($count/$total)*1000000;
		print ">$seq($norm)\n$seq\n";
	}
}
close INPUT;

# by Simon Moxon
