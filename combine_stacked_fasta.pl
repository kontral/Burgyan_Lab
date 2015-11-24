#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my @input = @ARGV;
if (!@input){
	print STDERR "$0: <input file - stacked fasta format>\n";
	print STDERR "No input file given - exiting\n";
	exit();
	}
my %seq_count;
foreach my $file (@input){
open (FILE, $file) or die "Can't open file $file\n";
	while (<FILE>){
		if ($_=~m/^>(\S+)\((\S+)\)/){
			my $sequence = $1;
			my $count = $2;
			chomp ($sequence, $count);
			$seq_count{$sequence}+=$count;
		}
	}
close FILE;
}
foreach my $element (keys %seq_count){
	chomp $element;
	my $count = $seq_count{$element};
	my $length = length($element);
    print ">$element($count)\n$element\n";
}

#by Simon Moxon
