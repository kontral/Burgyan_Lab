#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $input = shift;
if (!$input){
	print STDERR "$0: <input file - stacked fasta format>\n";
	print STDERR "No input file given - exiting\n";
	exit();
	}
my %sizes;
my $total = 0;
my $seqio_object = Bio::SeqIO->new(-file => $input, '-format' => 'Fasta');
while (my $seq = $seqio_object->next_seq){
my $sequence = $seq->seq();
my $id = $seq->id();
my $length = length($sequence);
my $abundance = 0;
	if ($id =~m/^\S+\((\S+)\)/){
	$abundance = $1;
	}
	else{
	print STDERR "Sequence file does not appear to be in the correct format - exiting\n";
	exit();
	}
$sizes{$length}+=$abundance;
$total+=$abundance;
}
foreach my $element (sort keys %sizes){
my $count = $sizes{$element};
my $percentage = ($count/$total)*100;
$percentage = sprintf("%.2f", $percentage);
print "$element\t$count\n";
}

#made by Simon Moxon
