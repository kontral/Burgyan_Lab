#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $fasta = shift;
my $patman = shift;

if ((!$fasta)||(!$patman)){
	print STDERR"$0: A script to detect 1 nt truncated miRNAs\n";
	print STDERR"$0 please enter stacked <fasta file> and a <PatMaN file> generated by using the entered fasta file\n";
	exit();
}

my %output;
my %seen;
my %full_length;

my $seqio_object = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
while (my $seq = $seqio_object->next_seq){
	my $id = $seq->id();
	my $sequence = $seq->seq();
	my $length = length ($sequence);
	$seen{$id}++;
	$full_length{$id} = $length;
}

open (INPUT, "$patman") or die "Couldn't open patman input file\n";
while (<INPUT>){
	if ($_=~m/^(\S+)\s*.*\t(\w+)\((\S+)\)\t(\d+)\t(\d+)\t(\S+)\t(\d+)/){
	my $name = $1;
	my $sequence = $2;
	my $abundance = $3;
	my $start = $4;
	my $end = $5;
	my $mismatch = $7;		
	my $truncated = $end+1;
	$sequence =~ tr/U/T/;
	$sequence =~ tr/a-z/A-Z/;
	chomp ($name, $start, $end, $abundance, $sequence);
	if ($mismatch == 0){
		unless ($seen{$sequence}){
			if (($start ==1)and($full_length{$name} == $end)){  
			$name .="-full_length";
			$output{$name} = $abundance;
			}
			elsif (($start ==1)and($full_length{$name} == $truncated)){
			$name .="-3p_truncated";
			$output{$name}=$abundance;
			}
			elsif (($start ==2)and($full_length{$name} == $end)){
			$name .="-5p_truncated";
			$output{$name}=$abundance;
			}
		}
	}
}
}
my %final;
foreach my $element (sort keys %full_length){
	my $mir_f =$element;
	$mir_f .="-full_length";
	if (!$output{$mir_f}){
		$output{$mir_f}=0;
	}
	$final{$element} = "$output{$mir_f}";

	my $mir3 =$element;
	$mir3 .="-3p_truncated";
	if (!$output{$mir3}){
		$output{$mir3}=0;
	}
	$final{$element} .= "\t$output{$mir3}";

	my $mir5 =$element;
	$mir5 .="-5p_truncated";
	if (!$output{$mir5}){
		$output{$mir5}=0;
	}
	$final{$element} .= "\t$output{$mir5}";
}

print "id\tfullength\t3'_truncated\t5'_truncated\n";
foreach my $f (sort keys %final){
		print "$f\t$final{$f}\n";
} 


#made by Levente Kontra
