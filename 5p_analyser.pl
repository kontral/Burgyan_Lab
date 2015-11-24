#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $input = shift;
if (!$input){
	print STDERR "$0: a script to analys 5' ending of all sRNA size classes each\n";
	print STDERR "$0 <INPUT sRNA file in stacked fasta format>\n";
	exit();
}
my %t;
my %a;
my %g;
my %c;
my %sizes;
my %output;
my $total = 0;
my $seqio_object = Bio::SeqIO->new(-file => $input, '-format' => 'Fasta');
while (my $seq = $seqio_object->next_seq){
	my $sequence = $seq->seq();
	$sequence =~ tr/U/T/;
	$sequence =~ tr/a-z/A-Z/;
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
	my $end_5p = substr($sequence, 0, 1);
	if (($end_5p eq 'T')|| ($end_5p eq 'U')){
		$t{$length}+=$abundance;
		$c{$length}+=0;
		$g{$length}+=0;
		$a{$length}+=0;
	}
	elsif ($end_5p eq 'C'){
		$c{$length}+=$abundance;
		$t{$length}+=0;
		$g{$length}+=0;
		$a{$length}+=0;
	}
	elsif ($end_5p eq 'G'){
		$g{$length}+=$abundance;
		$t{$length}+=0;
		$c{$length}+=0;
		$a{$length}+=0;
	}
	elsif ($end_5p eq 'A'){
		$a{$length}+=$abundance;
		$t{$length}+=0;
		$c{$length}+=0;
		$g{$length}+=0;
	}
	$sizes{$length}+=$abundance;
}
print "size\tU\tC\tG\tA\n";
my $i = 0;
foreach my $size (sort keys %sizes){
	$output{$i} .= $size;
		chomp $size;
	my $count = $sizes{$size};
	my $tpercentage = 0;
	foreach my $tvalue ($t{$size}){
		if (!$tvalue){
			$tvalue = 0;
			$tvalue = sprintf("%.2f", $tvalue);
			$output{$i} .= "\t$tvalue";
		}
		else{
			my $tpercentage = ($tvalue/$count)*100;
			$tpercentage = sprintf("%.2f", $tpercentage);
			$output{$i} .= "\t$tpercentage";
		}
	}
		my $cpercentage = 0;
	foreach my $cvalue ($c{$size}){
		if (!$cvalue){
			$cvalue = 0;
			$cvalue = sprintf("%.2f", $cvalue);
			$output{$i} .= "\t$cvalue";
		}
		else{
			my $cpercentage = ($cvalue/$count)*100;
			$cpercentage = sprintf("%.2f", $cpercentage);
			$output{$i} .= "\t$cpercentage";
		}
	}
	my $gpercentage = 0;
	foreach my $gvalue ($g{$size}){
		if (!$gvalue){
		$gvalue = 0;
		$gvalue = sprintf("%.2f", $gvalue);
		$output{$i} .= "\t$gvalue";
		}
		else{
		my $gpercentage = ($gvalue/$count)*100;
		$gpercentage = sprintf("%.2f", $gpercentage);
		$output{$i} .= "\t$gpercentage";
		}
	}
	my $apercentage = 0;
	foreach my $avalue ($a{$size}){
		if (!$avalue){
			$avalue = 0;
			$avalue = sprintf("%.2f", $avalue);
			$output{$i} .= "\t$avalue";
		}
		else{
			my $apercentage = ($avalue/$count)*100;
			$apercentage = sprintf("%.2f", $apercentage);
			$output{$i} .= "\t$apercentage";
		}
	}
	$i++;
}
foreach my $o (sort keys %output){
	print "$output{$o}\n";
}

# moded from create_PSWM_from_sRNAs.pl by Simon Moxon by Levente Kontra
