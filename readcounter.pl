#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage = <<__EOUSAGE__;

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    Omnitool:readcounter.pl
    ###############################################################

    Hi there! I count reads in a patman file. But for this I need a
    fasta file too. I do not coop well with multi fasta. If you want you
    could set the bucket size a.k.a. how many positins to be added
    into one value. I count both orientation.
           
    readcounter.pl (options) <patman input> <fasta input>

    -c            count covarage insted of 5p start position

    -b   NUM      set bucket size

    -s            Use output format 2. (Output only "+" strand,
		  each on to a new line.) 
  
    -a            Use output format 2. (Output only "-" strand,
		  each on to a new line.) 

    -m            Selects the first nucleotide to count from.

    -M            Selects the lest nucleotide to count to.

    Version 0.1 (alpha)
    I was written by Levente Kontra

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

__EOUSAGE__
    ;    

my %options;

#handeling input

getopts("cb:asm:M:",\%options);
my $patman = shift;
my $fasta = shift;

#safe gourd

if ((!$patman)||(!$fasta)){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Missing input\n";
    die $usage;
}
if (!$options{b}){$options{b}++;}
if(!$options{m}){$options{m}=0;}
else{$options{m}--;}

#reading in fasta
my %seq_lengths;
my $seqio_object = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
while (my $seq = $seqio_object->next_seq){
    my $id = $seq->id();
    my $length = $seq->length;
    $seq_lengths{$id}=$length;
    if (!$options{M}){
	$options{M} = $length;
    }
}
my %plus1;
my %minus1;

foreach my $element (keys %seq_lengths){
    my $length = $seq_lengths{$element};
    my $i = 0;
    while ($i < $length){
	push (@{$plus1{$element}}, "0");
	push (@{$minus1{$element}}, "0");
	$i++;
    }
}



open (INPUT1, "$patman") or die "Couldn't open patman input file\n";
while (<INPUT1>){
    if ($_=~m/^(\S+)\s*.*\t(\w+)\((\S+)\)\s+(\d+)\s+(\d+)\s+(\S+)\s+\d+/){
	my $name = $1;
	my $sequence = $2;
	my $abundance = $3;
	my $start = $4;
	my $end = $5;
	my $strand = $6;
	my $i = $start -1;
	my $j = $end -1;
	if ($options{c}){
	    if ($strand eq '+'){
		while ($i <= $end-1){
		    $plus1{$name}[$i]+=$abundance; 
		    $i++;
		}
	    }
	    else{
		while ($j >= $start-1){
		    $minus1{$name}[$j]+=$abundance; 
		    $j--;
		}
		
	    }
	}
	else {	
	    if ($strand eq '+'){
		$plus1{$name}[$i]+=$abundance; 
	    }
	    else{
		$minus1{$name}[$j]+=$abundance; 
	    }
	}
    }
    
}	
close INPUT1;


foreach my $locus (keys %plus1){
    my $id = $locus;
    my @plus_bucket;
    my @minus_bucket;
    my $i=1;
    my $j=0;
    my $k=0;
    
    foreach my $element (@{$plus1{$id}}){
	if ($i > $options{b}){
	    $i = 1;
	    $j++;
	}
	$plus_bucket[$j] +=$element;
	$minus_bucket[$j] +=$minus1{$id}[$k];
	$i++;
	$k++;
	}
    if ($options{s}){
	my $l=0;
	foreach my $position (@plus_bucket){
	    if (($l >$options{m})&&($l<=$options{M})){
		print "$position\n";
	    }
	    $l++;
	    if ($l >$options{M}){next;}
	}
    }
    elsif ($options{a}){
	my $l=0;
	foreach my $position (@minus_bucket){
	    if (($l >$options{m})&&($l<=$options{M})){
		print "$position\n";
	    }
	    $l++;
	    if ($l >$options{M}){next;}
	}
    }
    else{
	my $l=0;
	my $m=0;
	print "$id plus:\t";
	foreach my $position (@plus_bucket){
	    if (($l =$options{m})&&($l<=$options{M})){
		print "$position\t";
		$l++;
		if ($l >$options{M}){next;}
	    }
	}
	
	print "\n$id minus:\t";
	    foreach my $element (@minus_bucket){
		if (($m = $options{m})&&($m<=$options{M})){
		    print "$element\t";
		}
		$m++;
		if ($m >$options{M}){
		    next;
		}
	}
    }
}


