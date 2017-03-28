#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $usage = <<__EOUSAGE__;

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    normalizator.pl: the read count normalizator of omnitools
    #########################################################

    Hi there! I normalize read counts from stacked fasta.
    Also, I can cut off low abundance reads, of course I count
    them first. Any number of fasta file can be inputed.

    -c          Suppress fasta output

    -O   FILE   designste an output file (only works 
		with a single input file) 

    -p   FILE   add another Nr. fasta file that will be
                counted for normalization

    -l   NUM    set cut off limit for low abundance (default is 10)

    -n   NUM    normalize to (default is 1 000 000)

    -a   DIR    Output into a generated file(s) in this
                directory (generates names just like when
                multiple input is given)

    -s          suppress header for tables 

    -f   FILE   create a table with the normalization factors

    -e   FILE   use external normalization factors
                (like the ones made with -f)


    Version 0.1 (alpha)

    I was written by Levente Kontra
-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
    
__EOUSAGE__
    ;    


my %options;
getopts("a:p:O:l:n:scf:e:",\%options);
my $nameing;
#seting defaults
unless ($options{l}){$options{l} = 10;}
unless ($options{a}){$options{a} = "";}
unless ($options{n}){$options{n} = 1000000;}


my @input = @ARGV;
#error handeling
unless (@input) {
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Missing input\n";
   die $usage;
}
if (($#ARGV)&&($options{O})){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Cannot handel -o and multiple input at the same time.\n";
    die $usage;
}
if (($options{c})&&($options{O})){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Cannot handel -c and -O at the same time.\n";
    die $usage;
}
if (($options{c})&&(!$options{f})){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "There is no point in running if there is no output.\n";
    die $usage;
}

if (($options{e})&&(($options{p})||($options{l})||($options{n})||($options{f}))){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "If using external normalization I cannot handel -p, -l, -n or -f.\n";
    die $usage;
}

$SIG{'INT'} = sub {die "interrupted, I am aborting...\n"};

my $header;
my $factor;
#processing
my $norm_string;
my @normalization;
if ($options{e}){
    open (NORM, $options{e}) or die "Can't open normalization table input file $options{n}\n";
    while (<NORM>){
        if ($_=~/^normalization factor:\t(.*)/){
            $norm_string = $1;
        }
    }

    @normalization = split (/\t/, $norm_string);
}


my $p_total;
if ($options{p}){
    open (ADDON, "$options{p}") or die "Can't open input file $options{p}\n";
    while (<ADDON>){
	if ($_=~/^>(\w+)\((\S+)\)/){
	    my $count = $2;
	    chomp ($count);
	    $p_total+=$count;
	}
    }
    close ADDON;
}

my $i;
foreach my $input (@input){
    my $namer = $input;
    $namer =~ s/.*\///;
    $namer =~s/\.\w+$//;
    $header .= "\t$namer";    
    my $total = 0;
    my %seen;
    open (INPUT, "$input") or die "Can't open input file $input\n";
    while (<INPUT>){
        if ($_=~/^>(\w+)\((\S+)\)/){
	    my $count = $2;
	    chomp ($count);
	    $total+=$count;
	}
    }
    close INPUT;
    if ($options{p}){
	$total += $p_total;
    }
    my $factor_count;
    if ($total != 0){
	$factor_count = (1/$total)*$options{n};
    }
    else{
	$factor_count = 0;
    }
    $factor .= "\t$factor_count";
    open (INPUT, "$input") or die "Can't open input file $input\n";
    while (<INPUT>){
	if ($_=~/^>(\w+)\((\S+)\)/){
	    my $seq = $1;
	    my $count = $2;
	    chomp ($seq, $count);
	    if ($count >= $options{l}){
		my $norm;
		if ($options{e}){
		$norm = $count * $normalization[$i];
		}
		else{
		$norm = ($count/$total)*$options{n};
		}
		$seen{$seq} = $norm;
	    }
	}
    }
    close INPUT;
    unless($options{c}){    
	my @keys = sort { $seen{$b} <=> $seen{$a} } keys %seen;
	
# outputing fasta    
	if (($#ARGV)||($options{a})){
	    $input=~ s/\.fa/-normed.fa/;
	    $input=~ s/^.*\///;
	    open (MULTI_OUT, ">", "$options{a}/$input") or die "Can't create output file: $input\n";
	    foreach my $element (@keys){
		print MULTI_OUT ">$element\($seen{$element}\)\n$element\n";
	    }
	    close MULTI_OUT;
	}
	elsif ($options{O}){
	    open (NEW, ">","$options{O}") or die "Can't create output file: $options{O}\n";
	    foreach my $element (@keys){
		print NEW ">$element\($seen{$element}\)\n$element\n";
	    }
	    close NEW;
	}
    
	else {
	    foreach my $element (@keys){
		print ">$element\($seen{$element}\)\n$element\n";
	    }
	}
    }
    $i++;
}
if ($options{f}){
    open (FACT, ">", "$options{f}") or die "Can't create normalizing factor table outputfile file: $options{f}\n";
    unless ($options{s}){
        print FACT "file name:$header\n";
    }
    print FACT "normalization factor:$factor\n";
    close FACT;
}

