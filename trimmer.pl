#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage = <<__EOUSAGE__;

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-	
    
    trimmer.pl: the adaptor trimmer module of Omnitool
    ##################################################

    Hi there! Please input <fastq files> to remove reads without
    adapter seqs and trim adapter bases(cofigurable) from reads.
    Also, I can remove low complexity reads, reads that are
    too short or too long. The output is a stacked fasta file.
    I can handel multiple inputs. However in this case
    make sure that the input filenames end with .fq or .fastq.

    -o   FILE     designate an output file (only works
                  with a single input file) 

    -b   STRING   adapter sequances (default is TGGAATTC)

    -c   FILE     create a file with the total readcount

    -l   NUM      lower length limit (default is 16)

    -u   NUM      upper length limit (default is 28)

    -z            input is gzip-ed

    -a   DIR      Output into a generated file(s) in this
                  directory (generates names just like when
                  multiple input is given)

    Version 0.1 (alpha)
    I was written by Levente Kontra

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

__EOUSAGE__
    ;    


my %options;
getopts("a:o:c:b:l:u:z",\%options);

#seting defaults
unless ($options{b}){$options{b} = "TGGAATTC";}
unless ($options{l}){$options{l} = 16;}
unless ($options{u}){$options{u} = 28;}
unless ($options{a}){$options{a} = "./";}
   
my @input = @ARGV;
my $headder;
my $trim_out;
my $total_out;
my $perc_out;
#error handeling
unless (@input) {
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Missing input\n";
    die $usage;
}
if (($#ARGV)&&($options{o})){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Cannot handel -o and multiple input at the same time.\n";
    die $usage;
}
$SIG{'INT'} = sub { die "interrupted, I am aborting...\n"};

foreach my $input (@input){
    my %seqs;
    chomp $options{a};
    my $total_count = 0;
    my $trimmed_count = 0;
    my $sequences;
    if ($options{z}){
	$sequences = Bio::SeqIO->new(-format => 'fastq', -file => "/bin/gunzip -c $input |");
        $input=~ s/.gz//; 
    }
    else {
	$sequences = Bio::SeqIO->new(-format => 'fastq', -file => "$input");
    }
    while (my $seq = $sequences->next_seq ) {
	$total_count++;
	my $sequence = $seq->seq();
	my $id = $seq->id();
	if ($sequence =~m/^(\w+)($options{b}\w*)$/){
	    my $processed_seq = $1;
	    my $adaptor = $2;
	    my %bases;
	    my $proced_length = length $processed_seq;
	    if ($proced_length >= $options{l}){
		if($proced_length <=$options{u}){
		    $trimmed_count++;
		    foreach (split('',$processed_seq)) {
			$bases{$_} = 1 ;
		    }
		    my $bcount = keys %bases;
		    if ($bcount <3) { # low complexity sequence
			next;
		    }
		    else{
			$seqs{$processed_seq}++;
		    }
		}
	    }
	}
    }
    my $input_name = $input;
    $input_name  =~ s/.*\///;
    $input_name =~ s/\.\w*$//;
    $headder .= "\t$input_name";
    $trim_out .= "\t$trimmed_count";
    $total_out .= "\t$total_count";
    my @keys = sort { $seqs{$b} <=> $seqs{$a} } keys %seqs;
    
#formating output
    if (($#ARGV)||($options{a})){
	$input =~ s/^.*\///;
	$input =~ s/\.\w+$/.fa/;
	open (MULTI_OUT, ">", "$options{a}/$input") or die "Can't create output file: $input\n";
	foreach my $element (@keys){
	    my $count = $seqs{$element};
	    print MULTI_OUT ">$element\($count\)\n$element\n";
	}
	close MULTI_OUT;
    }	
    elsif (defined $options{o}){
	open (NEW, ">","$options{o}") or die "Can't create output file: $options{o}\n";
	foreach my $element (@keys){
	    my $count = $seqs{$element};
	    print NEW ">$element\($count\)\n$element\n";
	}
	close NEW;
    }
    
    else {
	foreach my $element (@keys){
	    my $count = $seqs{$element};
	    print ">$element\($count\)\n$element\n";
	}
    }
#    print STDERR "finnished with $input\n";
}

#count table
if (defined $options{c}){
    open (COUNTS, ">", "$options{c}") or die "Can't create output file: $options{c}\n";
    print COUNTS "file name:$headder\n";
    print COUNTS "raw readcount:$total_out\n";
    print COUNTS "trimmed read count:$trim_out\n";
    close COUNTS;
}

