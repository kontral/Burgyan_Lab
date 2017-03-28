#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage = <<__EOUSAGE__;

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    muti_hit_finder.pl: the multi genome matching read finder
                          of Omnitool
    #########################################################

    I analyse seqs present in two stacked fasta. Please input
    two stacked fasta files or list of them separated
    by an argument: ",". Like this:

"muti_hit_finder.pl in_1_pair1.fa in_2_pair1.fa ... in_N_pair1.fa , in1_pair2.fa in2_pair2.fa ... in_N_pair2.fa"

    It is important to enter the files in the same order before
    and after the comma! Also, please select at least one output type: 

    -F   output seqs that are only in the first file(s)
    
    -S   output seqs that are only in the second file(s)

    -B   output seqs that are in both files

    -C   output a count table of reads in both files


    I could change the file names if you would like me to:

    -f   STRING   output seqs that are only in the first file(s)
                  whit this label. Default is: "-only.fa"

    -s   STRING   output seqs that are only in the second file(s)
                  whit this label. Default is: "-only.fa"

    -b   STRING   output seqs that are in both files with this
                  label. Default is: "in-both.fa"
           
    -c   FILE     Create a count table file 
                  Default is: "./reads_in_both.txt"
    
    -n            Suppress header for the count table

    -a   STRING   set working directory (default is: current
                  work directory)


    Version 0.1 (alpha)
    I was written by Levente Kontra

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

__EOUSAGE__
    ;    

my %options;
getopts("FBCSf:s:b:c:na:", \%options);


# input process
my $k;
my @input1;
my @input2;

foreach my $separator (@ARGV){
    if (!$k){
	if ($separator ne ","){
	    push @input1, $separator;
	}
	else {
	    $k++;
	    next;
	}
    }
    else {
	push @input2, $separator;
    }
}


#failsafe 
if ($options{f}){$options{F} = 1;}
if ($options{s}){$options{S} = 1;}
if ($options{b}){$options{B} = 1;}
if ($options{c}){$options{C} = 1;}
if (!$options{a}){$options{a} = "./";}


if ((!@input1) || (!@input2)){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Missing input\n";
    die $usage;
}
if ((!$options{F}) && (!$options{B}) && (!$options{C}) && (!$options{S})){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "\nplease set at least one output.\n";    
    die $usage;
}
if (($options{n})&&(!$options{C})){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Only count tables can have a header. Option -n will not work without -C (or -c)\n";
    die $usage;
}


#seting defaults
if (!$options{f}){$options{f} = "-only.fa";}
if (!$options{s}){$options{s} = "-only.fa";}
if (!$options{b}){$options{b} = "in-both.fa";}
if (!$options{c}){$options{c} = "reads_in_both.txt";}

my $header = "sample name:";
my $count_table ="counts:";

# Read in input sequences to look for..
foreach my $input1 (@input1){
    my %both; 			
    my %seen;
    my %uniq2;
    my $input2 = shift @input2;
    my $total = 0;
    my $seqio_object = Bio::SeqIO->new(-file => $input1, '-format' => 'Fasta');
    while (my $seq = $seqio_object->next_seq){
	my $sequence = $seq->seq();
	my $id = $seq->id();
	if ($id =~m/\w+\((\S+)\)/){
	    my $abundance = $1;
	chomp ($sequence);
	$seen{$sequence}=$abundance;
	}
	else{
	    print STDERR "$input1 does not appear to be in the correct format - exiting\n";
	    die $usage;
	}
    }

# Read in sequence DB to count from..
    my $seqio_object2 = Bio::SeqIO->new(-file => $input2, '-format' => 'Fasta');
    while (my $seq = $seqio_object2->next_seq){
	my $sequence = $seq->seq();
	my $id = $seq->id();
	chomp ($sequence);
	if ($id =~m/\w+\((\S+)\)/){
	    my $abundance = $1;
	    $uniq2{$sequence} = $abundance;
	}  
	if ($seen{$sequence}) { 	
	    $both{$sequence} = $seen{$sequence};		
	    $total += $seen{$sequence};
	}				
    }

# creating output

    if ($options{C}){
	my $namer1 = $input1;
	$namer1 =~ s/.*\///;
	$namer1 =~ s/\.\w+$//;

	my $namer2 = $input2;
	$namer2 =~ s/.*\///;
	$namer2 =~ s/\.\w+$//;
	$header .= "\tin $namer1 and $namer2";
	$count_table .= "\t$total";
    }
    if ($options{F}){
	my @keys_uniq1 = sort { $seen{$b} <=> $seen{$a} } keys %seen;
	my $uniq_file1 = $input1;
	$uniq_file1 =~ s/.fa/$options{f}/g or die "Incorrect file ending:$input1\n";
	open (UNIQ1, ">", "$options{a}/$uniq_file1") or die "Can't create output file: $input1\n";
	foreach my $element (@keys_uniq1){
	    unless ($both{$element}){
		    print UNIQ1 ">$element\($seen{$element})\n$element\n";
	    }
	}
	close UNIQ1;
    }



    if ($options{S}){
	my @keys_uniq2 = sort { $uniq2{$b} <=> $uniq2{$a} } keys %uniq2;
	my $uniq_file2 = $input2;
	$uniq_file2 =~ s/.fa/$options{s}/g or die "Incorrect file ending:$input2\n";
	open (UNIQ2, ">","$options{a}/$uniq_file2") or die "Can't create output file: $input2\n";
	foreach my $element (@keys_uniq2){
	    unless ($both{$element}){
		    print UNIQ2 ">$element\($uniq2{$element})\n$element\n";
	    }
	}
	close UNIQ2;
    }




     
    if ($options{B}){
	my @keys_both = sort { $both{$b} <=> $both{$a} } keys %both;
	$input1=~ s/\w+.fa/$options{b}/g or die "Incorrect file ending:$input1\n";
	    open (MULTI_OUT, ">", "$options{a}$input1") or die "Can't create output file: $input1\n";
	foreach my $element (@keys_both){
	    print MULTI_OUT ">$element\($both{$element}\)\n$element\n";
	}
	close MULTI_OUT;
    }
}

if ($options{C}){
    open (SUM, ">", "$options{c}") or die "Can't create count summary file: $options{c}\n";
    unless ($options{n}){print SUM "$header\n";}
    print SUM "$count_table\n";
    close SUM;
}

