#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage = <<__EOUSAGE__;

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    3p_analyzator.pl: the advanced read length analyzer of omnitool
    ###############################################################

    Hi there! I analyse size distribution and 5 prime starting
    nucleotid of sRNA reads. Please enter any number of stacked
    fasta files.
    I output is a single PDF file with the plots. I know they ara
    not too pretty so I output a plain text files with a count table 
    for each input. So if you know a better way to draw them ...
           
    -s   plot each nucleotide in a different column (default is to
	 plot all 5 prime starting nucleotide into one column)

    -m   set minimum of read length to accept, also make column 
         if there are no reads (requires -M)

    -M   set maximum of read length to accept, also make column
         if there are no reads (requires -m)

    -T   remove table files in the end

    -P   keep temporary PDF files

    -R   keep temporary R script file


    Version 0.1 (alpha)
    I was written by Levente Kontra

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

__EOUSAGE__
    ;    

my %options;
getopts("o:sRPTm:M:",\%options);

if ($options{s}){
    $options{s} = ", beside=TRUE";
}
else {
    $options{s} ="";
}


my @input = @ARGV;

if (!@input){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Missing input\n";
    die $usage;
}
if ((($options{m})&&(!$options{M}))||((!$options{m})&&($options{M}))){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "If using size limitaiton both -m or -M options are required\n";
    die $usage;
}


mkdir "./pdf_dir";
mkdir "./tables_dir";
my $plot_elements_min;
my $plot_elements_max;
my $plot_max = 0;
#unifying steps
foreach my $unify (@input){
    my %abundance_by_length;
    my $seqio_object = Bio::SeqIO->new(-file => $unify, '-format' => 'Fasta');
    while (my $seq = $seqio_object->next_seq){
	my $sequence = $seq->seq();
	my $id = $seq->id();
	my $length = length($sequence);
	if ($options{m}){
	    if ($length < $options{m}){next;}
	    if ($length > $options{M}){next;}
	}
	if (!$plot_elements_min){
	    $plot_elements_min = $length;
	    $plot_elements_max = $length;
	}
	elsif ($plot_elements_max < $length){
	    $plot_elements_max = $length;
	}
	elsif ($plot_elements_min > $length){
	    $plot_elements_min = $length;
	}
	if ($id =~m/^\S+\((\S+)\)/){
	    my $abundance = $1;
	    $abundance_by_length{$length} += $abundance;
	}
	else{
	    print STDERR "Sequence file does not appear to be in the correct format - exiting\n";
	    exit();
	}
    }
    foreach my $element (keys %abundance_by_length){
	if($plot_max < $abundance_by_length{$element}){
	    $plot_max = $abundance_by_length{$element};
	}
    }
}

$plot_max = $plot_max * 1.1;
my $i;
if ($options{m}){
    $i = $options{M} - $options{m};
}
else{
    $ i= $plot_elements_max - $plot_elements_min;
}

foreach my $input (@input){
    my %end_a;
    if ($options{m}){
	my $j = $options{m};
	while ($j <= $options{M}){
	    $end_a{$j} = 0;
	    $j++;
	}
    }
    my %end_t;
    my %end_g;
    my %end_c;
    my %second_a;
    my %second_t;
    my %second_g;
    my %second_c;
    my %sizes;
    my %output;
    my $total = 0;
    my $j = $plot_elements_min;
    while ($j <= $plot_elements_max){
	$end_t{$j} = 0;
	$end_a{$j} = 0;
	$end_g{$j} = 0;
	$end_c{$j} = 0;
	$second_t{$j} = 0;
	$second_a{$j} = 0;
	$second_g{$j} = 0;
	$second_c{$j} = 0;
	$j++;
    }
#reading in sRNA fasta
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
# geting 5' starting nucleotide counts for each length
	my $end_3p = substr($sequence, -1);
	if (($end_3p eq 'T')|| ($end_3p eq 'U')){
	    $end_t{$length}+=$abundance;
	    $end_c{$length}+=0;
	    $end_g{$length}+=0;
	    $end_a{$length}+=0;
	}
	elsif ($end_3p eq 'C'){
	    $end_c{$length}+=$abundance;
	    $end_t{$length}+=0;
	    $end_g{$length}+=0;
	    $end_a{$length}+=0;
	}
	elsif ($end_3p eq 'G'){
	    $end_g{$length}+=$abundance;
	    $end_t{$length}+=0;
	    $end_c{$length}+=0;
	    $end_a{$length}+=0;
	}
	elsif ($end_3p eq 'A'){
	    $end_a{$length}+=$abundance;
	    $end_t{$length}+=0;
	    $end_c{$length}+=0;
	    $end_g{$length}+=0;
	}
	my $second_3p = substr($sequence, -2, 1);
	if (($second_3p eq 'T')|| ($second_3p eq 'U')){
	    $second_t{$length}+=$abundance;
	    $second_c{$length}+=0;
	    $second_g{$length}+=0;
	    $second_a{$length}+=0;
	}
	elsif ($end_3p eq 'C'){
	    $second_c{$length}+=$abundance;
	    $second_t{$length}+=0;
	    $second_g{$length}+=0;
	    $second_a{$length}+=0;
	}
	elsif ($end_3p eq 'G'){
	    $second_g{$length}+=$abundance;
	    $second_t{$length}+=0;
	    $second_c{$length}+=0;
	    $second_a{$length}+=0;
	}
	elsif ($end_3p eq 'A'){
	    $second_a{$length}+=$abundance;
	    $second_t{$length}+=0;
	    $second_c{$length}+=0;
	    $second_g{$length}+=0;
	}

	$sizes{$length}+=$abundance;

	
	
    }
#a formating step
#print "size\tU\tC\tG\tA\n";
#now lets count the percents...  
    foreach my $size (sort keys %end_a){  
	$output{$size} .= "$size nt";
	chomp $size;
	my $count = $sizes{$size};
# T
	if (!$end_t{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $t_each = ($end_t{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$t_each";
	}
# C
	if (!$end_c{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $c_each = ($end_c{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$c_each";
	}
# G
	if (!$end_g{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $g_each = ($end_g{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$g_each";
	}
# A
	if (!$end_a{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $a_each = ($end_a{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$a_each";
	}
    }
    $input =~ s/\.fa.*$// or die "Incorrect file ending:$input\n";
    $input =~ s/.*\///;
    open (TABLE, ">$input.3p_end.table") or die "Can't create table output file\n";  
    print TABLE "sizes\tU\tC\tG\tA\n";
    foreach my $o (sort keys %output){
	print TABLE "$output{$o}\n";
    }
    close TABLE;
    &plot($input, $i, $plot_max);
    system "mv $input.3p_end.table ./tables_dir";
    system "mv $input.pdf ./pdf_dir";
}    
system "pdf_merger.py ./pdf_dir plots.pdf";        #<-- portability issue

unless ($options{P}){
    system "rm -r pdf_dir";
}
if ($options{T}){
    system "rm -r tables_dir";
}


#creating temporary R script file
sub plot{
    my $plot_input = shift;
    my $pdf_size = shift;
    my $max =shift;
    open (NEW, ">$plot_input.plot") or die "Can't create temp file $plot_input.plot\n";
    print NEW "data<-read.table(\"$plot_input.3p_end.table\",sep=\"\\t\", header = TRUE, row.names = 1)\n";
    print NEW "pdf(\"$plot_input.pdf\", width = $pdf_size)\n";
    print NEW "barplot(t(as.matrix(data)), ylim=c(0,$max), main=\"$plot_input\", col=c(\"blue\",\"red\",\"green\",\"purple\")$options{s})\n";
    print NEW "legend(\"topright\", c(\"U\",\"C\",\"G\",\"A\"), fill=c(\"blue\",\"red\",\"green\",\"purple\"))\n";
    print NEW "dev.off()\n";
    print NEW "quit(\"no\")\n";
    system "R --vanilla --quiet --slave < $plot_input.plot >blabla.log";
    unless ($options{R}){
	unlink "$plot_input.plot";
    }
}
