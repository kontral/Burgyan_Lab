#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage = <<__EOUSAGE__;

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    5p_analyzator.pl: the advanced read length analyzer of omnitool
    ###############################################################

    Hi there! I analyse size distribution and 5 prime starting
    nucleotid of sRNA reads. Please enter any number of stacked
    fasta files.
    I output a single PDF file with the plots. I know they ara
    not too pretty so I output a plain text files with a count
    table for each input. So if you know a better way to draw them.
           
    -s            plot each nucleotide in a different column
                  (default is to plot all 5 prime starting 
		  nucleotide into one column)

    -m   NUM      set minimum of read length to accept, also make
                  column if there are no reads (requires -M)

    -M   NUM      set maximum of read length to accept, also make
                  column if there are no reads (requires -m)

    -T            remove table files in the end

    -P            keep temporary PDF files

    -R            keep temporary R script file

    -a   STRING   set working directory (default is: current 
                  work directory)


    Version 0.1 (alpha)
    I was written by Levente Kontra

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

__EOUSAGE__
    ;    

my %options;
getopts("o:sRPTm:M:a:",\%options);

if ($options{s}){
    $options{s} = ", beside=TRUE";
}
else {
    $options{s} ="";
}

unless($options{a}){
    $options{a} = "./";
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
    

mkdir "$options{a}/pdf_dir";
mkdir "$options{a}/tables_dir";
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
    my %t;
    my %a;
    if ($options{m}){
	my $j = $options{m};
	while ($j <= $options{M}){
	    $a{$j} = 0;
	    $j++;
	}
    }
    
    my %g;
    my %c;
    my %sizes;
    my %output;
    my $total = 0;
    my $j = $plot_elements_min;
    while ($j <= $plot_elements_max){
	$t{$j} = 0;
	$a{$j} = 0;
	$g{$j} = 0;
	$c{$j} = 0;
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
#a formating step
#print "size\tU\tC\tG\tA\n";
#now lets count the percents...  
    foreach my $size (sort keys %a){  
	$output{$size} .= "$size nt";
	chomp $size;
	my $count = $sizes{$size};
# T
	if (!$t{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $t_each = ($t{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$t_each";
	}
# C
	if (!$c{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $c_each = ($c{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$c_each";
	}
# G
	if (!$g{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $g_each = ($g{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$g_each";
	}
# A
	if (!$a{$size}){
	    $output{$size} .="\t0";
	}
	else{
	    my $a_each = ($a{$size}/$count)*$sizes{$size};
	    $output{$size} .= "\t$a_each";
	}
    }
    $input =~ s/.fa//g or die "Incorrect file ending:$input\n";
    $input =~ s/.*\///;
    open (TABLE, ">$options{a}/$input.table") or die "Can't create table output file\n";  
    print TABLE "sizes\tU\tC\tG\tA\n";
    foreach my $o (sort keys %output){
	print TABLE "$output{$o}\n";
    }
    close TABLE;
    &plot($input, $i, $plot_max);
    system "mv $options{a}/$input.table $options{a}/tables_dir/";
    system "mv $options{a}/$input.pdf $options{a}/pdf_dir";
}    
system "pdf_merger.py $options{a}/pdf_dir $options{a}/plots.pdf";        #<-- portability issue

unless ($options{P}){
    system "rm -r $options{a}/pdf_dir";
}
if ($options{T}){
    system "rm -r $options{a}/tables_dir";
}


#creating temporary R script file
sub plot{
    my $plot_input = shift;
    my $pdf_size = shift;
    my $max =shift;
    open (NEW, ">$options{a}/$plot_input.plot") or die "Can't create temp file $plot_input.plot\n";
    print NEW "data<-read.table(\"$options{a}/$plot_input.table\",sep=\"\\t\", header = TRUE, row.names = 1)\n";
    print NEW "pdf(\"$options{a}/$plot_input.pdf\", width = $pdf_size)\n";
    print NEW "barplot(t(as.matrix(data)), ylim=c(0,$max), main=\"$plot_input\", col=c(\"blue\",\"red\",\"green\",\"purple\")$options{s})\n";
    print NEW "legend(\"topright\", c(\"U\",\"C\",\"G\",\"A\"), fill=c(\"blue\",\"red\",\"green\",\"purple\"))\n";
    print NEW "dev.off()\n";
    print NEW "quit(\"no\")\n";
    system "R --vanilla --quiet --slave < $options{a}/$plot_input.plot >$options{a}/blabla.log";
    unless ($options{R}){
	unlink "$options{a}/$plot_input.plot";
    }
}
