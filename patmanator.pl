#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;

my $usage = <<__EOUSAGE__;

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    Patmanator.pl: the PatMaN file processor of Omnitool
    ####################################################

    Hi there! Please enter any number of <patman files>. By default
    I create stacked fasta files from PatMaN files. Multiple hits
    will result only a single read. Alternatively I can make a 
    matchtable but than I will need a fasta input too.

    It can be specified to output only sense (+) or antisense (-), the
    default is to output both; accaptable mismatch score of hits.
    A total match count table can be created in a separate file.

    -c               Suppress patman to fasta conversion

    -o   +|-         Define orientation (default is both).

    -s               Suppress the header for tables.

    -e   NUM         Define mismatches to accept (default is 0).

    -l   STRING      A label that is added to the tables.

    -i   FILE        Input a fasta file. Some option needs the
                     fasta file containanig the seqs to look for.
                     (Usually the database file of the patman.
                     Alternatively, any fasta file can be used to
                     filter the patman file(s).)
    
    -n   FILE        Input a file with the normalization factors 
                     to directly normalize the outputs.
                     (normalizer can create this for you.)


    -O   FILE        Specify a name of the fasta output file 

    -T   FILE        Create a file with the total readcount.

    -P   DIR         Create a table file and bar charts of 
                     positive and negative readcount ratio in
                     this directory.

    -H   FILE        Create a file with hit counts. (requires -i)

    -C   DIR         Create coverage plots (requires -i)

    -y   NUM         Fix the Y axis for the coverage plots

    -a   DIR         Output into a generated file(s) in this
                     directory (generates names just like when
                     multiple input is given)


    Version 0.1 (alpha)
    I was written by Levente Kontra

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


__EOUSAGE__
    ;    
my %options;
getopts("a:co:O:sel:T:P:H:C:i:y:n:",\%options);

my @input = @ARGV;


#safe guard
if (!@input){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Missing input\n";
    die $usage;
}
if (($#ARGV)&&($options{O})){
    print STDERR "\nERROR!!\n$0:\n";
    print STDERR "Cannot handel -o and multiple input at the same time.\n";
    die $usage;
}

if (($options{H})&&(!$options{i})){
    print STDERR "ERROR!!\n$0:";
    print STDERR "To create a hit count table -i option (fasta input) is required\n";
    die $usage;
}

if (($options{C})&&(!$options{i})){
    print STDERR "ERROR!!\n$0:";
    print STDERR "To create a coverage plot both a fasta input (-i) is required\n";
    die $usage;
}
if (($options{o})&&($options{P})){
    print STDERR "ERROR!!\n$0:";
    print STDERR "There is no point in counting positive-negative ratio if considering only one orientation\n";
    die $usage;
}
if (($options{O})&&($options{c})){
    print STDERR "ERROR!!\n$0:";
    print STDERR "There is no fasta to output if the converstion is canceled. I do not know what to do with options -O and -
c at the same time\\n";
    die $usage;
}
if (($options{a})&&($options{c})){
    print STDERR "ERROR!!\n$0:";
    print STDERR "There is no fasta to output if the converstion is canceled. I do not know what to do with options -a and -c at the same time\n";
    die $usage;
}
$SIG{'INT'} = sub {die "interrupted, I am aborting...\n"};


#defaults
unless ($options{e}){$options{e} = 0;}
unless ($options{l}){$options{l} ="";}
unless ($options{a}){$options{a} = "./fasta/";}

#creating direcoties
if ($options{P}){
    unless(-e "./pos-neg/"){
	system "mkdir pos-neg";
    }
}
if ($options{C}){
    unless(-e "./$options{C}/"){
	system "mkdir $options{C}";
    }
}
if ($#ARGV){
    unless($options{a}){
	unless(-e "./fasta/"){
	    system "mkdir fasta";
	}
    }
}
my $header;
my $total;
my $positive_ratio;
my $negative_ratio;
my %hits_table;
my %seq_lengths;

#reading in fasta
if ($options{i}){
    my $seqio_object = Bio::SeqIO->new(-file => $options{i}, '-format' => 'Fasta');
    while (my $seq = $seqio_object->next_seq){
	my $id = $seq->id();
	my $length = $seq->length;
	if ($options{H}){
	    $hits_table{$id} .= "$id";
	}
	if ($options{C}){
	    $seq_lengths{$id}=$length;
	}
    }
}    

#reading in normalization factors
my $norm_string;
my @normalization;
if ($options{n}){
    open (NORM, $options{n}) or die "Can't open normalization table input file $options{n}\n";
    while (<NORM>){
        if ($_=~/^normalization factor:\t(.*)/){
	    $norm_string = $1;
	}
    }
@normalization = split (/\t/, $norm_string);
}

#reading in patman input(s)
my $i=0;

my @head = @input;
foreach my $head (@head){
    $head =~ s/.*\///;
    $head =~ s/\.\w+$//;
    $header .= "\t$head";
}

foreach my $input (@input){
    my %plus1;
    my %minus1;
    my %hits;
    my %seen;
    my $matching;
    my $positive;
    my $negative;

    foreach my $element (keys %seq_lengths){
	my $length = $seq_lengths{$element};
	my $i = 0;
	while ($i < $length){
	    ${$plus1{$element}}[$i] = "0";
	    ${$minus1{$element}}[$i] = "0";
	    $i++;
	}
    }

    open (INPUT, "$input") or die "Couldn't open patman input file\n";
    my $namer = $input;
    $namer =~ s/.*\///;
    $namer =~ s/\.\w+$//;
    while (<INPUT>){
	if ($_=~m/^(\S+)\s*.*\t(\w+)\((\S+)\)\t(\d+)\t(\d+)\t(\S+)\t(\d+)/){
	    my $name = $1;
	    my $sequence = $2;
	    my $abundance = $3;
	    my $start = $4;
	    my $end = $5;
	    my $orientation = $6;
	    my $mismatch = $7;
	    if ($options{n}){
		$abundance *= $normalization[$i];
	    }
	    $sequence =~ tr/U/T/;
	    $sequence =~ tr/a-z/A-Z/;
	    chomp ($name, $sequence, $abundance, $start, $end, $orientation, $mismatch);
	    if ($mismatch != $options{e}){next;}
	    if (($options{o}) && ($orientation ne $options{o})){next;}

	    if ($options{H}){$hits{$name} += $abundance;}
	    unless ($seen{$sequence}){
		$matching += $abundance;
		$seen{$sequence} = $abundance;
		if ($options{P}){
		    if ($orientation eq "+"){$positive += $abundance;}
			    else {$negative += $abundance}
		}
	    }
	    if ($options{C}){
		my $i = $start -1;
		while ($i < $end-1){
		    if ($orientation eq '+'){
			$plus1{$name}[$i] += $abundance; 
		    }
		    else{
			$minus1{$name}[$i] += $abundance; 
		    }
		    $i++;
		}
	    }
	}
    }	
    if (!$matching){$matching = 0;}
    
    if ($options{T}){
    $total .= "\t$matching";
    }

    if ($options{P}){
	if (!$positive){$positive =0;}
	else{$positive = (($positive/$matching)*100);}
	$positive =sprintf("%.2f", $positive);
	$positive_ratio .= "\t$positive";
	if (!$negative){$negative =0;}
	else{$negative = (($negative/$matching)*100);}
	$negative =sprintf("%.2f", $negative);
	$negative_ratio .= "\t$negative";

	open (POSITIVE, ">", "plot.tmp") or die "Can't create temp positive ratio file\n";
	print POSITIVE "file name:\t$namer\n";
	print POSITIVE "$options{l} positive read ratio:\t$positive\n";
	print POSITIVE "$options{l} negative read ratio:\t$negative\n";
	close POSITIVE;
	
	open (PLOT, ">patmanator.plot") or die "Can't create temp plot file\n";
	print PLOT "data<-read.table(\"plot.tmp\",sep=\"\\t\", header = TRUE, row.names = 1)\n";
	print PLOT "pdf(\"$options{P}/$namer.pos-neg.pdf\", width = 2, height = 5)\n";
	print PLOT "par(mar=c(3.1, 2.1, 2.1, 6.1), xpd=TRUE)\n";
	print PLOT "barplot(apply(as.matrix(data),2,rev), ylim=c(0,100), col=c(\"blue\",\"red\"))\n";
	print PLOT "legend(\"topright\", inset=c(-3.3,0), c(\"positive\",\"negative\"), fill=c(\"red\",\"blue\"))\n";
	print PLOT "dev.off()\n";
	print PLOT "quit(\"no\")\n";
	system "R --vanilla --quiet --slave < patmanator.plot >blabla.log";
	unlink "plot.tmp";
	unlink "blabla.log";
	unlink "patmanator.plot";
    }

    if ($options{H}){
	foreach my $element (keys %hits_table){
	    unless ($hits{$element}){
		$hits_table{$element} .= "\t0";
	    }
	    else{
		    $hits_table{$element} .= "\t$hits{$element}";
	    }
	}    
	
    }    
    if ($options{C}){
	foreach my $locus (keys %plus1){
	    my $id = $locus;
	    chomp $id;
	    my $min = 0;
	    my $max = 0;
	    open (COV, ">cov.tmp") or die "Can't create tmp output file\n";
	    my $i = 1;
	    foreach my $element (@{$plus1{$id}}){
		my $p1 = $element;
		my $m1 = $minus1{$id}[$i-1];
		unless ($options{y}){
		    if ($p1 > $max){
			$max = $p1;
		    }
		    if ($m1 > $min){
			$min = $m1;
		    }
		}
		if ($m1 > 0){
		    $m1 = "\-$m1";
		}
		print COV "$i\t$p1\t$m1\n";
		$i++;
	    }
	    if ($options{y}){
		$min = $options{y};
		$max = $options{y};
	    }
	    close COV;
	    open (COVP, ">cov.plot") or die "Can't create temp plot file\n";
	    print COVP "data<-read.table(\"cov.tmp\",sep=\"\\t\")\n";
	    print COVP "pdf(\"./$options{C}/$namer-coverage.pdf\")\n";
	    print COVP "r.y<-c(data[,1],data[,2])\n";
	    print COVP "min<-min(-$min)\n";
	    print COVP "max<-max($max)\n";
	    print COVP "plot(data[,1],data[,2], ylim = c(min,max), type=\"h\", xlab=\"\", ylab=\"\", col=\"red\", lwd=\"0.1\", bg=\"4\", main=\"$id\")\n";
	    print COVP "par(new=TRUE)\n";
	    print COVP "plot(data[,1],data[,3], ylim = c(min,max), type=\"h\", col=\"4\", lwd=\"0.1\", bg=\"4\", xlab=\"Position\", ylab=\"Number of reads\")\n";
	    print COVP "dev.off()\n";
	    print COVP "quit(\"no\")\n";
	    system "R --vanilla --slave < cov.plot >blabla.log";
	    unlink "cov.plot";
	    unlink "blabla.log";
	    unlink "cov.tmp";
	}
    }



    unless ($options{c}){
	my @keys = sort { $seen{$b} <=> $seen{$a} } keys %seen;
	if (($#ARGV)||($options{a})){
	    $input=~ s/.patman$/.fa/ or die "Incorrect file ending:$input\n";
	    $input =~ s/^.*\///;
	    open (MULTI_OUT, ">", "$options{a}/$input") or die "Can't create output file: $input\n";
	    foreach my $element (@keys){
		print MULTI_OUT ">$element\($seen{$element}\)\n$element\n";
	    }
	    close MULTI_OUT;
	}
	elsif (defined $options{O}){
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

#outputing tables

if ($options{T}){
    open (SUM, ">", "$options{T}") or die "Can't create output summary file: $options{T}\n";
    unless ($options{s}){
	print SUM "file name:$header\n";
    }
    print SUM "$options{l} matching:$total\n";
    close SUM;
}

if ($options{P}){
    open (POSITIVE, ">", "$options{P}/pos-neg.table.txt") or die "Can't create positive-negative ratio file: $options{P}\n";
    unless ($options{s}){
	print POSITIVE "file name:$header\n";
    }
    print POSITIVE "$options{l} positive read ratio:$positive_ratio\n";
    print POSITIVE "$options{l} negative read ratio:$negative_ratio\n";
    close POSITIVE;
    unlink "blabla.log";
}

if ($options{H}){
    open (COUNTS, ">", "$options{H}") or die "Can't create output table file: $options{m}\n";
    unless ($options{s}){
	print COUNTS "file name:$header\n";
    }
    foreach my $element (keys %hits_table){
	print COUNTS "$hits_table{$element}\n";
    }
    close COUNTS;
}
