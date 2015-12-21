#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

mkdir "./plots";
my $input1 = shift;
my $input2 = shift;
my $input3 = shift;
my $fasta = shift;
my $scale = shift;
if ((!$input1)||(!$input2)|| (!$input3)||(!$fasta)){
	print STDERR "$0: a script to plot sRNA reads";
	print STDERR "$0 please enter <patman input file 1> <patman input file 2> <patman input file 3> and a stacked <fasta file> that was added into PatMaN as database. Custum scalebar langth can be specified, but it is optional <scale>\n";
	exit();
	}
#                             green                    blue                red
my %seq_lengths;	
my $seqio_object = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
while (my $seq = $seqio_object->next_seq){
	my $id = $seq->id();
	my $length = $seq->length;
	$seq_lengths{$id}=$length;
}
my %plus1;
my %minus1;
my %plus2;
my %minus2;
my %plus3;
my %minus3;
foreach my $element (keys %seq_lengths){
	my $length = $seq_lengths{$element};
	my $i = 0;
	while ($i < $length){
		push (@{$plus1{$element}}, "0");
		push (@{$minus1{$element}}, "0");
		push (@{$plus2{$element}}, "0");
		push (@{$minus2{$element}}, "0");
		push (@{$plus3{$element}}, "0");
		push (@{$minus3{$element}}, "0");	
		$i++;
	}
}

open (INPUT1, "$input1") or die "Couldn't open patman input file\n";
while (<INPUT1>){
	if ($_=~m/^(\S+)\s*.*\t(\w+)\((\S+)\)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/){
		my $name = $1;
		my $sequence = $2;
		my $abundance = $3;
		my $start = $4;
		my $end = $5;
		my $strand = $6;
		my $mm = $7;
		if ($mm == 0){
			my $i = $start -1;
			while ($i < $end-1){
				if ($strand eq '+'){
					$plus1{$name}[$i]+=$abundance; 
				}
				else{
					$minus1{$name}[$i]+=$abundance; 
				}
			$i++;
			}
		}	
	}	
}
close INPUT1;

open (INPUT2, "$input2") or die "Couldn't open patman input file\n";
while (<INPUT2>){
	if ($_=~m/^(\S+)\s*.*\t(\w+)\((\S+)\)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/){
		my $name = $1;
		my $sequence = $2;
		my $abundance = $3;
		my $start = $4;
		my $end = $5;
		my $strand = $6;
		my $mm = $7;
		if ($mm == 0){
			my $i = $start -1;
			while ($i < $end-1){
				if ($strand eq '+'){
					$plus2{$name}[$i]+=$abundance; 
				}
				else{
					$minus2{$name}[$i]+=$abundance; 
				}
			$i++;
			}
		}
	}	
}
close INPUT2;

open (INPUT3, "$input3") or die "Couldn't open patman input file\n";
while (<INPUT3>){
	if ($_=~m/^(\S+)\s*.*\t(\w+)\((\S+)\)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/){
		my $name = $1;
		my $sequence = $2;
		my $abundance = $3;
		my $start = $4;
		my $end = $5;
		my $strand = $6;
		my $mm = $7;
		if ($mm == 0){
			my $i = $start -1;
			while ($i < $end-1){
				if ($strand eq '+'){
					$plus3{$name}[$i]+=$abundance; 
				}
				else{
					$minus3{$name}[$i]+=$abundance; 
				}
				$i++;
			}
		}
	}	
}
close INPUT3;

foreach my $locus (keys %plus1){
	my $id = $locus;
	chomp $id;
	my $min = 0;
	my $max = 0;
	open (TMP, ">tmp.$$") or die "Can't create tmp output file\n";
	my $i = 1;
	foreach my $element (@{$plus1{$id}}){
		my $p1 = $element;
		my $m1 = $minus1{$id}[$i-1];
		my $p2 = $plus2{$id}[$i-1];
		my $m2 = $minus2{$id}[$i-1];
		my $p3 = $plus3{$id}[$i-1];
		my $m3 = $minus3{$id}[$i-1];		
			if (!$scale){
				if ($p1 > $max){
					$max = $p1;
				}
				if ($p2 > $max){
					$max = $p2;
				}
				if ($p3 > $max){
					$max = $p3;
				}		
				if ($m1 > $min){
					$min = $m1;
				}
				if ($m2 > $min){
					$min = $m2;
				}
				if ($m3 > $min){
					$min = $m3;
				}
			}
			else {
				$min = $scale;
				$max = $scale;
			}
		if ($m1 > 0){
			$m1 = "\-$m1";
		}
		if ($m2 > 0){
			$m2 = "\-$m2";
		}
		if ($m3 > 0){
			$m3 = "\-$m3";
		}		
		print TMP "$i\t$p1\t$m1\t$p2\t$m2\t$p3\t$m3\n";
		$i++;
	}
close TMP;
	
&plot($id, $max, $min);
unlink "tmp.$$";
}

sub plot{
my $id = shift;
my $max = shift;
my $min = shift;
	if ($id=~/(\S+)\/(\d+)\-(\d+)\D+/){
	my $chr = $1;
	my $start =$2;
	my $end = $3;
	chomp ($chr, $start, $end);
	$id ="$chr"."_"."$start"."_"."$end";
	}
open (NEW, ">plot.$$") or die "Can't create temp file plot.$$\n";
print NEW "data<-read.table(\"tmp.$$\",sep=\"\\t\")\n";
print NEW "pdf(\"$id.pdf\")\n";
print NEW "r.y<-c(data[,1],data[,2])\n";
print NEW "min<-min(-$min)\n";
print NEW "max<-max($max)\n";
print NEW "plot(data[,1],data[,2], ylim = c(min,max), type=\"l\", xlab=\"\", ylab=\"\", col=\"3\", lwd=\"1\", bg=\"4\", main=\"$id\")\n";
print NEW "lines(data[,1],data[,4], col=\"4\", lwd=\"1\", bg=\"3\")\n";
print NEW "lines(data[,1],data[,6], col=\"2\", lwd=\"1\", bg=\"2\")\n";
print NEW "par(new=TRUE)\n";
print NEW "plot(data[,1],data[,3], ylim = c(min,max), type=\"l\", col=\"3\", lwd=\"1\", bg=\"4\", xlab=\"Position\", ylab=\"Number of reads\")\n";
print NEW "lines(data[,1],data[,5], col=\"4\", lwd=\"1\", bg=\"3\")\n";
print NEW "lines(data[,1],data[,7], col=\"2\", lwd=\"1\", bg=\"2\")\n";
print NEW "dev.off()\n";
print NEW "quit(\"no\")\n";
system "R --vanilla --quiet < plot.$$";
unlink "plot.$$";
system "mv $id.pdf ./plots";
}

# made by Simon Moxon mostly, a bit moded by Levente Kontra
