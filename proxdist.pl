#!/usr/bin/perl
# dividing proximal and distal Sitalica introns
use strict; use warnings;
use FAlite;

die "usage: [cutoff] <fasta> \n" unless @ARGV == 2;

my ($cut) = $ARGV[0];
my ($limit) = 10000; # might be better to avoid poly-intronic genes. This is an arbitrary threshold, should be played with and validated.


open(FASTA, "<$ARGV[1]") or die;
open(PROX, ">$cut.prox") or die "error creating $cut.prox\n";
open(DIST, ">$cut.dist") or die "error creating $cut.dist\n";
my $n;
my $fasta = new FAlite(\*FASTA);
while (my $entry = $fasta->nextEntry) {
	my $hed = $entry->def;
	my $seq = $entry->seq;
	my ($start) = $entry->def =~ /(\d+)$/;
		if ($start < $cut)						    {print PROX "$hed\n$seq\n"}
		if ($start > $cut and $start < $ limit)  {print DIST "$hed\n$seq\n"}
		if ($start > $limit) 			  {$n++}
	};
	
#print "Excluded $n entries with TSS > $limit\n";	
close PROX;
close DIST;
close FASTA;