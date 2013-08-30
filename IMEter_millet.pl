#!usr/bin/perl
use strict; use warnings;
use FAlite;
use Getopt::Std;
use vars qw($opt_g);
our $opt_g; # alternative to our()
getopts('g');
die "Usage: [option] <FASTA> <DATA>" unless @ARGV == 2;


my ($fasta_file, $data_file) = @ARGV; 
my ($k) = $data_file =~ /(\d+)/;
my %data;	
	
open(DATA, $data_file) or die "error opening $data_file for processing\n";
while (<DATA>) {
	next if /^#/;
	my ($kmer, $freq) = split;
	$data{$kmer} = $freq;
}

open(my $in, "<$ARGV[0]") or die "error reading $ARGV[0] for reading";
my $fasta = new FAlite($in);
while (my $entry = $fasta->nextEntry) { #extracts header and sequence from file
	my $header = $entry->def;
	my $seq = $entry->seq;
	my ($dist) = $header =~ /(\d+)$/;
	my $sum = IMEsum($seq, %data);
	if ($opt_g) {print "$dist\t$sum\n";}
	else	    {print "$header\n $sum\n";}
}

close $in;


	

sub IMEsum {
	my ($seq, %data) = @_;
	
	my $sum = 0;
	my $k = 5;
	for (my $i = 0; $i < (length($seq) - ($k-1)); $i ++) {
		my $K_mer = substr($seq, $i, $k);
		if (exists $data{$K_mer}) {
			$sum += $data{$K_mer}
			
		}
	}
	return $sum;
}