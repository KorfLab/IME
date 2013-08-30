#!/usr/bin/perl
# logodds.pl v2.0
use strict; use warnings;
use FALite;

die "Usage: logodds.pl <proximal file> <distal file> <k>" unless @ARGV == 3;

# Declare variables, hashes
my ($proxfile, $distfile, $k) = @ARGV;
my %proxfreq = kmer_freq($proxfile, $k);
my %distfreq = kmer_freq($distfile, $k);

my $expected_kmers = 4 ** $k;
if (scalar keys %proxfreq != $expected_kmers or
	scalar keys %distfreq != $expected_kmers) {
	print STDERR "wordsize may be too large, missing kmers\n";
}

foreach my $kmer (keys %proxfreq) {
	next if not defined $distfreq{$kmer};
	printf "%s\t%.4f\n", $kmer, log2($proxfreq{$kmer} / $distfreq{$kmer});
}

# a check for empty kmer counts  
my $expected_kmers = 4 ** $k;
if (scalar keys %proxfreq != $expected_kmers or
	scalar keys %distfreq != $expected_kmers) {
	print STDERR "wordsize may be too large, missing kmers\n";
}

sub log2 {
	my ($n) = @_;
	return log($n)/log(2);
}

sub kmer_freq {
	my ($file, $k) = @_;
	
	open(my $fh, $file) or die "error opening @_ for reading";
	my $fasta = new FAlite($fh);
	my %count;
	my $total;
	while (my $entry = $fasta->nextEntry) {
		my $header = $entry->def;
		my $seq = $entry->seq;
		$seq =~ tr/acgt/ACGT/;
	
		for (my $i = 0; $i < length($seq) - $k + 1; $i++) {
			my $kmer = substr($seq, $i, $k);
			$count{$kmer}++;
			$total++;
		}
	}
	close $fh;
	
	my %freq;
	foreach my $kmer (keys %count) {$freq{$kmer} = $count{$kmer} / $total}
	
	return %freq;
}
