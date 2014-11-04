#!/usr/bin/perl
use strict; use warnings;
use FAlite;
use DataBrowser;

my $LIMIT = 2000;

my %c0;
my $c0_total = 0;

# set up blank data structure with pseudocounts
my @c1;
my @alph = qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT);
for (my $i = 0; $i < $LIMIT; $i++) {
	foreach my $nt (@alph) {
		$c1[$i]{$nt} = 1;
	}
}

# collect raw data
my $file = shift @ARGV;
open(IN, $file) or die;
my $fasta = new FAlite(\*IN);
while (my $entry = $fasta->nextEntry) {
	my ($offset) = $entry->def =~ m/COORDS=(\d+)-/;
	my $seq = $entry->seq;
	for (my $i = 0; $i < length($seq) -1; $i++) {
		my $kmer = substr($seq, $i, 2);
		my $pos = $i + $offset - 1;
		$pos = $LIMIT if $pos > $LIMIT;
		$c1[$pos]{$kmer}++;
		$c0{$kmer}++;
		$c0_total++;
	}
}
close IN;

# smooth data
my @c2;
my $w = 9;
for (my $i = $w; $i < @c1 - $w; $i++) {
	for (my $j = -$w; $j <= $w; $j++) {
		foreach my $kmer (keys %{$c1[$i+$j]}) {
			$c2[$i]{$kmer} += $c1[$i+$j]{$kmer};
		}
	}
}
$c2[0] = $c1[0];
$c2[@c1 -1] = $c1[@c1 -1];

# output frequencies
for (my $i = 0; $i < @c2; $i++) {
	print $i + 1, "\t";
	my $total = 0;
	foreach my $kmer (keys %{$c2[$i]}) {$total += $c2[$i]{$kmer}}
	foreach my $kmer (sort keys %{$c2[$i]}) {
		print "$kmer ";
		printf "%.3f", $c2[$i]{$kmer} / $total;
		my $ratio = sprintf("%.2f", ($c2[$i]{$kmer} / $total) / ($c0{$kmer} / $c0_total));
		if($ratio >= 4){
			print "+++ ";
		} elsif ($ratio >= 3){
			print "++ ";
		} elsif ($ratio >= 2){
			print "+ ";
		} elsif ($ratio <= 0.25){
			print "--- ";
		} elsif ($ratio <= 0.33333){
			print "-- ";
		} elsif ($ratio <= 0.5){
			print "- ";
		} else {
			print " ";
		}
	}

	my $dinuc = "GT";
	print "\t$dinuc ";
	if ($c2[$i]{$dinuc}){
		printf "%.2f ", $c2[$i]{$dinuc} / $total * 100 ;
	} else {
		print "0";
	}
	print "vs ";
	printf "%.2f ", $c0{$dinuc} / $c0_total * 100 ;

	print "\n";
}

