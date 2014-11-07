#!/usr/bin/perl
use strict; use warnings;
use FAlite;
use DataBrowser;

# how far away from TSS do we want to go?
# this is just for reporting purposes, will still count all kmers past this point
# to establish background dinucleotide frequencies for all introns
my $LIMIT = 1000;

die "Usage: $0 <input intron FASTA file> <window size for smoothing>" unless @ARGV == 2;
my ($file, $window) = @ARGV;

# counts of individual dinucleotides (irrespective of position)
my %c0;

# count of all dinucleotides (irrespective of position)
my $c0_total = 0;

# set up blank data structure with pseudocounts
# counting dinucleotides w.r.t. position from TSS
my @c1;
my @alph = qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT);
for (my $i = 0; $i < $LIMIT; $i++) {
	foreach my $nt (@alph) {
		$c1[$i]{$nt} = 1;
	}
}

# collect raw data
open(my $in, "<", $file) or die "Can't open $file";
my $fasta = new FAlite($in);
while (my $entry = $fasta->nextEntry) {

	# extract start position of intron (relative to TSS) to use as offset
	my ($intron_start_pos) = $entry->def =~ m/COORDS=(\d+)-/;
	my $seq = $entry->seq;
	# extract dinucleotides at each position
	DINUC: for (my $i = 0; $i < length($seq) -1; $i++) {
		my $kmer = substr($seq, $i, 2);

		# skip any dinucleotides with ambiguity codes
		next DINUC unless ($kmer =~ m/[ACGT]{2}/);

		# calculate position w.r.t. TSS but -1 to work with zero-index in array
		my $pos = $i + $intron_start_pos - 1;
		$pos = $LIMIT if $pos > $LIMIT;
		$c1[$pos]{$kmer}++;
		$c0{$kmer}++;
		$c0_total++;
	}
}
close $in;

# smooth data, place into new data structure
# take $window bases either side of current position
my @c2;

for (my $i = $window; $i < @c1 - $window; $i++) {
	for (my $j = -$window; $j <= $window; $j++) {
		foreach my $kmer (keys %{$c1[$i+$j]}) {
			$c2[$i]{$kmer} += $c1[$i+$j]{$kmer}; 
		}
	}
}

# not sure what point of these two lines are
# making just 1st and last elements of smoothed array equal to unsmoothed array?
#$c2[0] = $c1[0];
$c2[@c1 -1] = $c1[@c1 -1];


# output frequencies
for (my $i = 0; $i < @c2; $i++) {
	print $i + 1, "\t";

	# calculate total of …?
	my $total = 0;
	foreach my $kmer (keys %{$c2[$i]}) {$total += $c2[$i]{$kmer}}

	# loop over all dinucleotides
	foreach my $kmer (sort keys %{$c2[$i]}) {
		print "$kmer ";
		my $raw_dinuc_freq = ($c0{$kmer} / $c0_total);
		my $smoothed_dinuc_freq = $c2[$i]{$kmer} / $total;
#		printf "%.3f", $smoothed_dinuc_freq;

		# calculate ratio of smoothed dinucleotide frequency at current position from TSS
		# to raw dinucleotide frequency across all positions
		my $ratio = sprintf("%.2f", $smoothed_dinuc_freq / $raw_dinuc_freq);
#		printf "%.2f", $ratio;

		if($ratio >= 1.25 or $ratio <= 0.75){
			printf "%.2f", $ratio;
		} else{
			print "   ";
		}
		
		# flag ratios that are either doubled or halved with respect to background.
		if($ratio >= 2){
			print "+ ";
		} elsif ($ratio <= 0.5){
			print "- ";
		} else {
			print " ";
		}
	}
	print "\n";
}

