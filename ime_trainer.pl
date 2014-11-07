#!/usr/bin/perl
# Copyright (C) 2007-2014 Ian Korf and Keith Brandam
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

use strict;
use warnings 'FATAL' => 'all';
use Getopt::Std;
use vars qw($opt_h $opt_p $opt_q $opt_d $opt_P $opt_D $opt_k $opt_c $opt_i $opt_1 $opt_v $opt_x);
use FAlite;

my $cmdline = "$0 @ARGV";

getopts('hp:d:q:P:D:k:ci1vx');

die "
usage: ime_trainer.pl [params] <fasta file(s) with meta-data>
params:
  -p <int> proximal start cutoff using nt coordinates
  -q <int> proximal end cutoff using nt coordinates
  -d <int> distal cutoff using nt coordinates
  -P <int> proximal cutoff using intron position
  -D <int> distal cutoff using intron position
  -k <int> word size
  -c       require complete genes
  -i       require primary isoform
  -1       IMEter 1 style training
  -x       Exclude regions of introns that span -q and -p cutoff values
           i.e. if using -q, can use some sequence of introns which start before -q value
           but end after it. Likewise, for -p, can trim introns that start before value
           of -p but end after it.
  -v       verify only, do not produce output file
" if @ARGV == 0 or $opt_h;

my @file = @ARGV;

# check command-line options
die "-k is required\n" unless $opt_k;
my $cutoffs_ok;
if (($opt_p or $opt_d) and ($opt_P or $opt_D)) {
	die "Use either -p and -d or -P and -D but do not mix the two\n";
} elsif ($opt_p and not $opt_d){
	die "-p option must be paired with -d option (and optionally -q)\n";
} elsif ($opt_P and not $opt_D){
	die "-P option must be paired with -D option\n";
} elsif (not $opt_p and not $opt_P){
	die "Must specify -p and -d (and optionally -q), or -P and -D options\n";
} elsif ($opt_q and $opt_P){
	die "-q option should only be used with -p option\n";
} elsif ($opt_x and $opt_P){
	die "-x option should only be used with -p/-q options\n";
} else {warn "Options look good\n"}


# main loop

my %count;

# will have various log statistics, set a couple of them to zero
my %log;
$log{not_counted} = 0;

my %transcript_count;

foreach my $file (@ARGV) {
	open(my $fh, $file) or die;
	my $fasta = new FAlite($fh);
	SEQ: while (my $entry = $fasta->nextEntry) {
		my $seq = $entry->seq;
	
		# parse fasta header
		my ($uid, $type, $pos, $coor, $id1, $id2, $iso, $strct) = split(/\s/, $entry->def);
		($type) = $type =~ /TYPE=(\S+)/;
		($pos) = $pos =~ /POS=(\d+)/;
		my ($beg, $end) = $coor =~ /COORDS=(\d+)\-(\d+)/;
		($iso) = $iso =~ /ISOFORM=(\S+)/;
		($strct) = $strct =~ /STRUCTURE=(\S+)/;
		$transcript_count{$id1}++;
		
		# filtering
		if ($opt_c and $strct ne '5-3') {
			$log{skipped_incomplete}{$id1}++;
			next;
		}
		if ($opt_i and $iso ne 'primary') {
			$log{skipped_secondary}{$id1}++;
			next;
		}

		# classify intron
		my $class;
		
		# -q option to skip CG-rich regions in first part of transcript
		# and other dinucleotide biases, will handle this slightly differently
		# if using -x option later
		if ($opt_q and $beg <= $opt_q and not $opt_x){
			$class = 'not_counted';
		} elsif ($opt_q and $opt_x and ($end <= $opt_q)){
			# we can exclude some introns when using -q and -x if they end before
			# the value of -q, otherwise we will be counting at least part of the intron
			$class = 'not_counted';
		} elsif ($opt_p and $beg <= $opt_p) {
			$class = 'proximal';
		} elsif ($opt_P and $pos <= $opt_P) {
			$class = 'proximal';
		} elsif ($opt_d and $beg >= $opt_d) {
			$class = 'distal';
		} elsif ($opt_D and $pos >= $opt_D) {
			$class = 'distal'
		} else {
			$class = 'not_counted';
		}

		# keep track of how many introns we have seen in total
		$log{intron}++;

		# keep track of how many introns in each class
		$log{$class}++;
		
		# now want to count kmers in introns, but no real need to proceed if
		# $class isn't proximal or distal at this point (this means we don't
		# count non-ACGT kmers from uncounted introns)
		next SEQ unless ($class eq 'proximal' or $class eq 'distal');
		
		# now count kmers in our proximal or distal introns
		KMER: for (my $i = 0; $i < length($seq) - $opt_k + 1; $i++) {

			# if we are using -x option with -q, then we only want to look
			# at kmers that start *after* we pass the value of -q
			next KMER if ($class eq 'proximal' and $opt_q and $opt_x and ($i < $opt_q));


			# if we are using -x option with -p, then we only want to look
			# at kmers in intron up to the distance specified by -p. Can just ignore
			# any sequence after that point
			my $distance_from_TSS = $i + $beg;
			next SEQ if ($class eq 'proximal' and $opt_x and ($distance_from_TSS > $opt_p));

			
			# can now grab kmers and start counting
			my $kmer = substr($seq, $i, $opt_k);

			if ($kmer =~ /^[ACGT]+$/) {
 				$count{$kmer}{$class}++;
			} else {
				$log{skipped_kmer}{$kmer}++;
			}
		}
		
	}
}

# status report
no warnings;

# output section
print "# IME training log\n";
print "# ----------------\n";
print "#\n";
print "# ", `date`;
print "# build: $cmdline\n";
print "#\n";
print "# transcripts:         ", scalar keys %transcript_count, "\n";
print "# incomplete:          ", scalar keys %{$log{skipped_incomplete}}, "\n";
print "# secondary isoforms:  ", scalar keys %{$log{skipped_secondary}}, "\n";
print "#\n";
print "# introns counted:     ", $log{intron}, "\n";
print "# introns exluded:     ", $log{not_counted}, "\n";
print "#\n";
print "# proximal introns:    ", $log{proximal}, "\n";
print "# distal introns:      ", $log{distal}, "\n";
print "#\n";
print "# skipped kmers:       ";
my ($kerr) = sort {$log{skipped_kmer}{$b} <=> $log{skipped_kmer}{$a}} keys %{$log{skipped_kmer}};
if ($kerr) {
	print scalar keys %{$log{skipped_kmer}}, " maximum = ",
		"$kerr ($log{skipped_kmer}{$kerr})\n";
} else{
	print "0\n";
}
my $low_counts;
foreach my $kmer (keys %count) {
	
	# pseudocount everything to make sure we have no zero
	# counts later on (can't take a log of zero)
	$count{$kmer}{proximal}++;
	$count{$kmer}{distal}++;

	# keep track of those which are still 'low'
	if ($count{$kmer}{proximal} < 11 or 
		$count{$kmer}{distal}   < 11) {
		$low_counts++;
	}

}
if ($low_counts) {
	print "# WARNING: low k-mer counts in $low_counts k-mers, decrease -k\n";
}
exit if $opt_v;

my $ptotal;
my $dtotal;
foreach my $kmer (keys %count) {
	$ptotal += $count{$kmer}{proximal};
	$dtotal += $count{$kmer}{distal};
}
my %lod;
foreach my $kmer (keys %count) {
	$lod{$kmer} = log(($count{$kmer}{proximal}/$ptotal) / ($count{$kmer}{distal}/$dtotal)) / log(2), "\n";
}

print "#\n";
print "# IME parameters\n";
print "# --------------\n";


foreach my $kmer (sort {$lod{$b} <=> $lod{$a}} keys %lod) {
	print "$kmer\t$lod{$kmer}\n";
}
