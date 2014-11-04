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
use vars qw($opt_h $opt_p $opt_d $opt_P $opt_D $opt_k $opt_c $opt_i $opt_1 $opt_3 $opt_v);
use FAlite;

my $cmdline = "$0 @ARGV";

getopts('hp:d:P:D:k:ci13v');

die "
usage: ime_trainer.pl [params] <fasta file(s) with meta-data>
params:
  -p <int> proximal cutoff using nt coordinates
  -d <int> distal cutoff using nt coordinates
  -P <int> proximal cutoff using intron position
  -D <int> distal cutoff using intron position
  -k <int> word size
  -c       require complete genes
  -i       require primary isoform
  -1       IMEter 1 style training
  -3       IMEter 3 style training (decay)
  -v       verify only, do not produce output file
" if @ARGV == 0 or $opt_h;

my @file = @ARGV;

die "-k is required\n" unless $opt_k;
my $cutoffs_ok;
if    ($opt_p and $opt_d) {$cutoffs_ok = 1}
elsif ($opt_P and $opt_D) {$cutoffs_ok = 1}
else                      {die "pairs of -p and -d or -P and -D required\n"}

die "IMEter 3-style training not yet supported\n" if $opt_3;

# main loop

my %count;
my %log;
my %transcript_count;

foreach my $file (@ARGV) {
	open(my $fh, $file) or die;
	my $fasta = new FAlite($fh);
	while (my $entry = $fasta->nextEntry) {
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
		
		# fix to skip CG rich regions in first 250 bp
#		if ($beg <= 300){
#			$class = 'not_counted';
#		}
		
		if ($opt_p and $beg <= $opt_p) {
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
		$log{$class}++;
		
		# count intron
		$log{intron}++;
		for (my $i = 0; $i < length($seq) - $opt_k + 1; $i++) {
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

print "# transcripts:               ", scalar keys %transcript_count, "\n";
print "# transcripts not counted:   ", scalar keys %{$log{skipped_incomplete}}, "\n";

print "# introns counted:     ", $log{intron}, "\n";
print "# introns not counted: ", scalar keys %{$log{skipped_secondary}}, "\n";

print "# proximal introns:    ", $log{proximal}, "\n";
print "# distal introns:      ", $log{distal}, "\n";

print "# skipped kmers:       ";
my ($kerr) = sort {$log{skipped_kmer}{$b} <=> $log{skipped_kmer}{$a}} keys %{$log{skipped_kmer}};
if ($kerr) {
	print scalar keys %{$log{skipped_kmer}}, " maximum = ",
		"$kerr ($log{skipped_kmer}{$kerr})\n";
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
