#!/usr/bin/perl
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
my %genecount;
my %isocount;

foreach my $file (@ARGV) {
	open(my $fh, $file) or die;
	my $fasta = new FAlite($fh);
	while (my $entry = $fasta->nextEntry) {
		my $seq = $entry->seq;
	
		# parse fasta header
		my ($uid, $type, $pos, $coor, $id1, $id2, $iso, $strct)
			= split(/\s/, $entry->def);
		($type) = $type =~ /TYPE=(\S+)/;
		($pos) = $pos =~ /POS=(\d+)/;
		my ($beg, $end) = $coor =~ /COORDS=(\d+)\-(\d+)/;
		($iso) = $iso =~ /ISOFORM=(\S+)/;
		($strct) = $strct =~ /STRUCTURE=(\S+)/;
		$genecount{$id1}++;
		$isocount{$id2}++;
		
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
print STDERR "genes:               ", scalar keys %genecount, "\n";
print STDERR "isoforms             ", scalar keys %isocount, "\n";
print STDERR "genes not counted:   ", scalar keys %{$log{skipped_incomplete}}, "\n";

print STDERR "introns counted:     ", $log{intron}, "\n";
print STDERR "introns not counted: ", scalar keys %{$log{skipped_secondary}}, "\n";

print STDERR "proximal introns:    ", $log{proximal}, "\n";
print STDERR "distal introns:      ", $log{distal}, "\n";

print STDERR "skipped kmers:       ";
my ($kerr) = sort {$log{skipped_kmer}{$b} <=> $log{skipped_kmer}{$a}} keys %{$log{skipped_kmer}};
if ($kerr) {
	print STDERR scalar keys %{$log{skipped_kmer}}, " maximum = ",
		"$kerr ($log{skipped_kmer}{$kerr})\n";
}
my $low_counts;
foreach my $kmer (keys %count) {
	if (not defined $count{$kmer}{proximal} or
		not defined $count{$kmer}{proximal} or
		$count{$kmer}{proximal} < 10 or 
		$count{$kmer}{distal} < 10) {
		$low_counts++;
	}
}
if ($low_counts) {
	print STDERR "WARNING: low k-mer counts in $low_counts k-mers, increase -k\n";
}
exit if $opt_v;

# output section
print "# ", `date`;
print "# build: $cmdline\n";
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
foreach my $kmer (sort {$lod{$b} <=> $lod{$a}} keys %lod) {
	print "$kmer\t$lod{$kmer}\n";
}
