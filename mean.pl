#!/usr/bin/perl
#
# contig_stats.pl
#
# a quick script to calculate the mean and N50 sizes for a file of
# contig or scaffold sequences
#
# by Keith Bradnam
#
# Last updated by: $Author$
# Last updated on: $Date$
##############################################

use strict;
use warnings;
use FAlite;

my @files = glob("*intron.fa");


foreach my $file (@files){
	warn "Processing $file\n";
	my $species = $file;
	$species =~ s/_IME_intron.fa//;

	my ($n, $running_length) = (0, 0);
	# parse info from FASTA file
	open(my $in,"$file") or die "Can't read from $file\n";
	my $fasta = new FAlite($in);

	# loop through each sequence in target file
	while(my $entry = $fasta->nextEntry){
		my $seq = $entry->seq;
		$n++;
		$running_length += length($seq);
	}	
	
	close($in);
	my $average_length = sprintf("%.1f", $running_length / $n * 100);
	print "$species\t$n\t$average_length\n";
}
exit(0);