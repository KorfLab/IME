#!/usr/bin/perl
#
# fastify_IME_data.pl
#
# A script to process output from proc_phytozome.pl and produce FASTA formatted versions of intron sequences
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use Getopt::Std;

use vars qw($opt_h $opt_t $opt_s);
getopts('ht:s:');

die "
usage: fastify_IME_data.pl -f <type> <FASTA definition stub> <IME data file>

type must be either: i (intron), e (exon), 5 (5' UTR), or 3 (3' UTR)

The stub should be some text that will start off the FASTA definition line (e.g. 'IME_Ath').
Unique numbers will be appended to the stub (following an undersore character) and all other
meta information will follow after a space character.				
	
" if @ARGV != 1 or $opt_h or !$opt_t or !$opt_s;

my $desired_type;
$desired_type = "intron"          if ($opt_t eq 'i');
$desired_type = "exon"            if ($opt_t eq 'e');
$desired_type = "five_prime_UTR"  if ($opt_t eq '5');
$desired_type = "three_prime_UTR" if ($opt_t eq '3');

my ($file) = @ARGV;

my $stub = $opt_s;
my $n = 0;

open(my $in, "<", $file) or die "Can't open $file";

while (<$in>) {
	chomp;
	next if /^#/;
	my ($phytozome_id, $local_id, $isoform, $structure, $ortholog, $type, $beg, $end, $position, $total, $seq) = split;
	next if ($type ne $desired_type);
	$n++;
	my $tidied_seq = tidy_seq($seq);
	print ">${stub}_$n TYPE=$type POS=${position}/$total COORDS=$beg-$end ";
	print "ID1=$phytozome_id ID2=$local_id ISOFORM=$isoform STRUCTURE=$structure ORTHOLOG=$ortholog\n";
	print "$tidied_seq\n";
}
close $in;


sub tidy_seq{
    my ($seq) = @_;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;

    my ($output_seq) = "";
    my (@seq2) = "";
    my ($start,$end);

    @seq2 = split(//,$seq);
    my $length = @seq2;
    my ($to_add) = int($length/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    $output_seq .= substr($seq,$start);
    return ($output_seq);
}
