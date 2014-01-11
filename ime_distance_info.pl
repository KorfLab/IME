#!/usr/bin/perl
# Copyright (C) 2014 Ian Korf and Keith Brandam
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
use warnings FATAL => 'all';
use FAlite;
use List::Util qw(sum);
use Getopt::Std;
use vars qw($opt_h $opt_w $opt_c $opt_i);


# set up some default values
my $window_size = 200;
my $cutoff = 5000;
my $imeter = "imeter.pl";

getopts('hw:c:i:');
$window_size = $opt_w if $opt_w;
$cutoff      = $opt_c if $opt_c;
$imeter      = $opt_i if $opt_i;

die "
usage: $0.pl [options] <directory name containing IME intron sequence & parameter files>

options:
  -w       window size [default $window_size]
  -c       cutoff distance from TSS [default $cutoff]
  -i       version of IMEter program to use [default $imeter]
  -h       this help

1) Script will look in the specified directory for IME intron files ('*.IME_intron.fa')
and IME parameter files ('*.IME_intron.params'). For each intron file, it will 
run the IMEter (imeter.pl) on the specified file using the specified parameters
to produce a '*.IME_intron.ime' file containing IMEter scores for each intron.

2) Script will then take the '*.IME_intron.ime' file and for each intron, will extract
the distance from the transcription start site (TSS) from the FASTA file and append this
information to a new output file ('*.IME_intron.ime2').

3) Finally, the script will produce a third output file ('*.IME_intron.ime3') that averages
all intron scores in binned windows of distance from the TSS. The default parameter
takes windows of 200 bp from the TSS. In addition to calculating the average IMEter score
for all introns that start in that window, the script outputs 95% confidence limits for
the mean value.
  
" if @ARGV != 1 or $opt_h;


# process all IME intron files in current directory
my ($DIR) = @ARGV;

my @introns = glob("$DIR/*IME_intron.fa");

foreach my $file (@introns){
	my $param = $file;
	$param =~ s/fa$//;
	$param .= "params";
	warn "Processing file $file\n";
	
	# skip parameter files that are altogether missing
	next if (not -e $param);
	
	# skip parameter files that don't have k-mer information
	# these will be less than ~500 bytes
	if (-s $param < 500){
		warn "Error: no parameter data for $file. Skipping...\n";
		next;
	}

	# prepare output file names
	my ($out_file, $out_file1, $out_file2, $out_file3, $out_file4);
	
	$out_file = "$file";
	$out_file =~ s/fa$//;

	$out_file1 = $out_file . "ime1";
	$out_file2 = $out_file . "ime2";
	$out_file3 = $out_file . "ime3";
	$out_file4 = $out_file . "ime4";
	
	# can skip if all files exist
	next if (-e $out_file1 and -e $out_file2 and -e $out_file3 and -e $out_file4);

	# no need to make this output file if the next output already exists.
	unless (-e $out_file1){
		# will run IMEter against each combination of IME sequence and parameter file
		my $command = "$imeter -m $param $file > $out_file1";
		warn "Running: $command\n";
		system($command) && die "Can't run $command\n";					
	}

	# Grab distance from TSS information in FASTA file, store in hash
	my %tss;
	
	open(my $ffh, "$file") or die "Can't open $file\n";
	my $fasta = new FAlite($ffh);
	while (my $entry = $fasta->nextEntry) {
		my ($id, $start) = $entry->def =~ />(\S+?)\s.*COORDS=(\S+)-/;
		$tss{$id} = $start;
	}
	close $ffh;


	# need to track all IMEter scores_by_pos
	# 1) track by position
	# 2) track all scores (to work out percentile data)
	my %scores_by_pos;
	my @all_scores;
	

	open(my $in,  "<", $out_file1) or die "can't read from $out_file1\n";
	open(my $out, ">", $out_file2) or die "can't write to $out_file2\n";
	while(<$in>){
		chomp;

		my ($seq, $score_v1, $score_v2) = split(/\t/);
		print $out "$seq\t$tss{$seq}\t$score_v1\t$score_v2\n";

		# only tracking IMEter v2 scores for percentile information
		push(@all_scores, $score_v2);
	
		# convert distance to a bin number (based on window size)
		# and add score to an array of scores at this bin
		my $pos = int($tss{$seq}/$window_size) + 1;
		push(@{$scores_by_pos{$pos}}, $score_v2);
	}
	close($in);
	close($out);


	open($out, ">", $out_file3) or die "can't write to $out_file3\n";

	# print header
	print $out "Bin\tStart\tEnd\tMean_score\tN\tSum_score\tStd_err\tLower_CI\tUpper_CI\n";

	# loop over score data and print results
	foreach my $bin (sort {$a <=> $b} keys %scores_by_pos){

		my $sum = sum @{$scores_by_pos{$bin}};
		my $n = scalar @{$scores_by_pos{$bin}};
		my $mean = sprintf("%.3f", $sum / $n);
		$sum = int($sum);
		# need two alternative for calculating STDDEV
		# can't use n-1 if $n is 1, fudge it to leave n as 1.
		my $std_dev;
		if($n == 1){
			$std_dev = sqrt(sum(map {($_ - $mean) ** 2} @{$scores_by_pos{$bin}}) / ($n));
		} else{
			$std_dev = sqrt(sum(map {($_ - $mean) ** 2} @{$scores_by_pos{$bin}}) / ($n-1));
		}

		# want to calculate standard error of the mean and
		# show 95% confidence limits
		my $std_err  = sprintf("%.3f",$std_dev/sqrt($n));
		my $lower_CI = sprintf("%.3f", $mean - (1.96 * $std_err));		
		my $upper_CI = sprintf("%.3f", $mean + (1.96 * $std_err));

		my $min_coord = $window_size * ($bin-1) + 1;
		my $max_coord = $min_coord + $window_size -1;

		# stop printing output if we go too far into the transcript
		last if ($min_coord >= $cutoff);
		print $out "$bin\t$min_coord\t$max_coord\t$mean\t$n\t$sum\t$std_err\t$lower_CI\t$upper_CI\n";
	}
	close($out);

	@all_scores = sort {$a <=> $b} @all_scores;

	my $n = scalar @all_scores;
	open($out, ">", $out_file4) or die "can't write to $out_file4\n";

	my $percentile_target = 1;
	for (my $i = 0; $i < @all_scores; $i++){

		my $percent_of_scores = sprintf("%.2f", (($i + 1) / $n) * 100); 
		if($percent_of_scores >= $percentile_target){
			print $out "$percentile_target\t$all_scores[$i]\n";
			$percentile_target++;
		}
	}
	close($out);

}

exit;