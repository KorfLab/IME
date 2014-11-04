#!/usr/bin/perl
#
# ime_explore_parameter_space.pl
#
# A script to run the ime_trainer.pl script under many different conditions
# and test results against known Arabidopsis introns and their de facto
# IMEter scores with IMEter v2.0. This is to address two things:
#
# 1) Does our new training regime with Phytozome data give broadly comparable 
# IMEter scores to what we did previously using TAIR data and slightly 
# different filtering rules?
#
# 2) Are there different parameter combinations that better explain known
# expression increase values than what we used previously?
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use vars qw($opt_h);
use Getopt::Std;

getopts('h');

die "
usage: $0 [params] <At intron file> <test intron file (with expression scores)> <output TSV filename>

First file should be produced by ime_trainer.pl 
Second file should be taken from dbIME
Third argument is final name of output file to create

params:
  -h       this help
" if @ARGV != 3 or $opt_h;

my ($intron_file, $test_seq_file, $table_filename)  = @ARGV;

my $script1 = "ime_trainer.pl";
#my $script2 = "the_imeter.pl";
my $script2 = "imeter.pl";

print "\n\n";

# vary size of k used for kmer counting
for (my $k = 3; $k < 8; $k++){

	# just some examples of distance parameters that could be used. Substitute your own here
	foreach my $distance_metric (qw(position14 position15 coordinate200_400 coordinate300_300 coordinate400_400)){

		# require complete or partial gene structures 
		foreach my $structure (qw(complete partial)){

			# only use primary splice isoform, or use all isoforms (secondary)
			foreach my $isoform (qw(primary secondary)){

				print "Parameters: k=$k distance_method=$distance_metric ";
				print "structure_requirement=$structure isoforms=$isoform\n";
				print "----------------------------------------------------------------------------------\n";
				
				my $command = "$script1 -k $k ";
				my $out_file = "Ath_IME_k$k";
				if($distance_metric =~ m/coordinate\d+_\d+/){
					my ($p, $d) = $distance_metric =~ m/(\d+)_(\d+)/;
					$command .= "-p $p -d $d ";
					$out_file .= "_coords${p}_$d";
				} 				
				elsif($distance_metric =~ m/position\d\d/){
					my ($p, $d) = $distance_metric =~ m/position(\d)(\d)/;
					$command .= "-P $p -D $d ";
					$out_file .= "_position${p}$d";						
				} else{
					die "error: $distance_metric";
				}

				if($structure eq "complete"){
					$command .= "-c ";
					$out_file .= "_complete";
				} else{
					$out_file .= "_incomplete";				
				}

				if($isoform eq "primary"){
					$command .= "-i ";
					$out_file .= "_primary";
				} else{
					$out_file .= "_secondary";				
				}
								
				$out_file .= ".params";

				$command .= "$intron_file > $out_file";
				
				# can now run command to train IMEter with some particular set of
				# parameters. Don't run if param file already exists.
				if (-e $out_file){
					warn "$out_file already exists, skipping...\n";
				} else{
					warn "Running: $command\n\n";
					system($command) && die "Can't run $command\n";				
				}
				
				
				# now run actual IMEter script to score the test set of introns
				# using supplied parameter file
				my $out_file_2 = $out_file;

				$out_file_2 =~ s/\.params//;
				$out_file_2 .= ".scores";
				
				$command = "$script2 -m $out_file $test_seq_file > $out_file_2";

				if(-e $out_file_2){
					warn "$out_file_2 already exists, skipping...\n";
				} else{
					warn "Running: $command\n";	
					system("echo \"Intron\t$out_file\" >> $out_file_2") && die "Can't echo";
					system($command) && die "Can't run $command\n";				

				}

				print "\n\n";
			}
		}	
	}
}

# now we have all scores in separate files, we can combine them into one main table

# how many introns were there in test file?
my $intron_count = `grep -c \">\" $test_seq_file`; chomp $intron_count;

my @ime_score_files = glob("*.scores");

# hash to store data
my %scores;

foreach my $file (@ime_score_files){
	warn "Processing $file\n";
	my $key = $file;
	$key =~ s/\.scores//;
	
	open(my $in, "<", $file) or die "Can't read $file\n";
	
	while(<$in>){
		chomp;
		my ($intron, $score_v1, $score_v2) = split;
#		print "Intron $intron has an IMEter v2.0 score of $score_v2\n";
		push(@{$scores{$key}}, $score_v2);
	}
	close($in);
}

# now print table

open(my $out, ">", $table_filename) or die "can't write to $table_filename\n";
# first print header row
foreach my $key (sort keys %scores){
	print $out "$key\t";
}
print $out "\n";

# now print remaining rows

for (my $i=0; $i < $intron_count; $i++){
	foreach my $key (sort keys %scores){
		print $out "${$scores{$key}}[$i]\t";
	}
	print $out "\n";
}
close($out);


