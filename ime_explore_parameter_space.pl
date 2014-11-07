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
use List::Util qw(sum);
use Getopt::Std;
 
getopts('h');

die "
usage: $0 [params] <At intron file> <test intron file> <file of expression values>

First file should be parameter file produced by ime_trainer.pl 
Second file should be taken from dbIME
Optional third argument should be a file containing expression values for each intron in second file

params:
  -h       this help
" if (@ARGV < 2 or @ARGV > 3 or $opt_h);

my ($intron_file, $test_seq_file, $expression_file)  = @ARGV;

my $script1 = "ime_trainer.pl";
my $script2 = "imeter.pl";

print "\n\n";

# set up list of various distance metrics that we can test
my @distance_metrics;

# integer positions for introns
push @distance_metrics, qw(position12 position13 position14 position15 position16 position17 position18);
push @distance_metrics, qw(position23 position34 position45);

# absolute distances for -p and -d parameters
push @distance_metrics, qw(coordinate100_400 coordinate200_400 coordinate300_400 coordinate400_400);
push @distance_metrics, qw(coordinate125_450 coordinate125_500 coordinate125_550);
push @distance_metrics, qw(coordinate100_600 coordinate200_600 coordinate300_600 coordinate400_600 coordinate500_600 coordinate600_600);
push @distance_metrics, qw(coordinate100_800 coordinate300_800 coordinate500_800 coordinate700_800);

# include offset for -q option to exclude introns very close to TSS
push @distance_metrics, qw(coordinate50_400_400 coordinate75_400_400 coordinate125_400_400);
push @distance_metrics, qw(coordinate150_400_400 coordinate175_400_400 coordinate200_400_400);


push @distance_metrics, qw(coordinate125_300_300 coordinate125_500_500 coordinate125_600_600);
push @distance_metrics, qw(coordinate125_700_700 coordinate125_500_700 coordinate125_400_900);
push @distance_metrics, qw(coordinate125_300_900 coordinate125_300_1000 coordinate125_400_1000);
push @distance_metrics, qw(coordinate125_500_1000 coordinate125_600_1200 coordinate125_750_1500);
push @distance_metrics, qw(coordinate150_500_1000 coordinate150_1000_1500 coordinate150_1000_2000);

push @distance_metrics, qw(coordinate100_400_400 coordinate200_400_400 coordinate300_400_400);
push @distance_metrics, qw(coordinate100_500_500 coordinate200_500_500 coordinate300_500_500);
push @distance_metrics, qw(coordinate100_600_600 coordinate200_600_600 coordinate300_600_600);
push @distance_metrics, qw(coordinate100_700_700 coordinate200_700_700 coordinate300_700_700);


my $parameter_set = 0;

# will store all correlation results in a hash
my %correlations;



# vary size of k used for kmer counting
for (my $k = 3; $k <= 7; $k++){

	# just some examples of distance parameters that could be used. Substitute your own here
	# if using coordinate system with 2 numbers, these will be used for values of -p and -d
	# if using coordinate system with *3* numbers, these will be used for values of -q, -p and -d 
	# where -q must be less than -p
	foreach my $distance_metric (@distance_metrics){
	
		# require complete or partial gene structures 
		foreach my $structure (qw(complete partial)){

			# only use primary splice isoform, or use all isoforms (secondary)
			foreach my $isoform (qw(primary secondary)){

				# finally want a way of using or not using -x option for clipping intron sequences
				foreach my $x_option (qw(yes no)){
				
					$parameter_set++;
					
					print "Parameter set $parameter_set: k=$k distance_method=$distance_metric ";
					print "structure_requirement=$structure isoforms=$isoform x_option = $x_option\n";
					print "----------------------------------------------------------------------------------\n";
				
					my $command = "$script1 -k $k ";
					my $out_file = "Ath_IME_k$k";
					if($distance_metric =~ m/coordinate\d+_\d+_\d+/){
						my ($q, $p, $d) = $distance_metric =~ m/(\d+)_(\d+)_(\d+)/;
						$command .= "-q $q -p $p -d $d ";
						$out_file .= "_coords${q}_${p}_$d";
					} 				

					elsif($distance_metric =~ m/coordinate\d+_\d+/){
						my ($p, $d) = $distance_metric =~ m/(\d+)_(\d+)/;
						$command .= "-p $p -d $d ";
						$out_file .= "_coords0_${p}_$d";
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

					# only use x-option with some of the distance parameters
					if($x_option eq "yes" and $distance_metric =~ m/coordinate/){
						$command .= "-x ";
						$out_file .= "_x";
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

					# can now calculate correlation for this parameter set
					# but only if we have expression data
					if ($expression_file){
					
						# make correlations for IMEter v1 and v2
						my $correlation1 = correlation($out_file_2, 'v1');
						my $correlation2 = correlation($out_file_2, 'v2');

						print "Correlation (r): v1 = $correlation1, v2 = $correlation2\n";
						
					}
					print "\n\n";
				}
			}
		}	
	}
}


# if we have been tracking expression scores we can print out all of the results
# sorted by correlation
exit(0) unless ($expression_file);

my $outfile = "IME_correlations.tsv";
open(my $out, ">", $outfile) or die "Can't write to $outfile\n";
foreach my $correlation (sort {$correlations{$a} <=> $correlations{$b}} keys %correlations){
	print $out "$correlation\t$correlations{$correlation}\n";
}
close($out);

exit(0);


####################################
#
# SUBROUTINES
#
####################################


sub correlation{

	my ($imeter_score_file, $version) = @_;
	my $prefix = $imeter_score_file;
	$prefix =~ s/\.scores//;
	my @expression_values = read_expression_data($expression_file);

	my @scores;

	
	open(my $in, "<", $imeter_score_file) or die "Can't open $imeter_score_file\n";
	while(my $line = <$in>){
		chomp($line);
		my ($id, $v1, $v2) = split(/\s+/, $line);
		push(@scores, $v1) if ($version eq 'v1'); 
		push(@scores, $v2) if ($version eq 'v2'); 
	}
	
	close($in);

	# now calculate pearson correlation coefficient from IMEter v2 scores
	my $n = @expression_values;
	my $mean1 = sum(@expression_values) / $n;
	my $mean2 = sum(@scores) / $n;
	
	my $covariance = 0;
	my $variance1 = 0;
	my $variance2 = 0;
	
	for (my $i = 0; $i < @expression_values; $i++){
		$covariance += (($expression_values[$i] - $mean1) * ($scores[$i] - $mean2));
		$variance1  += (($expression_values[$i] - $mean1)**2);
		$variance2  += (($scores[$i] - $mean2)**2);	
	}
	
	if($variance1 == 0 or $variance2 == 0){
		$correlations{$prefix} = 0.000;
		return(sprintf("%.3f",0));		
	}
	my $r = $covariance / (sqrt($variance1) * sqrt($variance2));

	# store correlation in hash and return formatted correlation
	my $index = $prefix . "_" . $version;
	$correlations{$index} = $r;
	return(sprintf("%.3f",$r));		
	
}

sub read_expression_data{
	my ($file) = @_;
	
	my @values;
	open(my $in, "<", $file) or die "Can't open $file\n";
	while(my $line = <$in>){
		chomp($line);
		push(@values, $line);
	}
	
	close($in);
	return(@values);
}
# now we have all scores in separate files, we can combine them into one main table

exit(0);
