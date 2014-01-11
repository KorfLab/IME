# IME #

A new code base for the IME project.


## Background to proc_phytozome.pl script ##

Script assumes that you have downloaded all data in a specific version of Phytozome to a 
local directory. E.g. download ftp://ftp.jgi-psf.org//pub/compgen/phytozome/v9.0 
to /Volumes/Scratch/Phytozome/v9.0.

The top level v9.0 directory should contain subdirectories for each species in Phytozome.  
The proc_phytozome.pl script runs on the top-level directory and will process each species
it comes across in the subdirectories (though it skips the 'early_release' directory.  
For each species-level subdirectory it will find the necessary files that it needs 
(a genome file, a GFF3 file, and a Phytozome annotation text file).


## Output of the proc_phytozome.pl script ##

By default the script outputs a separate FASTA-formatted file for exons, introns and UTR
exons (5' and 3') from all genes that pass a set of criteria (to rule out bad gene 
annotations). It also logs to a separate how many errors were found with the annotations. 

The FASTA header of the output file stores a lot of meta information. E.g.

	>IME_Acoerulea_1 TYPE=intron POS=1/1 COORDS=1594-3461 ID1=22022986 ID2=Aquca_023_00143.1 
	ISOFORM=primary STRUCTURE=5-3 ORTHOLOG=AT1G18740.1

These fields correspond to:
				
1. IME_Acoerulea_1 - a unique FASTA ID formed by the script, the syntax can be changed by
using the -s (stub) option. Successive features have an incremented numerical suffix. This 
makes all features have a unique FASTA ID across all species in Phytozome.
2. TYPE=intron - a simple flag to indicate intron, exon, five_prime_utr, or three_prime_utr
3. POS=1/1 - in this case, this indicates 1st intron out of 1 intron total (for the 
current transcript)
4. COORDS=1594-3461 - coordinates of feature relative to transcription start site (TSS).
5. ID1=22022986 - the unique Phytozome transcript ID for each mRNA (taken from column 9 of 
GFF3 field)
6. ID2=Aquca_023_00143.1 - the transcript ID from source database (tracked by Phytozome 
in a different file)
7. ISOFORM=primary - type of transcript: primary (longest isoform) or secondary
8. STRUCTURE=5-3 - transcript structure: does it have a 5' UTR, a 3' UTR, both, or neither?
9. ORTHOLOG=AT1G18740.1 - best A. thaliana ortholog (from reciprocal BLASTP by Phytozome). 
For A. thaliana, the best O. sativa ortholog is shown. NA is shown if no ortholog.


Example of error/statistics output file:

	GENOME: /Volumes/Scratch/Phytozome/v9.0//Acoerulea
	GENES: 41063
	GENES_WITH_ERRORS: 3587
	GENES_KEPT: 37476
	EXONS: 274774
	CDS:   256568
	5'UTR: 41063
	3'UTR: 41063
	COMPLETE: 41063
	Error types:
			max_3utr: 1696
			max_5utr: 135
			min_cds: 69
			min_intron: 483
			non_canonical_splice: 68
			start_not_found: 1053
			stop_not_found: 241
	Splice sites:
			GT..AG  230560
			GC..AG  2921
			AT..AC  162
			G..G    16
			AT..AA  6
			T..T    4
			TA..GT  3
			C..C    3
			AT..GG  3
			AA..TC  2
			A..A    2
			TC..TA  1
        

## Suggested usage for proc_phytozome.pl script ##

Minimal usage is just to specify a directory and redirect output to a file 

	./proc_phytozome.pl /Volumes/Scratch/Phytozome/v9.0/
	
Verbose mode (-v) will give you information about the scripts's progress. E.g.

	proc_phytozome.pl -v /Volumes/Scratch/Phytozome/v9.0/
	Processing data for Acoerulea
		Step 1: processing GFF file
		Step 2: creating intron information from exon data
		Step 3: processing genome sequence file
		Step 4: looking for errors in annotations
		Step 5: extracting extra information from Phytozome annotation file
		Step 6: final output
			Creating Acoerulea_IME_exon.fa
			Creating Acoerulea_IME_intron.fa
			
May want to skip some problem species by using the -i option:

	proc_phytozome.pl -v -i "Mdomestica Pvirgatum Rcommunis" /Volumes/Scratch/Phytozome/v9.0/
	
Any species matching the text in the space-separated list (provided by -i option) will
be skipped. Some species with very large number of sequences in their genomes take a long
time to process.


## Training IMEter from an sequence file ##

With the files produced by proc_phytozome.pl, we can now generate an IMEter parameter 
file. These files simply contain the log-odd ratios of a set of k-mers in relation to their
presence in a) a set of proximal introns (or exons, 5' UTRs etc) and b) a set of distal
introns.

There are many different ways of defining proximal and distal (1st introns vs 5th introns, 
introns that start < 400 nt from the TSS vs introns that start > 600 nt etc). The ime_trainer.pl
script can make these considerations. You can specify:

+ k (the word size to use)
+ proximal cutoff using nt coordinates
+ distal cutoff using nt coordinates
+ proximal cutoff using intron position
+ distal cutoff using intron position
+ require complete transcripts (those with both 5' and 3' UTRs annotated)
+ only use introns from primary isoform (defined by Phytozome as the longest transcript)

These options let you a train an IMEter based on your specific preferences, or based on
limitations in the data. E.g. in less well annotated genomes, there may be fewer annotated
UTRs, which means that coordinate based cutoffs may be less accurate (because you may not
know the true start of the transcript).

We trained the latest version of the IMEter like so:

	ime_trainer.pl -k 5 -p 400 -d 400 -c Athaliana_IME_intron.fa > Ath_IME_k5_400_400_complete.params

The first three parameters (k = 5, proximal and distal cutoff = 400 nt from TSS) were the
same as used in the last published version of the IMEter. Previously we used data from 
TAIR and slightly different filtering criteria. Here, the -c option restricts the script to
only train from transcripts with 5' and 3' UTRs. The lack of the -i option means that
we are not relying on just the primary isoform. If a gene has multiple transcripts, all our
used to train the IMEter:


The final parameter file will look something like this:

	# IME training log
	# ----------------
	#
	# Fri Jan 10 11:10:51 PST 2014
	# build: /Users/keith/Work/Code/IME/ime_trainer.pl -c -k 5 -p 400 -d 400 Crubella_IME_intron.fa
	# transcripts:               20017
	# transcripts not counted:   9082
	# introns counted:     59702
	# introns not counted: 0
	# proximal introns:    5929
	# distal introns:      53773
	# skipped kmers:       169 maximum = NNNNN (7818)
	#
	# IME parameters
	# --------------
	CGATT   1.05876003163932
	CGATC   1.03513727248767
	TCGAT   0.989253683124108
	CGACG   0.958234172791253
	GATCG   0.953133431986251
	TCGCG   0.906243637469971
	CGGCG   0.896922851358586
	TTCGA   0.872332381929915
	CCGAT   0.865947978284506
	TCCGA   0.843627296804883

After some statistics about the training run, each line shows the log-odds ratio of the 
frequency of the specified kmer in the proximal data set vs the distal data set. This 
parameter file is required by the main IMEter script.


## Generating many parameter files ##

If you want to try lots of different parameter combinations, you can use the 
ime_explore_parameter_space.pl script. This is a wrapper script around the ime_trainer.pl
script and changes several of the parameters in a set of nested loops. For each final
combination of parameters, it runs ime_trainer.pl and produces an output file.

It then runs the main imeter.pl script (see below) using that parameter file to score a set 
of wild-type Arabidopsis intron sequences for which we know (experimentally) how much they
increase expression (relative to an intronless control). E.g.

	ime_explore_parameter_space.pl Athaliana_IME_intron.fa db_IME_all_WT_introns.fa final_IME_scores_all_WT.tsv
	
+ The first argument is the input data to train from.
+ The second argument is the file of 15 introns for which we know how much the intron
increases expression (if at all)
+ The third argument is the name of a final output file in which to put the IMEter scores
for all introns in the test file, under all parameter combinations.


## Running the IMEter ##

You always need to use a suitable set of IME parameters with the imeter.pl script:
	
	imeter.pl -m Ath_IME_k5_400_400_complete.params db_IME_all_WT_introns.fa

There are many options to further change how the IMEter scores a sequence:

	Options:
	  -w <int>     window size nt    [default 7]
	  -s <int>     step size  nt     [default 1]
	  -d <int>     donor sequence to clip [default 5]
	  -a <int>     acceptor sequence to clip [default 10]
	  -g <int>     minimum gap allowed between high scoring windows [default 5]
	  -c <float>   threshold value to decide what is a high scoring window [default 1.2]
	  -f <int>     weighting factor to penalize peaks that are further away [default 200]
	  -m <file>    IMEter parameter file [default is to use embedded pentamers]
	  -r           calculate score for reverse strand
  	  -p           A file of IMEter scores at each percentile (to include in output)
  	  -o		   print GFF info of each peak
	  
	  
	  
	  
## Calculating scores at different distances from TSS ##

There is another wrapper script (ime_distance_info.pl) that runs the imeter.pl script to 
produce some more useful information about IMEter scores at different distances from
the TSS. This script does 3 things:

1. Script will look in the specified directory for IME intron files ('*.IME_intron.fa')
and IME parameter files ('*.IME_intron.params'). For each intron file, it will 
run the IMEter (imeter.pl) on the specified file using the specified parameters
to produce a '*.IME_intron.ime1' file containing IMEter scores for each intron.
2. Script will then take the '*.IME_intron.ime1' file and for each intron, will extract
the distance from the transcription start site (TSS) from the FASTA file and include this
information to a new output file ('*.IME_intron.ime2').
3. The script will produce a third output file ('*.IME_intron.ime3') that averages
all intron scores in binned windows of distance from the TSS. The default parameter
takes windows of 200 bp from the TSS. In addition to calculating the average IMEter score
for all introns that start in that window, the script outputs 95% confidence limits for
the mean value.
4. Finally, the script calculates what IMEter v2.1 score is at each percentile. This
information is written to the fourth file ('*.IME_intron.ime4'). This is useful to give
an overview of what scores are 'high' or 'low'.
	
Suggested usage:

	ime_distance_info.pl -w 100 -c 2000 -i ~/bin/imeter.pl ime_data_directory

The first two output files may be of less use and can probably be discarded. The 3rd file
is another tab-separated values file, that will look like this:

	Bin     Start   End     Mean_score      N       Sum_score       Std_err Lower_CI        Upper_CI
	1       1       100     4.39    2625    11528   0.09    4.21    4.57
	2       101     200     4.56    4512    20552   0.07    4.42    4.70
	3       201     300     3.28    6106    20013   0.05    3.18    3.38
	4       301     400     2.57    7209    18497   0.04    2.49    2.65
	5       401     500     2.14    7008    15023   0.03    2.08    2.20
	6       501     600     1.83    6712    12266   0.03    1.77    1.89
	7       601     700     1.64    6772    11102   0.03    1.58    1.70
	8       701     800     1.55    6831    10619   0.03    1.49    1.61
	9       801     900     1.50    6560    9836    0.03    1.44    1.56
	10      901     1000    1.55    6476    10033   0.03    1.49    1.61
	11      1001    1100    1.54    6484    9996    0.03    1.48    1.60
	12      1101    1200    1.52    6162    9371    0.03    1.46    1.58
	13      1201    1300    1.58    6102    9664    0.03    1.52    1.64
	14      1301    1400    1.56    5612    8754    0.03    1.50    1.62
	15      1401    1500    1.55    5390    8364    0.03    1.49    1.61
	16      1501    1600    1.63    5120    8359    0.03    1.57    1.69
	17      1601    1700    1.59    4785    7620    0.03    1.53    1.65
	18      1701    1800    1.67    4525    7553    0.04    1.59    1.75
	19      1801    1900    1.60    4196    6695    0.04    1.52    1.68
	20      1901    2000    1.62    3951    6419    0.04    1.54    1.70
	
This shows the average IMEter score for the first 20 windows (size = 100 bp) from the TSS.
This data — for A. thaliana — shows that IMEter scores drop very quickly after the first 
400 bp or so, and then reach a uniform level. The average scores are fairly low due to 
presence of thousands of introns with zero IMEter scores.	

The end of the fourth output file, will look something like this (for A. thaliana):

	90      7.54
	91      7.98
	92      8.50
	93      9.17
	94      9.96
	95      10.96
	96      12.13
	97      13.67
	98      16.05
	99      19.71
	100     39.61

This can be simply be interpreted as follows:

+ Any IMEter score (v2.1) greater than 39.61 puts that intron into the 100th percentile of
all intron scores. 
+ A score between 7.98 and 8.49 would put an intron into the 91st percentile.
+ Note, this file does not show the absolute highest score observed. 
+ These files can be included as an argument to the imeter.pl script (with the -p) option
to include this information as an extra column in the IMEter output.