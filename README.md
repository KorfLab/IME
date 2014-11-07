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
			There are 275 sequences present in the genome for this species
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
be skipped. Some species with very large number of sequences in their genomes may take a long
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

The first three parameters (k = 5, proximal and distal cutoff = 400 nt from TSS) were the same as used in the last published version of the IMEter. Previously we used data from  TAIR and slightly different filtering criteria. Here, the -c option restricts the script to only train from transcripts with 5' and 3' UTRs. The lack of the -i option means that we are not relying on just the primary isoform. If a gene has multiple transcripts, all our used to train the IMEter:


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

I used a simple Perl wrapper script to just call ime_trainer.pl many times:

```perl
	#!/usr/bin/perl
	use strict;
	use warnings FATAL => 'all';

	foreach my $file (glob("*_IME_intron.fa")){
		my $outfile = $file;
		$outfile =~ s/fa$/params/;
		my $command = "ime_trainer.pl -k 5 -p 400 -d 400 -c $file > $outfile";
		print "Running $command\n";
		system($command) && die "Can't run $command\n";
	}
```



## Exploring parameter space ##

If you want to try lots of different parameter combinations, you can use the 
ime_explore_parameter_space.pl script. This is a wrapper script around the ime_trainer.pl
script and changes several of the parameters in a set of nested loops. For each final
combination of parameters, it runs ime_trainer.pl and produces an output parameter file.

It then runs the main imeter.pl script (see below) using that parameter file to score a set 
of wild-type Arabidopsis intron sequences for which we know (experimentally) how much they
increase expression (relative to an intronless control). E.g.

	ime_explore_parameter_space.pl Athaliana_IME_intron.fa db_IME_Rose_WT_introns.fa final_IME_scores.tsv
	
+ The first argument is the input intron sequence data to train from.
+ The second argument is a file of 15 wild-type introns for which we know how much the intron
increases expression (if at all).
+ The third argument is the name of a final output file in which to put the IMEter scores
for all introns in the test file, under all parameter combinations.

By default, the file will test various sizes of k (3–7), different distances from the TSS,
use of partial or complete gene structures, and whether all splice isoforms of a gene 
are included in training, or just the primary isoform. It is expected that you would 
change the script to add/remove parameters as desired. In particular, change this line:

	foreach my $distance_metric (qw(position14 position15 coordinate200_400 coordinate300_300 coordinate400_400)){

The defaults are that it will test 1st vs 4th introns, 1st vs 5th, and introns that
start before 200, 300, or 400 bp from the TSS against introns that start after 400, 300, 
or 400 bp from the TSS.

The output file can be opened with Excel or R and subjected to correlation analysis to see
which parameter sets best explain the known levels of expression increase (which are stored
in the FASTA headers of the db_IME_Rose_WT_introns.fa file).



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
	  -i <pattern> Ignore any kmers from the scoring process if they match pattern
	  
	  
	  
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


# Getting data from Phytozome v10.0.4

1. First [set up Globus](http://genome.jgi.doe.gov/help/download.jsf#globus) account as per instructions.
2. Set up your Mac/PC to have an endpoint
3. Using Globus web interface, make a connection between endpoints
	4. Target 1: jgi#portal
	5. Path 1: By_Organism_Name/P/PhytozomeV10_[PhytozomeV10]/ 
	5. Target 2: kbradnam#bioinformatics
	6. Path 2: /Volumes/Phytozome/v10.0.4
	7. Give transfer a name (no periods), e.g. v_10_0_4_transfer
	8. Start transfer and wait



# Average intron length per species #

Wrote quick little script to process all of the *IME_intron.fa files to calculate mean intron length:

```bash
./mean.pl *_IME_intron.fa > mean_lengths.tsv
cat mean_lengths.tsv
Acoerulea       208945  467.4
Alyrata 77360   196.8
Athaliana       166960  163.9
Bdistachyon     135503  387.5
Brapa   162127  206.7
Cclementina     149331  350.3
Cpapaya 76379   453.6
Creinhardtii    149208  269.7
Crubella        110086  174.9
Csativus        139693  452.4
Csinensis       197920  352.3
Csubellipsoidea_C-169   62478   284.1
Egrandis        186570  439.3
Fvesca  128426  413.7
Gmax    330915  492.3
Graimondii      404366  352.5
Lusitatissimum  169894  242.2
Mesculenta      123412  420.5
Mguttatus_v1.1  107829  285.3
Mpusilla_CCMP1545       7992    189.0
Mpusilla_RCC299 4459    160.7
Mtruncatula     135596  437.2
Olucimarinus    1910    171.8
Osativa 167200  410.0
Ppatens 150985  275.1
Ppersica        99942   326.4
Ptrichocarpa    346711  382.2
Pvulgaris       145176  480.0
Sbicolor        109431  412.6
Sitalica        142458  338.6
Slycopersicum   113855  536.9
Smoellendorffii 76214   101.7
Stuberosum      133500  577.7
Tcacao  173226  485.0
Thalophila      126919  178.4
Vcarteri        81353   411.8
Vvinifera       114501  725.7
Zmays   187657  487.0

```


# Testing effect of excluding first region of transcript from IMEter training #

Modified ime_trainer.pl to accept a new -q option. This is used in conjuction with -p (proximal option). This new option allows an initial offset to be specified (in bp). Regions of introns that start anywhere between 1 bp from TSS up to the value provided by -q will be excluded. E.g. if -q is set to 100 and a 200 bp intron starts at 75 bp from the TSS, then kmers from the first 25 bp of the intron will not be included in the proximal data set.

Also trying a new -x option to exclude sequence from proximal introns that occur past the position specified by -p. E.g. if we set -p to 400, then by default we classify a 200 bp intron that starts at position 390 from the TSS as proximal, even though only 10 bp of the intron occurs in that range. The -x option just excludes all bases past the value of -p.

Now let's see what difference this makes as we increase -q (for A. thaliana). First
no value at all (just showing some of the parameter file, with only the top 10 kmers):

```bash
# build: ./ime_trainer.pl -k 5 -p 400 -d 400 -c Athaliana_IME_intron.fa
# introns counted:     138144
# bases clipped:       0
# proximal introns:    15353
# distal introns:      122791
CGCCG	1.22953648512629
CGATC	1.05461245354295
CGATT	1.03625398431636
CGGCG	0.982759636835548
TCGAT	0.964288118654301
TCCGA	0.892261809503596
TCGCG	0.890909071438347
GATCG	0.872558511203674
TAGGG	0.867923670376819
GGGTT	0.861602930056838
```

Now try -q = 50 (without -x):

```bash
# build: ./ime_trainer.pl -k 5 -p 400 -d 400 -c -q 50 Athaliana_IME_intron.fa
# introns counted:     138144
# bases clipped:       0
# proximal introns:    15353
# distal introns:      122791
CGCCG	1.06322721405146
GGGTT	0.977503611146258
TAGGG	0.952333887461336
AGGGT	0.927545796090563
GATCG	0.9193202969594
CGATT	0.914550077442204
CGGCG	0.883751912446206
CGACG	0.858184508444106
CGATC	0.837855159370795
CTGGG	0.831612794751938
```

This changes some of the top kmers that appear (and we see the non-CG-kmers: GGGTT, TAGGG, and AGGGT). Now we add the -x option into the mix (without -q):

```bash
# build: ./ime_trainer.pl -k 5 -p 400 -d 400 -c -x Athaliana_IME_intron.fa
# introns counted:     138144
# bases clipped:       1698391
# proximal introns:    15353
# distal introns:      122791
CGCCG	1.60524841755586
CGATC	1.59753071910004
CGATT	1.49984536467771
CGGCG	1.42480019742612
TCGAT	1.41630270342777
TCGCG	1.37613923833504
TCCGA	1.35693456044754
CGCGA	1.33072771219721
CCGAT	1.29872260925767
CTCCG	1.2769067957431

```bash

Addition of -x option increases log-odds scores quite a lot and changes the top 10 kmers again. So how does -q 50 and -x look:

```bash
# build: ./ime_trainer.pl -k 5 -p 400 -d 400 -c -q 50 -x Athaliana_IME_intron.fa
# introns counted:     138144
# bases clipped:       1698391
# proximal introns:    15353
# distal introns:      122791
CGATT	1.55043311963481
CGCCG	1.52927608597663
CGATC	1.52721337339612
GATCG	1.52563512967403
TCGCG	1.48513682662658
CGGCG	1.48287996440789
CGCGA	1.48097976431641
TCGAT	1.42604652272866
GGGTT	1.38024903705267
CTGGG	1.33097685014413

```bash

Try to see how -q and -x look along side the default results

```bash
Position	Kmer	Default score		With -q=50 and -x
#1			CGCCG	1.22953648512629	1.52927608597663 #2
#2			CGATC	1.05461245354295	1.52721337339612 #3
#3			CGATT	1.03625398431636	1.55043311963481 #1
#4			CGGCG	0.982759636835548	1.48287996440789 #6
#5			TCGAT	0.964288118654301	1.42604652272866 #8
#6			TCCGA	0.892261809503596
#7			TCGCG	0.890909071438347	1.48513682662658 #5
#8			GATCG	0.872558511203674	1.52563512967403 #4
#9			TAGGG	0.867923670376819
#10			GGGTT	0.861602930056838	1.38024903705267 #9
			CGCGA						1.48097976431641 #7
			CTGGG						1.33097685014413 #10

```

So changing the training program in this way generally boosts all log odds scores, but in relative terms it reduces the extra boost that was given to CGCCG. It also makes CGATT the highest kmer.

Now redo all of this with -q = 100. Will skip to the headline comparison:

```bash
Position	Kmer	Default score		With -q=100 and -x
#1			CGCCG	1.22953648512629	1.60986472320887 #1
#2			CGATC	1.05461245354295	
#3			CGATT	1.03625398431636	1.46511447708335 #7
#4			CGGCG	0.982759636835548	1.42560435903824 #10
#5			TCGAT	0.964288118654301	
#6			TCCGA	0.892261809503596
#7			TCGCG	0.890909071438347	1.4331938769954  #8
#8			GATCG	0.872558511203674	1.54163035481501 #2
#9			TAGGG	0.867923670376819	1.50126320812673 #4
#10			GGGTT	0.861602930056838	1.42600510249279 #9
			CGCGA						1.53087160432573 #3
			CGACG						1.48728417779348 #5
			AGGGT						1.47154824508    #6
```

So with -q at 100 bp, this makes more drastic changes compared to the old default values. What about -q=150?

```bash
Position	Kmer	Default score		With -q=150 and -x
#1			CGCCG	1.22953648512629	
#2			CGATC	1.05461245354295	1.35676232597733 #9
#3			CGATT	1.03625398431636	1.35395792474312 #10
#4			CGGCG	0.982759636835548	1.44595258982649 #4
#5			TCGAT	0.964288118654301	
#6			TCCGA	0.892261809503596
#7			TCGCG	0.890909071438347	1.38560982568778 #7
#8			GATCG	0.872558511203674	1.399509085695   #5
#9			TAGGG	0.867923670376819	1.46102377834578 #3
#10			GGGTT	0.861602930056838	
			CGCGA						1.61491451043138 #1
			CGACG						1.57132708389913 #2
			TTCGA						1.39159401248468 #6
			TCGAA						1.38258112063352 #8
```


Some of these kmers seems to dance up and down the list as you change the value of q. Now the only thing to do is to see whether any of these parameter options do a better job of explaining known variation in expression enhancing abilities of introns that we have tested. Will make 7 parameter files:

```bash
./ime_trainer.pl -k 5 -p 400 -d 400 -c  Athaliana_IME_intron.fa  > At.params
./ime_trainer.pl -k 5 -p 400 -d 400 -c  -q 150 -x Athaliana_IME_intron.fa  > At_q150_x.params
./ime_trainer.pl -k 5 -p 400 -d 400 -c  -q 100 -x Athaliana_IME_intron.fa  > At_q100_x.params
./ime_trainer.pl -k 5 -p 400 -d 400 -c  -q 50 -x Athaliana_IME_intron.fa  > At_q50_x.params
./ime_trainer.pl -k 5 -p 400 -d 400 -c  -q 50  Athaliana_IME_intron.fa  > At_q50.params
./ime_trainer.pl -k 5 -p 400 -d 400 -c  -q 100  Athaliana_IME_intron.fa  > At_q100.params
./ime_trainer.pl -k 5 -p 400 -d 400 -c  -q 150  Athaliana_IME_intron.fa  > At_q150.params
```

Now for each of these parameter files we will run the IMEter aginst the test set of 15 wildtype introns:

```bash
./imeter.pl -m At_q150.params db_IME_Rose_WT_introns.fa
```

After some experimentation, q=125 without -x option, and without implicit clipping with -q option, gives better results (as determined by final correlation):

```bash
Position	Kmer	Default score		With -q=125
#1			CGCCG	1.22953648512629	1.06633217188795 #1
#2			CGATC	1.05461245354295	0.876360044391628 #3
#3			CGATT	1.03625398431636	0.883135549901786 #2
#4			CGGCG	0.982759636835548	0.791173522227393 #8
#5			TCGAT	0.964288118654301	0.805111240054373 #6
#6			TCCGA	0.892261809503596	0.738486076412956 #12
#7			TCGCG	0.890909071438347	0.859387867265468 #4
#8			GATCG	0.872558511203674	0.719652200853388 #15
#9			TAGGG	0.867923670376819	0.794060471074419 #7
#10			GGGTT	0.861602930056838	0.819760174260988 #5
#11 		TTCGA   0.853614225351389	0.747575198743333 #10
#12			CGCGA   0.847351833937376	0.784740501033 #9
#13			CGTCG   0.830087669049207
#14			CGACG   0.829707674657411
#15			AGGGT   0.822745275109967	0.744162866759895 #11
#16			ATCGA   0.811812039338427	0.713130865185136 #16
#17			CGAAT   0.810244108724747	0.729560512894274 #13
#18			CCGAT   0.806683868159494	0.675179865369885 #18
#19			TCGAA   0.789184771054985	0.726436511010194 #14
#20			AATCG   0.787173377488944	0.70467421301563  #17
			CTCCG   					0.665565846584791 #19
			CTGGG   					0.66020922614789  #20
```

But I find the enriched kmers iun this set to be so similar to the original set, and CGCCG still scores so much higher than other kmers.

