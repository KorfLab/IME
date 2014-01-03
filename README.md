# IME #

A new code base for the IME project.


## Background to proc_phytozome.pl script ##

Script assumes that you have downloaded all data in a specific version of Phytozome to a local directory. E.g.
download ftp://ftp.jgi-psf.org//pub/compgen/phytozome/v9.0 to /Volumes/Scratch/Phytozome/v9.0.

The top level v9.0 directory should contain subdirectories for each species in Phytozome. The proc_phytozome.pl
script runs on one of these subdirectories at a time. It will search these directories to find the necessary files
that it needs (a genome file and a GFF3 file).


## Output of the proc_phytozome.pl script ##

By default the script outputs sequence information in a custom format for exons, introns and UTRs from all genes
that pass a set of criteria (to rule out bad gene annotations). It also logs how many errors were found with the 
annotations. 

Example of sequence output that is sent to STDOUT (only showing three lines):
	
	24093359        LOC_Os10g26570.1        exon    2231    2326    9       14      AT2G22480.1     ATTGTGCTTACACTGGAAAAATATGGAGTGAAAAACATTGTTGGAATACAGCATGGTTTCCGTGGATTTTTTGAGGATCATTTAGCAGAAGTGCCA
	24093359        LOC_Os10g26570.1        intron  541     630     1       13      AT2G22480.1     GTAAGGCACCTTTCCCTTCATTAGCTTTCTGCCACACAAGTTCTAGCAGTTCATGTTGTCTCAGCATTGATTTCTGATGGAACTTCGCAG
	24093359        LOC_Os10g26570.1        five_prime_UTR  3544    3577    1       2       AT2G22480.1     GGTAGTTGCATGATCGGCAACATCAATCTGAAGA

These fields correspond to:

1. Unique Phytozome transcript ID for each mRNA (taken from column 9 of GFF3 field)
2. Transcript ID from source database (tracked by Phytozome in a different file)
3. Type of sequence (either exon, intron, five_prime_UTR, or three_prime_UTR)
4. Start coordinate of feature (relative to TSS)
5. End coordinate of feature (relative to TSS)
6. Number of that feature for the current transcript (e.g. intron 1, intron 2 etc)
7. Total number of this feature for the current transcript. 
8. Best A. thaliana ortholog (from reciprocal BLASTP by Phytozome). For A. thaliana, the best O. sativa ortholog is shown. NA is shown if no ortholog.
9. Sequence

Example of error/statistics output (sent to STDERR):

	GENOME: /Volumes/Scratch/Phytozome/v9.0/Cpapaya/
	GENES: 27793
	EXONS: 112604
	CDS:   112604
	5'UTR: 27793
	3'UTR: 27793
	COMPLETE: 27793
	ERRORS: 4592
	KEPT: 23201
	Error types:
			internal_stop: 10
			max_intron: 105
			min_cds: 1384
			min_intron: 678
			non_canonical_splice: 3
			start_not_found: 1129
			stop_not_found: 2943
	Splice sites:
			GT..AG  84205
			GC..AG  601
			GT..GG  1
			GT..CA  1
			GG..AG  1
			TA..AG  1
			CT..AG  1
        

## Suggested usage for proc_phytozome.pl script ##

Minimal usage is just to specify a directory and redirect output to a file 

	./proc_phytozome.pl /Volumes/Scratch/Phytozome/v9.0/Cpapaya/ > Cpapaya.txt
	
Require genes to have 5' and 3' UTRs:

	./proc_phytozome.pl -u /Volumes/Scratch/Phytozome/v9.0/Cpapaya/ > Cpapaya_with_UTRs.txt
