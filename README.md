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

