## Exploring parameter space ##

	ime_explore_parameter_space.pl Athaliana_IME_intron.fa db_IME_all_WT_introns.fa final_IME_scores_all_WT.tsv

Tried 264 different parameter combinations. Some observations:

1) Best correlations by k value (3–8):

Ath_IME_k6_coords300_300_complete_primary 0.846142953
Ath_IME_k7_coords300_300_incomplete_secondary 0.840725093
Ath_IME_k4_coords300_300_incomplete_primary 0.837650818
Ath_IME_k5_coords300_300_incomplete_primary 0.837562068
Ath_IME_k3_coords300_300_incomplete_primary 0.818041013


2) Best correlations by intron positions (1v2 through to 1v5)

Ath_IME_k8_position12_complete_primary 0.820445601
Ath_IME_k8_position15_complete_secondary 0.818214099
Ath_IME_k8_position12_incomplete_primary 0.816309691
Ath_IME_k8_position13_incomplete_primary 0.811909421
Ath_IME_k8_position15_incomplete_primary 0.810679028


3) Best correlations for distance (fixing k=6, complete, primary)

Ath_IME_k6_coords300_300_complete_primary 0.846142953
Ath_IME_k6_coords200_400_complete_primary 0.828969165
Ath_IME_k6_coords400_400_complete_primary 0.826970673
Ath_IME_k6_coords300_500_complete_primary 0.805049303
Ath_IME_k6_coords500_500_complete_primary 0.799511682
Ath_IME_k6_coords400_1000_complete_primary 0.764247735
Ath_IME_k6_coords500_1500_complete_primary 0.736115214

Then tried exploring many more parameters that involved changing distance (664 parameters
in total). Top 10 parameter combinations:

Ath_IME_k7_coords75_75_complete_primary 0.864809105
Ath_IME_k6_coords125_125_complete_primary 0.861917072
Ath_IME_k6_coords50_50_complete_secondary 0.8618146
Ath_IME_k6_coords75_75_complete_primary 0.859450469
Ath_IME_k6_coords125_125_incomplete_secondary 0.859295113
Ath_IME_k6_coords125_125_complete_secondary 0.85889568
Ath_IME_k6_coords100_200_complete_secondary 0.858740162
Ath_IME_k6_coords150_150_complete_primary 0.858696257
Ath_IME_k7_coords50_50_complete_primary 0.858660677
Ath_IME_k5_coords125_125_complete_secondary 0.858627399


Doesn't seem to help to use two different cutoffs for proximal and distal:

Ath_IME_k5_coords150_150_complete_secondary 0.844814716
Ath_IME_k5_coords150_200_complete_secondary 0.843852888
Ath_IME_k5_coords150_250_complete_secondary 0.837171073
Ath_IME_k5_coords150_300_complete_secondary 0.823089078
Ath_IME_k5_coords150_350_complete_secondary 0.818899223
Ath_IME_k5_coords150_450_complete_secondary 0.817451999
Ath_IME_k5_coords150_400_complete_secondary 0.816109251
Ath_IME_k5_coords150_500_complete_secondary 0.812512429

Ath_IME_k5_coords200_200_complete_secondary 0.829438923
Ath_IME_k5_coords200_250_complete_secondary 0.822225746
Ath_IME_k5_coords200_300_complete_secondary 0.819685147
Ath_IME_k5_coords200_350_complete_secondary 0.814450373
Ath_IME_k5_coords200_400_complete_secondary 0.811879036
Ath_IME_k5_coords200_450_complete_secondary 0.807568859

You see mostly similar results for other values of K and other permutations of 
complete/incomplete and primary/secondary, though not always.

At this point I decided to explore many more distance permutations, to see whether an
even closer cutoff than 150/150 would improve things. Here are the best 10 correlations 
from 920 different combinations of parameters:

Ath_IME_k7_coords60_60_complete_primary 0.882694056
Ath_IME_k7_coords70_70_complete_primary 0.870923488
Ath_IME_k6_coords70_70_complete_primary 0.868670409
Ath_IME_k7_coords75_75_complete_primary 0.864809105
Ath_IME_k6_coords140_140_complete_primary 0.863643425
Ath_IME_k7_coords110_110_complete_primary 0.862574621
Ath_IME_k6_coords125_125_complete_primary 0.861917072
Ath_IME_k6_coords50_50_complete_secondary 0.8618146
Ath_IME_k7_coords120_120_complete_primary 0.861680914
Ath_IME_k6_coords130_130_complete_primary 0.861677519

So it seems that an even closer distance for the proximal cutoff improves things, but not 
too close. If we keep other variables constant and just look at cutoffs which use the
same distance for proximal and distal, this is what we see (results organized by distance
not by score):

Ath_IME_k7_coords20_20_complete_primary 0.6249875 (27 proximal introns)
Ath_IME_k7_coords30_30_complete_primary 0.694340496 (91 proximal introns)
Ath_IME_k7_coords40_40_complete_primary 0.834509502 (171 proximal introns)
Ath_IME_k7_coords50_50_complete_primary 0.858660677 (295 proximal introns)
Ath_IME_k7_coords60_60_complete_primary 0.882694056
Ath_IME_k7_coords70_70_complete_primary 0.870923488
Ath_IME_k7_coords75_75_complete_primary 0.864809105
Ath_IME_k7_coords80_80_complete_primary 0.855402744
Ath_IME_k7_coords90_90_complete_primary 0.846729249
Ath_IME_k7_coords100_100_complete_primary 0.838913918
Ath_IME_k7_coords110_110_complete_primary 0.862574621
Ath_IME_k7_coords120_120_complete_primary 0.861680914
Ath_IME_k7_coords125_125_complete_primary 0.857973515
Ath_IME_k7_coords130_130_complete_primary 0.857834323
Ath_IME_k7_coords140_140_complete_primary 0.841992991
Ath_IME_k7_coords150_150_complete_primary 0.851001023
Ath_IME_k7_coords170_170_complete_primary 0.840904986
Ath_IME_k7_coords175_175_complete_primary 0.842199618
Ath_IME_k7_coords180_180_complete_primary 0.829647351
Ath_IME_k7_coords190_190_complete_primary 0.837363673
Ath_IME_k7_coords200_200_complete_primary 0.835334444
Ath_IME_k7_coords300_300_complete_primary 0.828851407
Ath_IME_k7_coords350_350_complete_primary 0.830569894
Ath_IME_k7_coords400_400_complete_primary 0.829410215
Ath_IME_k7_coords500_500_complete_primary 0.816330004


At this point it seems that we are just training the IMEter to pick up CG dinucleotides 
in early introns that may be all due to CpG islands. E.g. the top kmers in the parameter
set with the highest correlation (Ath_IME_k7_coords60_60_complete_primary) look like this:

ACGCGCG 5.75977379373248
GGGCGCG 5.75977379373248
GCGCGGC 5.75977379373248
CGACGCG 5.34473629445363
GGCGCGC 5.17481129301132
CGGGCGC 5.17481129301132
CCGGACG 5.02280819956627
CGCCCCG 4.95241887167487

But we know that the IMEter motif isn't CG rich. This gave the idea of trying to run
the IMEter in a mode where we just ignore any kmer which contains CGs. I added a new
-i option to the IMEter to allow you to ignore kmers that match a specific pattern. E.g.

	imeter.pl -m Ath_IME_k5_coords400_400_complete_primary.params -i CG db_IME_Rose_WT_introns.fa

In this instance, any kmers that contain CG do not count to the IMEter (v1 or v2) score. 
So we can repeat the above search of various IMEter parameters and ask whether we find
any combination that still explains the known expression values. This produces a new top 10:

Ath_IME_k6_coords125_125_incomplete_secondary_noCG 0.834920823
Ath_IME_k4_coords250_250_incomplete_primary_noCG 0.829442815
Ath_IME_k4_coords160_160_incomplete_primary_noCG 0.826446838
Ath_IME_k7_coords70_70_complete_primary_noCG 0.822882076
Ath_IME_k4_coords130_130_incomplete_secondary_noCG 0.821503393
Ath_IME_k7_coords110_110_complete_secondary_noCG 0.821022716
Ath_IME_k7_coords100_250_complete_secondary_noCG 0.819715794
Ath_IME_k4_coords175_175_incomplete_secondary_noCG 0.817962676
Ath_IME_k7_coords180_180_complete_secondary_noCG 0.817873316
Ath_IME_k7_coords90_90_complete_secondary_noCG 0.816591278

These correlation efficients are all still better than the original IMEter v2.0 results.
They suggest that when you remove the effect of CG dinucleotides (due to CpG), you can 
still train an IMEter and that maybe the optimal distance cutoff are not quite as close
as when you include CG dinculeotides in the kmers. When we inspect the kmers present in 
the top parameter file above (Ath_IME_k6_coords125_125_incomplete_secondary_noCG) and then
remove any CG-based kmers, we are left with the following most enriched kmers.

CTAGGG  1.27900001561784
TAGGGT  1.27546359297961
AGGGTT  1.26590072925284
TCTGGG  1.22092351236576
CTGGGT  1.20692200043185
ACCCTA  1.13320279230771
GGGTTT  1.08994571927468
TTAGGG  1.06475602235191


So maybe Tali's motif is the thing to be hunting after all (along with the TAGGG motif).
Back to my old script to find which introns (from Phytozome v9.0) have most copies of 
this Motif (using log-odds threshold 3)

	where_is_motif.pl -motif ~/Work/IME_data/Motifs/tali_motif1.xms -target Athaliana_IME_intron.fa  -bg ati -threshold 3 -mcount | grep -v ": [0-9]$"

This gave a list of 143 introns which have at least 10 copies of Tali's motif. Then tried
finding those with at least 30% of the sequence being occupied by the motif:

	where_is_motif.pl -motif ~/Work/IME_data/Motifs/tali_motif1.xms -target Athaliana_IME_intron.fa  -bg ati -threshold 3 -mdensity | grep -v "0.000%" | grep "[3-9][0-9]\.[0-9][0-9][0-9]%" > tali_motif_high_percent.txt

This gave me 423 introns. Now to find the overlap between them:

	$ cut -f 1 -d " "  tali_motif_top_counts.txt  > to_match
	$ grep -f to_match tali_motif_high_percent.txt
	>IME_Athaliana_30348 motif_density: 48/88 54.545%
	>IME_Athaliana_32355 motif_density: 104/334 31.138%
	>IME_Athaliana_64030 motif_density: 104/334 31.138%
	>IME_Athaliana_105386 motif_density: 108/345 31.304%
	>IME_Athaliana_105387 motif_density: 108/348 31.034%
	>IME_Athaliana_161877 motif_density: 96/291 32.990%

Now work out what genes these introns are from:

	grep -f to_match tali_motif_high_percent.txt | sed 's/ .*//' > tali_motif_best_candiates.txt
	grep -f tali_motif_best_candiates.txt Athaliana_IME_intron.fa
	>IME_Athaliana_30348 TYPE=intron POS=1/9 COORDS=445-532 ID1=19644782 ID2=AT4G16143.1 ISOFORM=primary STRUCTURE=5-3 ORTHOLOG=LOC_Os05g06350.1
	>IME_Athaliana_32355 TYPE=intron POS=1/1 COORDS=93-426 ID1=19645209 ID2=AT4G05050.1 ISOFORM=primary STRUCTURE=5-3 ORTHOLOG=LOC_Os06g46770.1
	>IME_Athaliana_64030 TYPE=intron POS=1/3 COORDS=59-391 ID1=19651751 ID2=AT1G07940.2 ISOFORM=secondary STRUCTURE=5-3 ORTHOLOG=LOC_Os03g08050.1
	>IME_Athaliana_105386 TYPE=intron POS=1/1 COORDS=56-400 ID1=19660213 ID2=AT3G29360.1 ISOFORM=primary STRUCTURE=5-3 ORTHOLOG=LOC_Os03g55070.1
	>IME_Athaliana_105387 TYPE=intron POS=1/1 COORDS=55-402 ID1=19660214 ID2=AT3G29360.2 ISOFORM=secondary STRUCTURE=5-3 ORTHOLOG=LOC_Os03g55070.1
	>IME_Athaliana_161877 TYPE=intron POS=1/7 COORDS=249-539 ID1=19672283 ID2=AT5G36000.1 ISOFORM=primary STRUCTURE=- ORTHOLOG=LOC_Os01g11990.1
	
	grep -f tali_motif_best_candiates.txt tali_motif_top_counts.txt
	>IME_Athaliana_32355 motif_count: 10
	>IME_Athaliana_64030 motif_count: 10
	>IME_Athaliana_105386 motif_count: 11
	>IME_Athaliana_105387 motif_count: 11
	>IME_Athaliana_161877 motif_count: 10

Some spurious matches (as a result of using grep) and two variants of the same gene, so a manual tidy up leaves us:	

	>IME_Athaliana_32355 motif_count: 10 motif_density: 104/334 31.138% POS=1/1 COORDS=93-426 ID1=19645209 ID2=AT4G05050.1
	UBQ11 - Ubiquitin
	
	>IME_Athaliana_64030 motif_count: 10 motif_density: 104/334 31.138%  POS=1/3 COORDS=59-391 ID1=19651751 ID2=AT1G07940.2
	GTP binding Elongation factor Tu family protein
	
	>IME_Athaliana_105386 motif_count: 11 motif_density: 108/345 31.304% POS=1/1 COORDS=56-400 ID1=19660213 ID2=AT3G29360.1
	UGD2 - Encodes one of four UDP-glucose dehydrogenase UGD) genes. 
	
	>IME_Athaliana_161877 motif_count: 10 motif_density: 96/291 32.990% POS=1/7 COORDS=249-539 ID1=19672283 ID2=AT5G36000.1
	Unknown


## Motif vs CG ##
Presence of CG dinucleotides alone gives a very high r^2 value for the test set of introns

Intron	Original IMEter	Expression increase	%CG
dbIMEintron1	15.73	4.3	0.92
dbIMEintron12	53.95	12.5	4.62
dbIMEintron13	4.96	1.4	1.01
dbIMEintron14	29.55	12.3	2.91
dbIMEintron15	6.75	2.2	0.33
dbIMEintron33	2.99	1.1	0.44
dbIMEintron34	12.09	4	0.88
dbIMEintron35	31.95	7.2	2.4
dbIMEintron36	4.69	1.8	0.8
dbIMEintron37	11.52	0.6	1.75
dbIMEintron38	8.36	1.2	1.05
dbIMEintron46	41.02	10.2	3.57
dbIMEintron47	37.53	4.1	2.33
dbIMEintron48	14.91	4.9	0.57
dbIMEintron76	2.89	5.7	1.02
			0.826425979
