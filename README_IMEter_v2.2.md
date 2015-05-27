# IMEter v2.2

Work in progress in improving IMEter v2.1 and leveraging better quality annotated data from Phytozome v10.



## Getting data from Phytozome v10.0.4

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
CGCCG   1.22953648512629
CGATC   1.05461245354295
CGATT   1.03625398431636
CGGCG   0.982759636835548
TCGAT   0.964288118654301
TCCGA   0.892261809503596
TCGCG   0.890909071438347
GATCG   0.872558511203674
TAGGG   0.867923670376819
GGGTT   0.861602930056838
```

Now try -q = 50 (without -x):

```bash
# build: ./ime_trainer.pl -k 5 -p 400 -d 400 -c -q 50 Athaliana_IME_intron.fa
# introns counted:     138144
# bases clipped:       0
# proximal introns:    15353
# distal introns:      122791
CGCCG   1.06322721405146
GGGTT   0.977503611146258
TAGGG   0.952333887461336
AGGGT   0.927545796090563
GATCG   0.9193202969594
CGATT   0.914550077442204
CGGCG   0.883751912446206
CGACG   0.858184508444106
CGATC   0.837855159370795
CTGGG   0.831612794751938
```

This changes some of the top kmers that appear (and we see the non-CG-kmers: GGGTT, TAGGG, and AGGGT). Now we add the -x option into the mix (without -q):

```bash
# build: ./ime_trainer.pl -k 5 -p 400 -d 400 -c -x Athaliana_IME_intron.fa
# introns counted:     138144
# bases clipped:       1698391
# proximal introns:    15353
# distal introns:      122791
CGCCG   1.60524841755586
CGATC   1.59753071910004
CGATT   1.49984536467771
CGGCG   1.42480019742612
TCGAT   1.41630270342777
TCGCG   1.37613923833504
TCCGA   1.35693456044754
CGCGA   1.33072771219721
CCGAT   1.29872260925767
CTCCG   1.2769067957431

```bash

Addition of -x option increases log-odds scores quite a lot and changes the top 10 kmers again. So how does -q 50 and -x look:

```bash
# build: ./ime_trainer.pl -k 5 -p 400 -d 400 -c -q 50 -x Athaliana_IME_intron.fa
# introns counted:     138144
# bases clipped:       1698391
# proximal introns:    15353
# distal introns:      122791
CGATT   1.55043311963481
CGCCG   1.52927608597663
CGATC   1.52721337339612
GATCG   1.52563512967403
TCGCG   1.48513682662658
CGGCG   1.48287996440789
CGCGA   1.48097976431641
TCGAT   1.42604652272866
GGGTT   1.38024903705267
CTGGG   1.33097685014413

```bash

Try to see how -q and -x look along side the default results

```bash
Position    Kmer    Default score       With -q=50 and -x
#1          CGCCG   1.22953648512629    1.52927608597663 #2
#2          CGATC   1.05461245354295    1.52721337339612 #3
#3          CGATT   1.03625398431636    1.55043311963481 #1
#4          CGGCG   0.982759636835548   1.48287996440789 #6
#5          TCGAT   0.964288118654301   1.42604652272866 #8
#6          TCCGA   0.892261809503596
#7          TCGCG   0.890909071438347   1.48513682662658 #5
#8          GATCG   0.872558511203674   1.52563512967403 #4
#9          TAGGG   0.867923670376819
#10         GGGTT   0.861602930056838   1.38024903705267 #9
            CGCGA                       1.48097976431641 #7
            CTGGG                       1.33097685014413 #10

```

So changing the training program in this way generally boosts all log odds scores, but in relative terms it reduces the extra boost that was given to CGCCG. It also makes CGATT the highest kmer.

Now redo all of this with -q = 100. Will skip to the headline comparison:

```bash
Position    Kmer    Default score       With -q=100 and -x
#1          CGCCG   1.22953648512629    1.60986472320887 #1
#2          CGATC   1.05461245354295    
#3          CGATT   1.03625398431636    1.46511447708335 #7
#4          CGGCG   0.982759636835548   1.42560435903824 #10
#5          TCGAT   0.964288118654301   
#6          TCCGA   0.892261809503596
#7          TCGCG   0.890909071438347   1.4331938769954  #8
#8          GATCG   0.872558511203674   1.54163035481501 #2
#9          TAGGG   0.867923670376819   1.50126320812673 #4
#10         GGGTT   0.861602930056838   1.42600510249279 #9
            CGCGA                       1.53087160432573 #3
            CGACG                       1.48728417779348 #5
            AGGGT                       1.47154824508    #6
```

So with -q at 100 bp, this makes more drastic changes compared to the old default values. What about -q=150?

```bash
Position    Kmer    Default score       With -q=150 and -x
#1          CGCCG   1.22953648512629    
#2          CGATC   1.05461245354295    1.35676232597733 #9
#3          CGATT   1.03625398431636    1.35395792474312 #10
#4          CGGCG   0.982759636835548   1.44595258982649 #4
#5          TCGAT   0.964288118654301   
#6          TCCGA   0.892261809503596
#7          TCGCG   0.890909071438347   1.38560982568778 #7
#8          GATCG   0.872558511203674   1.399509085695   #5
#9          TAGGG   0.867923670376819   1.46102377834578 #3
#10         GGGTT   0.861602930056838   
            CGCGA                       1.61491451043138 #1
            CGACG                       1.57132708389913 #2
            TTCGA                       1.39159401248468 #6
            TCGAA                       1.38258112063352 #8
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
Position    Kmer    Default score       With -q=125
#1          CGCCG   1.22953648512629    1.06633217188795 #1
#2          CGATC   1.05461245354295    0.876360044391628 #3
#3          CGATT   1.03625398431636    0.883135549901786 #2
#4          CGGCG   0.982759636835548   0.791173522227393 #8
#5          TCGAT   0.964288118654301   0.805111240054373 #6
#6          TCCGA   0.892261809503596   0.738486076412956 #12
#7          TCGCG   0.890909071438347   0.859387867265468 #4
#8          GATCG   0.872558511203674   0.719652200853388 #15
#9          TAGGG   0.867923670376819   0.794060471074419 #7
#10         GGGTT   0.861602930056838   0.819760174260988 #5
#11         TTCGA   0.853614225351389   0.747575198743333 #10
#12         CGCGA   0.847351833937376   0.784740501033 #9
#13         CGTCG   0.830087669049207
#14         CGACG   0.829707674657411
#15         AGGGT   0.822745275109967   0.744162866759895 #11
#16         ATCGA   0.811812039338427   0.713130865185136 #16
#17         CGAAT   0.810244108724747   0.729560512894274 #13
#18         CCGAT   0.806683868159494   0.675179865369885 #18
#19         TCGAA   0.789184771054985   0.726436511010194 #14
#20         AATCG   0.787173377488944   0.70467421301563  #17
            CTCCG                       0.665565846584791 #19
            CTGGG                       0.66020922614789  #20
```

But I find the enriched kmers iun this set to be so similar to the original set, and CGCCG still scores so much higher than other kmers.


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
