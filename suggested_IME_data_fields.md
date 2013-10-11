
## Mandatory fields ##

Data to be tab-separated

1. Unique ID - e.g. dbIME000001
2. Type: Intron or 5' UTR or CDS exon
3. Length
4. Distance from TSS
5. Position: intron 1, intron 2 etc.
6. Number of introns/UTRs in gene
7. Gene ID (e.g. TAIR ID)
8. Species name
9. Sequence

## Other possible fields ##

10. Confirmed by transcript data? Just a yes/no field which might help us extract 'reliable' introns later on. Many Phytozome species have transcript data.
11. Annotated? We might want to store things like 5' UTR introns which are not currently in the Phytozome annotations but which can be inferred from transcript data.
12. Phytozome release details (e.g. v9.0)
13. Flag to mark whether this is a gene of interest (highly expressed or conserved?)
14. Flag to mark whether gene is involved in alternative splicing
15. Arabidopsis thaliana Gene ID for ortholog (Phytozome provides this information in ancillary file (e.g. Acoerulea_195_annotation_info.txt.gz)


All to end up in an SQLite database. Could also consister one column of 'tags' for various metadata, rather than using separate columns.