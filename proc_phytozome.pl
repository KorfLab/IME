#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use IK;
use Getopt::Std;
use vars qw($opt_h $opt_m $opt_M $opt_5 $opt_3 $opt_c $opt_C);
getopts('hm:M:5:3:c:Cu');
use DataBrowser;

my $MIN_INTRON = 35;
my $MAX_INTRON = 9999;
my $MAX_5UTR   = 3;
my $MAX_3UTR   = 2;
my $MIN_CDS    = 180;
my $CDS_MOD3   = 0;

die "
usage: proc_phytozome.pl [options] <parent directory>

gene-level errors (will be omitted from output)
  -m <int> minimum intron size [$MIN_INTRON]
  -M <int> maximum intron size [$MAX_INTRON]
  -5 <int> maximum # 5'UTR exons [$MAX_5UTR]
  -3 <int> maximum # 3'UTR exons [$MAX_3UTR]
  -c <int> minimum CDS length [$MIN_CDS]
  -C CDS length must be a multiple of 3
" if @ARGV != 1 or $opt_h;

$MIN_INTRON = $opt_m if $opt_m;
$MAX_INTRON = $opt_M if $opt_M;
$MAX_5UTR   = $opt_5 if $opt_5;
$MAX_3UTR   = $opt_3 if $opt_3;
$MIN_CDS    = $opt_c if $opt_c;
$CDS_MOD3   = $opt_C;

my ($DIR) = @ARGV;

# will process 3 files: 
# 1) a FASTA file with the genome assembly/sequence
# 2) gene_exons GFF3 file with coordinates of mRNAs, exons, and UTRs
# 3) Phytozome's annotation file which has transcript names and A. thaliana ortholog assignments

my $FASTA      = `ls $DIR/assembly/*fa.gz`;                    chomp $FASTA;
my $GFF        = `ls $DIR/annotation/*gene_exons.gff3.gz`;     chomp $GFF;
my $ANNOTATION = `ls $DIR/annotation/*annotation_info.txt.gz`; chomp $ANNOTATION;


die "can't find FASTA"       unless -e $FASTA;
die "can't find GFF"         unless -e $GFF;
die "can't find ANNOTATION"  unless -e $ANNOTATION;


# part 1: get all the genes (actually only the longest variant at each locus)
open(my $gfh, "gunzip -c $GFF |") or die;
my %gene;
while (<$gfh>) {
	next if /^#/;
	my ($seq, $so, $fea, $beg, $end, $sco, $str, $frm, $group) = split;
	my ($id) = $group =~ /ID=PAC:(\d+)/;
	if ($fea eq 'mRNA') {
		my ($name) = $group =~ /Name=(\S+?);/;
		$gene{$id}{name} = $name;
		$gene{$id}{longest} = 1 if $group =~ /longest=1/;
	}
	next if not defined $id;
	$gene{$id}{chrom} = $seq;
	$gene{$id}{strand} = $str;
	push @{$gene{$id}{$fea}}, {
		beg => $beg,
		end => $end,
	};
}
close $gfh;

foreach my $id (keys %gene) {
	foreach my $feature ('exon', 'three_prime_UTR', 'five_prime_UTR') {
		if (not defined $gene{$id}{$feature}) {
			$gene{$id}{$feature} = [];
		} else {
			my @stuff = sort {$a->{beg} <=> $b->{beg}} @{$gene{$id}{$feature}};
			$gene{$id}{$feature} = \@stuff;
		}
	}
	
	$gene{$id}{tss} = $gene{$id}{exon}[0]{beg};
}


# part 2: create introns and cds_length attributes
foreach my $id (keys %gene) {
	my @exon = @{$gene{$id}{exon}};
	my @intron;
	for (my $i = 1; $i < @exon; $i++) {
		push @intron, {
			beg => $exon[$i-1]{end} + 1,
			end => $exon[$i  ]{beg} - 1,
		}
	}
	$gene{$id}{intron} = \@intron;
	
	my $cds_length;
	foreach my $exon (@{$gene{$id}{CDS}}) {
		my $len = $exon->{end} - $exon->{beg} + 1;
		$cds_length += $len;
	}
	
	$gene{$id}{cds_length} = $cds_length;	
}


# part 3: read fasta file to add sequence-based attributes

# use a secondary hash to track genes that have had their sequences extracted
# this will get smaller as each gene is processed and this will speed things up a lot
my %tmp_gene = %gene;

open(my $ffh, "gunzip -c $FASTA |") or die;
my $fasta = new FAlite($ffh);
while (my $entry = $fasta->nextEntry) {
	my ($chrom) = $entry->def =~ /^>(\S+)/;
		
#	foreach my $id (keys %gene) {
	foreach my $id (keys %tmp_gene) {

		next unless $gene{$id}{chrom} eq $chrom;
		
		# add sequence attributes
		my $strand = $gene{$id}{strand};
		foreach my $type ('intron', 'exon', 'five_prime_UTR', 'three_prime_UTR') {
			foreach my $f (@{$gene{$id}{$type}}) {
				$f->{seq} = extract_seq($f, $strand, $entry);
			}
		}
		
		# build cds sequence
		my $cds_seq;
		foreach my $cds (sort {$a->{beg} <=> $b->{beg}} @{$gene{$id}{CDS}}) {
			$cds_seq .= substr($entry->{SEQ}, $cds->{beg}-1,
				$cds->{end} - $cds->{beg} + 1);
		}
		if ($gene{$id}{strand} eq '-') {
			$cds_seq =~ tr[ACGTRYMKWSBDHV]
			              [TGCAYRKMWSVHDB];
			$cds_seq = reverse $cds_seq;
		}
		$gene{$id}{cds_seq} = $cds_seq;
		
		# translate
		$gene{$id}{protein} = IK::translate($cds_seq);
		
		# can remove the gene ID from tmp_gene as we won't need to look at it again
		delete $tmp_gene{$id};	
	}
}
close $ffh;

sub extract_seq {
	my ($feature, $strand, $entry) = @_;
	
	my $seq = substr($entry->{SEQ}, $feature->{beg} -1,
		$feature->{end} - $feature->{beg} + 1);
	if ($strand eq '-') {
		$seq =~ tr[ACGTRYMKWSBDHV]
		          [TGCAYRKMWSVHDB];
		$seq = reverse $seq;
	}
	
	return $seq;
}


# part 4: identify error conditions

my %canonical = (
	"GT..AG" => 1,
	"GC..AG" => 1,
	"AT..AC" => 1,
);

foreach my $id (keys %gene) {
	
	# intron length
	foreach my $intron (@{$gene{$id}{intron}}) {
		my $len = $intron->{end} - $intron->{beg} +1;
		if    ($len < $MIN_INTRON) {$gene{$id}{error}{min_intron}++}
		elsif ($len > $MAX_INTRON) {$gene{$id}{error}{max_intron}++}
	}
	
	# 5' and 3' UTR count
	my $n5 = defined @{$gene{$id}{five_prime_UTR}}
		? @{$gene{$id}{five_prime_UTR}}
		: 0;
	my $n3 = defined @{$gene{$id}{three_prime_UTR}}
		? @{$gene{$id}{three_prime_UTR}}
		: 0;
	$gene{$id}{error}{max_5utr}++ if $n5 > $MAX_5UTR;
	$gene{$id}{error}{max_3utr}++ if $n3 > $MAX_3UTR;
		
	# exon count
	if (not defined @{$gene{$id}{exon}}) {$gene{$id}{error}{no_exons}++}
	
	# cds length
	if ($gene{$id}{cds_length} < $MIN_CDS) {$gene{$id}{error}{min_cds}++}
	if ($gene{$id}{cds_length} % 3 != 0 and $CDS_MOD3)
		{$gene{$id}{error}{cds_mod3}++}
	
	# multiple stop codons
	my $protein = $gene{$id}{protein};
	my $stop_count = $protein =~ tr/*/*/;
	if ($protein =~ /\*$/) {$stop_count--}
	if ($stop_count != 0) {$gene{$id}{error}{internal_stop} = $stop_count}
	
	# incomplete start/stop
	if (substr($protein,  0, 1) ne 'M') {$gene{$id}{error}{start_not_found} = 1}
	if (substr($protein, -1, 1) ne '*') {$gene{$id}{error}{stop_not_found} = 1}
		
	# splice sites
	foreach my $intron (@{$gene{$id}{intron}}) {
		my $don = substr($intron->{seq}, 0, 2);
		my $acc = substr($intron->{seq}, -2);
		my $type = "$don..$acc";
		$gene{$id}{splices}{$type}++;
		if (not defined $canonical{$type}) {
			$gene{$id}{error}{non_canonical_splice}++;
		}
	}

}


# part 5: extract relevant info from annotation file
open(my $afh, "gunzip -c $ANNOTATION |") or die;
while (<$afh>) {
	next if /^#/;
	
  	# extract the 3 fields we need
  	# 1. Phytozome ID
  	# 2. transcript ID/name from source database (e.g. LOC_Os07g46980.1 for rice)
  	# 3. Best hit for A. thaliana ortholog (may not be present)
  	
	my @f = split(/\t/);
	my ($id ,$transcript_ID, $at_ortholog) = ($f[0], $f[2], $f[10]);

	# at_ortholog may not be present, so set to NA
	$at_ortholog = "NA" if not defined $at_ortholog;
	$at_ortholog = "NA" if $at_ortholog eq "";

	# add to master hash
	$gene{$id}{local_id}    = $transcript_ID;
	$gene{$id}{at_ortholog} = $at_ortholog;
}
close $afh;




# part 6: the big table
my ($skipped, $kept);
my $imeID = 0;
foreach my $id (sort keys %gene) {
	if ($gene{$id}{error}) {
		$skipped++;
		next;
	} else {$kept++}
	
	my $tss = $gene{$id}{tss};
		
	my @type = qw(exon intron five_prime_UTR three_prime_UTR);
	foreach my $type (@type) {
		my $ftotal = @{$gene{$id}{$type}};
		for (my $i = 0; $i < @{$gene{$id}{$type}}; $i++) {
			my $f = $gene{$id}{$type}[$i];

			# need to print whether this isoform is the primary version (longest)
			# otherwise, call it secondary. Might use this info later on in downstream steps
			my $isoform = "secondary";
			$isoform = "primary" if $gene{$id}{longest};
			
			# also ned to flag whether this feature belongs to a complete transcript (with 5' & 3' UTR)
			# or is missing one or the other
			my $transcript_status;
			my $n5 = @{$gene{$id}{five_prime_UTR}};
			my $n3 = @{$gene{$id}{three_prime_UTR}};
			warn "n5 = $n5\n";
			if ($n5 and $n3){
				$transcript_status = "5-3";								
			} elsif($n5){
				$transcript_status = "5-";					
			} elsif($n3){
				$transcript_status = "-3";		
			} else{
				$transcript_status = "-";			
			}
			
			print join("\t",
				$id,
				$gene{$id}{local_id},
				$isoform,
				$transcript_status,
				$gene{$id}{at_ortholog},
				$type,
				$f->{beg} +1 - $tss,
				$f->{end} +1 - $tss,
				$i + 1,
				$ftotal,
				$f->{seq},
				
			), "\n";
		}
	}
}

# part 7: some summary stats
my $genes = scalar keys %gene;
my %error;
my ($utr5, $utr3, $complete, $exons, $cds, $errors, %splices);
foreach my $id (keys %gene) {
	my $u5 = exists $gene{$id}{five_prime_UTR};
	my $u3 = exists $gene{$id}{three_prime_UTR};
	$exons += @{$gene{$id}{exon}};
	$cds += @{$gene{$id}{CDS}};
	$utr5++ if $u5;
	$utr3++ if $u3;
	$complete++ if $u5 and $u3;
	
	if (exists $gene{$id}{error}) {
		$errors++; 

		foreach my $etype (keys %{$gene{$id}{error}}) {
			$error{$etype}++;
		}
	}
	foreach my $type (keys %{$gene{$id}{splices}}) {
		$splices{$type} += $gene{$id}{splices}{$type};
	}
}


print STDERR "
GENOME: $DIR
GENES: $genes
EXONS: $exons
CDS:   $cds
5'UTR: $utr5
3'UTR: $utr3
COMPLETE: $complete
ERRORS: $errors
KEPT: $kept
";
print STDERR "Error types:\n";
foreach my $etype (sort keys %error) {
	print STDERR "\t$etype: $error{$etype}\n";
}
print STDERR "Splice sites:\n";
foreach my $type (sort {$splices{$b} <=> $splices{$a}} keys %splices) {
	print STDERR "\t$type\t$splices{$type}\n";
}



