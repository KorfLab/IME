#!/usr/bin/perl
# Copyright (C) 2007-2014 Ian Korf, Genís Parra, and Keith Brandam
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

use strict;
use warnings 'FATAL' => 'all';
use File::Basename;
use FAlite;
use IK;
use Getopt::Std;
use vars qw($opt_h $opt_m $opt_M $opt_5 $opt_3 $opt_c $opt_C $opt_s $opt_v $opt_i $opt_1);
getopts('hm:M:5:3:c:Cs:vi:1');
use DataBrowser;

my $MIN_INTRON = 35;
my $MAX_INTRON = 9999;
my $MAX_5UTR   = 3;
my $MAX_3UTR   = 2;
my $MIN_CDS    = 180;
my $CDS_MOD3   = 0;
my $VERBOSE;
my $STUB;

die "
usage: proc_phytozome.pl [options] <parent Phytozome directory>

gene-level errors (these will be omitted from output files)
  -m <int> minimum intron size [$MIN_INTRON]
  -M <int> maximum intron size [$MAX_INTRON]
  -5 <int> maximum # 5'UTR exons [$MAX_5UTR]
  -3 <int> maximum # 3'UTR exons [$MAX_3UTR]
  -c <int> minimum CDS length [$MIN_CDS]
  -C CDS length must be a multiple of 3
  -s <text> stub prefix name for inclusion in FASTA header
  -1 Specify a species-specific directory (default is to work through all species)
  -v verbose mode. Shows progress through each stage.
  -i ignore species that match pattern \"x\" (separate multiple species by spaces)

Minimum usage example: proc_phytozome /path/to/Phytozome/v9.0/
  
" if @ARGV != 1 or $opt_h;

$MIN_INTRON = $opt_m if $opt_m;
$MAX_INTRON = $opt_M if $opt_M;
$MAX_5UTR   = $opt_5 if $opt_5;
$MAX_3UTR   = $opt_3 if $opt_3;
$MIN_CDS    = $opt_c if $opt_c;
$STUB       = $opt_s if $opt_s;
$VERBOSE    = $opt_v if $opt_v;
$CDS_MOD3   = $opt_C;

my ($DIR) = @ARGV;


my @to_ignore;
if ($opt_i){
	@to_ignore = split(/\s+/, $opt_i);
}


# are we processing all species, or just one?
my @species_directories;

if($opt_1){
	@species_directories = "$DIR";
	$species_directories[0] =~ s/\/$//;
} else{
	@species_directories = glob("$DIR*");
}

SPECIES: foreach my $species_dir (@species_directories){

	my ($species, $directory) = fileparse($species_dir);

	# skip early release data
	next if $species eq "early_release";
	
	warn "Processing data for $species\n"; 
	
	# skip any species to ignore
	foreach my $ignore (@to_ignore){
		if($ignore =~ m/$species/){
			warn "\tthis is is on the ignore list, skipping...\n";
			next SPECIES;		
		}
	}
	
	# skip to next species if files exist?
	if(-e "${species}_IME_exon.fa" and 
	   -e "${species}_IME_intron.fa" and 
	   -e "${species}_IME_five_prime_UTR.fa" and 
	   -e "${species}_IME_three_prime_UTR.fa" and 
	   -e "${species}_IME_errors.txt"){
		warn "\tFiles already exist for this species, skipping...\n";
	   	next;
	}

	$STUB = "IME_$species" if (not $opt_s);


	# will process 3 files for each species: 
	# 1) a FASTA file with the genome assembly/sequence
	# 2) gene_exons GFF3 file with coordinates of mRNAs, exons, and UTRs
	# 3) Phytozome's annotation file which has transcript names and A. thaliana ortholog assignments

	my $FASTA      = `ls $directory/$species/assembly/*fa.gz`;                    chomp $FASTA;
	my $GFF        = `ls $directory/$species/annotation/*gene_exons.gff3.gz`;     chomp $GFF;
	my $ANNOTATION = `ls $directory/$species/annotation/*annotation_info.txt.gz`; chomp $ANNOTATION;

	die "can't find FASTA"       unless -e $FASTA;
	die "can't find GFF"         unless -e $GFF;
	die "can't find ANNOTATION"  unless -e $ANNOTATION;

	process_species($species, $directory, $FASTA, $GFF, $ANNOTATION);
}

exit;


####################################
#
# SUBROUTINES
#
####################################


sub process_species{
	my ($species, $directory, $FASTA, $GFF, $ANNOTATION) = @_;


	# Step 1: get all the transcript details from the GFF file
	warn "\tStep 1: processing GFF file\n" if $VERBOSE;

	my %transcripts;

	open(my $gfh, "gunzip -c $GFF |") or die;
	while (<$gfh>) {
		next if /^#/;
		my ($seq, $so, $fea, $beg, $end, $sco, $str, $frm, $group) = split;
		my ($id) = $group =~ /ID=PAC:(\d+)/;
		if ($fea eq 'mRNA') {
			my ($name) = $group =~ /Name=(\S+?);/;
			$transcripts{$id}{name} = $name;
			$transcripts{$id}{longest} = 1 if $group =~ /longest=1/;
		}
		next if not defined $id;
		$transcripts{$id}{chrom} = $seq;
		$transcripts{$id}{strand} = $str;
		push @{$transcripts{$id}{$fea}}, {
			beg => $beg,
			end => $end,
		};
	}
	close $gfh;

	foreach my $id (keys %transcripts) {
		foreach my $feature ('exon', 'three_prime_UTR', 'five_prime_UTR') {
			if (not defined $transcripts{$id}{$feature}) {
				$transcripts{$id}{$feature} = [];
			} else {
				my @stuff = sort {$a->{beg} <=> $b->{beg}} @{$transcripts{$id}{$feature}};
				$transcripts{$id}{$feature} = \@stuff;
			}
		}
	
		$transcripts{$id}{tss} = $transcripts{$id}{exon}[0]{beg};
	}



	# Step 2: create introns and cds_length attributes
	warn "\tStep 2: creating intron information from exon data\n" if $VERBOSE;

	foreach my $id (keys %transcripts) {
		my @exon = @{$transcripts{$id}{exon}};
		my @intron;
		for (my $i = 1; $i < @exon; $i++) {
			push @intron, {
				beg => $exon[$i-1]{end} + 1,
				end => $exon[$i  ]{beg} - 1,
			}
		}
		$transcripts{$id}{intron} = \@intron;
	
		my $cds_length;
		foreach my $exon (@{$transcripts{$id}{CDS}}) {
			my $len = $exon->{end} - $exon->{beg} + 1;
			$cds_length += $len;
		}
	
		$transcripts{$id}{cds_length} = $cds_length;	
	}



	# Step 3: read fasta file to add sequence-based attributes
	warn "\tStep 3: processing genome sequence file\n" if $VERBOSE;

	# use a secondary hash to track genes that have had their sequences extracted
	# this will get smaller as each gene is processed and this will speed things up a lot
	my %tmp_transcripts = %transcripts;
	open(my $ffh, "gunzip -c $FASTA |") or die;
	my $fasta = new FAlite($ffh);

	while (my $entry = $fasta->nextEntry) {
		my ($chrom) = $entry->def =~ /^>(\S+)/;
		foreach my $id (keys %tmp_transcripts) {
			next unless $transcripts{$id}{chrom} eq $chrom;
		
			# add sequence attributes
			my $strand = $transcripts{$id}{strand};
			foreach my $type ('intron', 'exon', 'five_prime_UTR', 'three_prime_UTR') {
				foreach my $f (@{$transcripts{$id}{$type}}) {
					$f->{seq} = extract_seq($f, $strand, $entry);
				}
			}
		
			# build cds sequence
			my $cds_seq;
			foreach my $cds (sort {$a->{beg} <=> $b->{beg}} @{$transcripts{$id}{CDS}}) {
				$cds_seq .= substr($entry->{SEQ}, $cds->{beg}-1,
					$cds->{end} - $cds->{beg} + 1);
			}
			if ($transcripts{$id}{strand} eq '-') {
				$cds_seq =~ tr[ACGTRYMKWSBDHV]
							  [TGCAYRKMWSVHDB];
				$cds_seq = reverse $cds_seq;
			}
			$transcripts{$id}{cds_seq} = $cds_seq;
		
			# translate
			$transcripts{$id}{protein} = IK::translate($cds_seq);
		
			# can remove the gene ID from tmp_transcripts as we won't need to look at it again
			delete $tmp_transcripts{$id};	
		}
	}
	close $ffh;


	# Step 4: identify error conditions
	warn "\tStep 4: looking for errors in annotations\n" if $VERBOSE;

	my %canonical = (
		"GT..AG" => 1,
		"GC..AG" => 1,
		"AT..AC" => 1,
	);

	foreach my $id (keys %transcripts) {
	
		# intron length
		foreach my $intron (@{$transcripts{$id}{intron}}) {
			my $len = $intron->{end} - $intron->{beg} +1;
			if    ($len < $MIN_INTRON) {$transcripts{$id}{error}{min_intron}++}
			elsif ($len > $MAX_INTRON) {$transcripts{$id}{error}{max_intron}++}
		}
	
		# 5' and 3' UTR count
		my $n5 = @{$transcripts{$id}{five_prime_UTR}}
			? @{$transcripts{$id}{five_prime_UTR}}
			: 0;
		my $n3 = @{$transcripts{$id}{three_prime_UTR}}
			? @{$transcripts{$id}{three_prime_UTR}}
			: 0;

		$transcripts{$id}{error}{max_5utr}++ if $n5 > $MAX_5UTR;
		$transcripts{$id}{error}{max_3utr}++ if $n3 > $MAX_3UTR;
		
		# exon count
		if (not @{$transcripts{$id}{exon}}) {$transcripts{$id}{error}{no_exons}++}
	
		# cds length
		if ($transcripts{$id}{cds_length} < $MIN_CDS) {$transcripts{$id}{error}{min_cds}++}
		if ($transcripts{$id}{cds_length} % 3 != 0 and $CDS_MOD3)
			{$transcripts{$id}{error}{cds_mod3}++}
	
		# multiple stop codons
		my $protein = $transcripts{$id}{protein};
		my $stop_count = $protein =~ tr/*/*/;
		if ($protein =~ /\*$/) {$stop_count--}
		if ($stop_count != 0) {$transcripts{$id}{error}{internal_stop} = $stop_count}
	
		# incomplete start/stop
		if (substr($protein,  0, 1) ne 'M') {$transcripts{$id}{error}{start_not_found} = 1}
		if (substr($protein, -1, 1) ne '*') {$transcripts{$id}{error}{stop_not_found} = 1}
		
		# splice sites
		foreach my $intron (@{$transcripts{$id}{intron}}) {
			my $don = substr($intron->{seq}, 0, 2);
			my $acc = substr($intron->{seq}, -2);
			my $type = "$don..$acc";
			$transcripts{$id}{splices}{$type}++;
			if (not defined $canonical{$type}) {
				$transcripts{$id}{error}{non_canonical_splice}++;
			}
		}
	}




	# Step 5: extract relevant info from annotation file
	warn "\tStep 5: extracting extra information from Phytozome annotation file\n" if $VERBOSE;

	open(my $afh, "gunzip -c $ANNOTATION |") or die;
	while (<$afh>) {
		next if /^#/;
	
		# extract the 3 fields we need
		# 1. Phytozome ID
		# 2. transcript ID/name from source database (e.g. LOC_Os07g46980.1 for rice)
		# 3. Best hit for A. thaliana ortholog (may not be present)
	
		my @f = split(/\t/);
		
		# may have to treat these fields slightly differently for some species
		# because - frustratingly - not all species have the same info in Phytozome
		my ($id ,$transcript_ID, $ortholog);
		
		if($species eq 'Csubellipsoidea_C-169'){
			($id ,$transcript_ID, $ortholog) = ($f[0], $f[1], $f[8]);
		} elsif($species eq 'Mpusilla_CCMP1545'){
			($id ,$transcript_ID, $ortholog) = ($f[0], $f[1], $f[8]);
		} elsif($species eq 'Mpusilla_RCC299'){
			($id ,$transcript_ID, $ortholog) = ($f[0], $f[1], $f[8]);
		} else{
			($id ,$transcript_ID, $ortholog) = ($f[0], $f[2], $f[10]);	
		}

		# add to master hash
		$transcripts{$id}{local_id} = $transcript_ID;
		$transcripts{$id}{ortholog} = $ortholog;
	}
	close $afh;




	# Step 6: the big table
	warn "\tStep 6: final output\n" if $VERBOSE;

	foreach my $type (qw(exon intron five_prime_UTR three_prime_UTR)) {
		# need counter for FASTA header
		my $n = 0;
		
		# prepare output file
		my $output_file = "${species}_IME_$type.fa";
		warn "\t\tCreating $output_file\n" if $VERBOSE;

		open(my $out, ">", $output_file) or die "Cannot write to $output_file";

		foreach my $id (sort keys %transcripts) {
	
			# skip genes with errors
			next if ($transcripts{$id}{error});

			my $tss = $transcripts{$id}{tss};
			my $ftotal = @{$transcripts{$id}{$type}};

			for (my $i = 0; $i < @{$transcripts{$id}{$type}}; $i++) {
				my $f = $transcripts{$id}{$type}[$i];

				# need to print whether this isoform is the primary version (longest)
				# otherwise, call it secondary. Might use this info later on in downstream steps
				my $isoform = "secondary";
				$isoform = "primary" if $transcripts{$id}{longest};
			
				# also ned to flag whether this feature belongs to a complete transcript (with 5' & 3' UTR)
				# or is missing one or the other. Capture 'structure' into a variable.
				my $structure;
				my $n5 = @{$transcripts{$id}{five_prime_UTR}};
				my $n3 = @{$transcripts{$id}{three_prime_UTR}};

				if ($n5 and $n3){
					$structure = "5-3";								
				} elsif($n5){
					$structure = "5-";					
				} elsif($n3){
					$structure = "-3";		
				} else{
					$structure = "-";			
				}

				$n++;
				my $tidied_seq = tidy_seq($f->{seq});
				my $beg = $f->{beg} + 1 - $tss;
				my $end = $f->{end} + 1 - $tss;
				my $position = $i + 1;


				# local ID and ortholog information may not be present in annotation file 
				# or may be missing in annotation file (but present in GFF file)
				$transcripts{$id}{local_id} = "NA" if not defined $transcripts{$id}{local_id};
				$transcripts{$id}{local_id} = "NA" if             $transcripts{$id}{local_id} eq "";
				$transcripts{$id}{ortholog} = "NA" if not defined $transcripts{$id}{ortholog};
				$transcripts{$id}{ortholog} = "NA" if             $transcripts{$id}{ortholog} eq "";

				print $out ">${STUB}_$n TYPE=$type POS=${position}/$ftotal COORDS=$beg-$end ";
				print $out "ID1=$id ";
				print $out "ID2=$transcripts{$id}{local_id} ";
				print $out "ISOFORM=$isoform STRUCTURE=$structure ORTHOLOG=$transcripts{$id}{ortholog}\n";
				print $out "$tidied_seq\n";
			}
		}
		close($out);
	}



	my $error_file = "${species}_IME_errors.txt";
	warn "\t\tCreating $error_file\n" if $VERBOSE;
	open(my $out, ">", $error_file) or die "Cannot write to $error_file";

	# Step 7: some summary stats
	my $transcripts = scalar keys %transcripts;
	my %error;
	my ($utr5, $utr3, $complete, $exons, $cds, $errors, %splices) = (0, 0, 0, 0, 0, 0);
	my $kept = 0;

	foreach my $id (keys %transcripts) {

#		my $u5 = exists $transcripts{$id}{five_prime_UTR};
#		my $u3 = exists $transcripts{$id}{three_prime_UTR};
		my $u5 = @{$transcripts{$id}{five_prime_UTR}};
		my $u3 = @{$transcripts{$id}{three_prime_UTR}};
		
		$exons += @{$transcripts{$id}{exon}};
		$cds += @{$transcripts{$id}{CDS}};
		$utr5++ if $u5;
		$utr3++ if $u3;
		$complete++ if $u5 and $u3;
	
		if (exists $transcripts{$id}{error}) {
			$errors++; 

			foreach my $etype (keys %{$transcripts{$id}{error}}) {
				$error{$etype}++;
			}
		} else{
			$kept++;
		}
		foreach my $type (keys %{$transcripts{$id}{splices}}) {
			$splices{$type} += $transcripts{$id}{splices}{$type};
		}
	}

	print $out "
GENOME: ${directory}$species
TRANSCRIPTS: $transcripts
TRANSCRIPTS_WITH_ERRORS: $errors
TRANSCRIPTS_RETAINED: $kept
EXONS: $exons
CDS:   $cds
TRANSCRIPTS_WITH_5'UTR: $utr5
TRANSCRIPTS_WITH_3'UTR: $utr3
TRANSCRIPTS_WITH_BOTH_UTRS: $complete
";

	print $out "Error types:\n";
	foreach my $etype (sort keys %error) {
		print $out "\t$etype: $error{$etype}\n";
	}
	print $out "Splice sites:\n";
	foreach my $type (sort {$splices{$b} <=> $splices{$a}} keys %splices) {
		print $out "\t$type\t$splices{$type}\n";
	}
	close($out);

}

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




# simple routine for tidying up FASTA file to print 60 characters per line
sub tidy_seq{
    my ($seq) = @_;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;

    my ($output_seq) = "";
    my (@seq2) = "";
    my ($start,$end);

    @seq2 = split(//,$seq);
    my $length = @seq2;
    my ($to_add) = int($length/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    $output_seq .= substr($seq,$start);
    return ($output_seq);
}


__END__



















####################################
#
# SUBROUTINES
#
####################################
