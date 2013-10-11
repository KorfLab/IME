#!/usr/bin/perl
use strict;
use warnings 'FATAL' => 'all';
use FAlite;
use IK;
use Getopt::Std;
use vars qw($opt_h $opt_m $opt_M $opt_5 $opt_3 $opt_c $opt_C $opt_f $opt_o $opt_u);
getopts('hm:M:5:3:c:Cf:o:u');
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
  -u UTRs are required
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
my $UTRs_REQ = $opt_u;
my $FASTA_OUT = $opt_f;
my $TABLE_OUT = $opt_o;

my ($DIR) = @ARGV;
my $FASTA = `ls $DIR/assembly/*fa.gz`;                chomp $FASTA;
my $GFF   = `ls $DIR/annotation/*gene_exons.gff3.gz`; chomp $GFF;
die "can't find FASTA" unless -e $FASTA;
die "can't find GFF"   unless -e $GFF;

# part 1: get all the genes (actually only the longest)
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
open(my $ffh, "gunzip -c $FASTA |") or die;
my $fasta = new FAlite($ffh);
while (my $entry = $fasta->nextEntry) {
	my ($chrom) = $entry->def =~ /^>(\S+)/;
	foreach my $id (keys %gene) {
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
	
	if ($UTRs_REQ) {
		if ($n5 == 0 or $n3 == 0) {$gene{$id}{error}{UTRs_missing}++}
	}
	
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

# part 5: the big table
my $skipped;
my $imeID = 0;
foreach my $id (sort keys %gene) {
	if ($gene{$id}{error}) {
		$skipped++;
		next;
	}
	
	my $tss = $gene{$id}{tss};
		
	my @type = qw(exon intron five_prime_UTR three_prime_UTR);
	foreach my $type (@type) {
		my $ftotal = @{$gene{$id}{$type}};
		for (my $i = 0; $i < @{$gene{$id}{$type}}; $i++) {
			my $f = $gene{$id}{$type}[$i];
			print join("\t",
				$id,
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

print STDERR "skipped $skipped genes due to errors\n";
#browse(\%gene);


__END__
# part Y: some summary stats
my $genes = scalar keys %gene;
my ($utr5, $utr3, $complete, $exons, $cds, $errors, $kept, %splices);
foreach my $id (keys %gene) {
	my $u5 = exists $gene{$id}{five_prime_UTR};
	my $u3 = exists $gene{$id}{three_prime_UTR};
	$exons += @{$gene{$id}{exon}};
	$cds += @{$gene{$id}{CDS}};
	$utr5++ if $u5;
	$utr3++ if $u3;
	$complete++ if $u5 and $u3;
	$errors++ if exists $gene{$id}{error};
	foreach my $type (keys %{$gene{$id}{splices}}) {
		$splices{$type} += $gene{$id}{splices}{$type};
	}
}
$kept = @keep;

print "
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
print "Splice sites:\n";
foreach my $type (sort {$splices{$b} <=> $splices{$a}} keys %splices) {
	print "\t$type\t$splices{$type}\n";
}

