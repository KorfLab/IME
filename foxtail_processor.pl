#!/usr/bin/perl

use strict; use warnings;
use FAlite;

die "usage: $0 <fasta> <gff>\n" unless @ARGV == 2;


# part 1: get all the genes (actually only the longest)
my %gene;
my %short;

open(GFF, $ARGV[1]) or die;
while (<GFF>) {
	next if /^#/;
	my ($seq, $so, $fea, $beg, $end, $sco, $str, $frm, $group) = split;
	my ($id) = $group =~ /ID=PAC:(\d+)/;
	if ($fea eq 'mRNA' and $group =~ /longest=0/) {$short{$id} = 1}
	next if not defined $id;
	$gene{$id}{chrom} = $seq;
	$gene{$id}{strand} = $str;
	push @{$gene{$id}{$fea}}, {
		beg => $beg,
		end => $end,
	};
}
close GFF;
my $genes = keys %gene;
print STDERR "found $genes total genes\n";

# keep genes with introns and both UTRs
foreach my $id (keys %gene) {
	my $exons = @{$gene{$id}{exon}};
	my $check5 = exists $gene{$id}{five_prime_UTR};
	my $check3 = exists $gene{$id}{three_prime_UTR};
	delete $gene{$id} unless $exons > 1 and $check5 and $check3;
	delete $gene{$id} if $short{$id};
}

$genes = keys %gene;
print STDERR "kept $genes genes\n";

# part 2: read fasta file and dump introns
open(FASTA, $ARGV[0]) or die;
my $fasta = new FAlite(\*FASTA);
while (my $entry = $fasta->nextEntry) {
	my ($chrom) = $entry->def =~ /^>(\S+)/;
	print STDERR "processing chromosome $chrom";
	foreach my $id (keys %gene) {
		next unless $gene{$id}{chrom} eq $chrom;
		print STDERR ".";
	}
	print STDERR "\n";
}
close FASTA;




__END__
foreach my $id (sort keys %gene) {
	my $utr5 = exists $gene{$id}{five_prime_UTR} ? 1 : 0;
	my $utr3 = exists $gene{$id}{three_prime_UTR} ? 1 : 0;
	print "gene $id, chrom $gene{$id}{chrom}, strand $gene{$id}{strand} $utr5 $utr3\n";	
#	foreach my $exon (@{$gene{$id}{exon}}) {
#		print "\t", $exon->{beg}, " ", $exon->{end}, "\n";
#	}
}
