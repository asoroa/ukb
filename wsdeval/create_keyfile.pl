#!/usr/bin/perl

# Convert ukb 'raw' output to an output suitable for WSDEval (basically, by
# converting synset representation from Wnet database id's to lexicagraphic
# codes).

use strict;
use File::Basename;

sub usage {
	my $cmd = basename $0;
	my $msg = <<"USG";
Usage: $0 id2lc.map keyfile.id.txt > keyfile.lc.txt
USG
	$msg .= join("\n", @_) if @_;
	chomp($msg);
	die $msg."\n";
}

&usage("Too few parameters") unless @ARGV == 2;
my $idmap_fname = $ARGV[0];
my $key_fname = $ARGV[1];
my $id2lc = &read_map($idmap_fname);
my @RES;
my $l_n=0;
open(my $fh, $key_fname) or die "Can't open $key_fname:$!\n";
while (<$fh>) {
	$l_n++;
	chomp;
	next if /^\!\!/;
	my @l = split (/\s+/, $_);
	shift @l;					# remove context
	my $instid = shift @l;
	my $lemma = pop @l;
	pop @l;						# remove "!!" comment mark
	my @IDS;
	foreach my $id (@l) {
		my $lcv = $id2lc->{$id};
		next unless defined $lcv;
		my @good_ids = grep { /^$lemma%/ } @{$lcv};
		next unless scalar @good_ids;
		push @IDS, @good_ids;
	}
	next unless scalar @IDS;
	push (@RES, "$instid ". join(" ", @IDS));
}
print join("\n", sort @RES);
print "\n";

sub read_map {
	my $map_fname = shift;
	# read map file

	my $M = {};
	open(my $fh, $map_fname) or die "Can't open $map_fname:$!\n";
	while (<$fh>) {
		chomp;
		my ($id, @lcv) = split(/\s+/, $_);
		$M->{$id} = \@lcv;
	}
	return $M;
}
