#!/usr/bin/perl

# Convert database ids of UKB output to lexicographic codes
use strict;
use FindBin qw($Bin);
use File::Basename;

sub usage {
	my $cmd = basename $0;
	my $msg = <<"USG";
Usage: $0 path_to_wn/index.sense > id2lc.map

USG
	$msg .= join("\n", @_) if @_;
	chomp($msg);
	die $msg."\n";
}

&usage("Too few parameters") unless @ARGV;
my $id2lc = &create_lc2id($ARGV[0]);

while (my ($id, $lcv) = each %{ $id2lc }) {
	print "$id\t".join(" ", sort keys %{ $lcv })."\n";
}

# create a 1:N map from sense id to lexicographic code
sub create_lc2id {

  my $wfname = shift;
  # lexicographic code pos:
  #  %1 -> noun
  #  %2 -> verb
  #  %3 -> adj
  #  %4 -> adv
  #  %5 -> adj
  my @Aklpos = ('n', 'v', 'a', 'r', 'a');
  my $lc2id = {};
  open(my $fh, $wfname) || die "Can't open $wfname:$!\n";
  while (<$fh>) {
	chomp;
	my ($lc, $id) = split(/\s+/, $_);
	die "Error in lex. code, line $.\n" unless $lc =~ /%(\d):/;
	my $lc_pos = $Aklpos[$1 - 1];
	next unless $lc_pos;
	die "$lc redefined! in line $." if defined $lc2id->{$lc};
	$lc2id->{"$id-$lc_pos"}->{$lc} = 1;
  }
  return $lc2id;
}
