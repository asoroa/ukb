#!/usr/bin/perl

# wnindex2ukbdict.pl
# given an index.sense file, generate a file with format:
#
#  word syn1-pos:freq ...
#
# @@ 07/11 also, the order of the synsets is consistent with 3rd
#          column of index.sense file (for every pos)

use strict;
use Getopt::Std;

my %opts;

getopt('', \%opts);

my $opt_v = $opts{'v'};

die "Usage: $0 [-v] index.sense\n\t-v\t\tBe verbose\n
     example: perl $0 path_to_wordnet/dict/index.sense > dict.txt\n" unless scalar @ARGV;

my $idx_sense_file = shift @ARGV;

my %H;

open(my $fh, $idx_sense_file) || die "Can't open $idx_sense_file:$!\n";
while (<$fh>) {
  my ($kl, $offset, $sense_n, $f) = split(/\s+/, $_);

  die "Bad index in line $.\n" unless $sense_n;
  $sense_n--; # start at 0
  my ($word, $pos) = &split_kl($kl);
  my $cid_str = "$offset-$pos";
  $cid_str .= $f ? ":$f" : ":0";
  my $wref = $H{$word};
  if (defined($wref)) {
    warn "Synset redefined: $kl in line $.\n" if $opt_v and defined $wref->{$pos}->[$sense_n];
    $wref->{$pos}->[$sense_n] = $cid_str;
  } else {
    $H{$word}->{$pos}->[$sense_n] = $cid_str;
  }
}


foreach my $hword (sort keys %H) {
  print "$hword";
  my $hk = $H{$hword};
  foreach my $hpos ( keys %{$hk} ) {
    print " ";
    print join (" ", @{ $hk->{$hpos} });
  }
  print "\n";
}

sub split_kl {
  my $kl=shift;

  my @Aklpos = ('n', 'v', 'a', 'r', 'a');

  my ($w, $aux) = split(/\%/, $kl);
  my $pos = (split(/\:/, $aux))[0];
  return ($w, $Aklpos[$pos-1]);
}
