#!/usr/bin/perl

# Create contexts for at least 20 words.
# usage:
#       ./ctx20words.pl Data/s3aw_in.txt > Data/s3aw_ctx20.txt

use strict;
use IO::File;
use Getopt::Std;

my %opts;

getopt('s', \%opts);

my $size = $opts{s} ? $opts{s} : 20;

die "Usage: ./ctx20words.pl -m size file_in\n" unless scalar @ARGV;
die "Error: size has to be at least 5\nUsage: ./ctx20words.pl -m size file_in\n" if scalar $size < 5;


my @Doc;
&read_txt($ARGV[0], \@Doc);


my $i = 0;
my $flipflop = 0;

my $Doc_n = scalar(@Doc);

# Problem: if the number of words is less than $size, this doesn't work

die "Text has less than 20 words!\n" if &w_count(\@Doc) < $size;

my ($pre, $post);

my $last_pre = 0;
my $last_post = 1;
my @curs;

my $ctx_id = 0;

while ($i < $Doc_n) {
  $pre = $i;
  $post = $i + 1;
  my $n = scalar(@{$Doc[$i]});
  $flipflop = 1;
  while ($n < $size) {
    if ($flipflop) {
      # decrement $pre if possible
      if ($pre != 0) {
	$pre--;
	$n+= scalar(@{$Doc[$pre]});
      }
    } else {
      # increment $post if possible
      if ($post != $Doc_n) {
	$n+= scalar(@{$Doc[$post]});
	$post++;
      }
    }
    $flipflop = !$flipflop;
  }
  &create_ctx(\$ctx_id, \@Doc, $i, $pre, $post);
  $i++;
}
#&create_ctx(\$ctx_id, \@Doc, \@curs, $last_pre, $last_post);

sub create_ctx {

  my ($idRef, $D, $cur, $pre, $post) = @_;

  my @res;

  for (my $i = $pre; $i < $post; $i++) {
    if ($i == $cur) {
      # input contexts can have fourth value of 1 or 0. If so, mantain it.
      my @aux;
      foreach my $item (@{$D->[$i]}) {
	my @f = split(/\#/, $item);
	if (scalar(@f) == 3) {
	  push @aux, "$item#1";
	} else {
	  push @aux, $item;
	}
      }
      push @res, join(" ", @aux);
    } else {
      # input contexts can have fourth value of 1 or 0. Anyway, set it to 0
      my @aux;
      foreach my $item (@{$D->[$i]}) {
	my @f = split(/\#/, $item);
	if (scalar(@f) == 3) {
	  push @aux, "$item#0";
	} else {
	  push @aux, "$f[0]#$f[1]#$f[2]#0"
	}
      }
      push @res, join(" ", @aux);
    }
  }

  printf("ctx%04d\n", $$idRef);
  $$idRef++;
  print join(" ", @res);
  print "\n";
}


sub read_txt {

  my ($file_in, $aref) = @_;

  my $fh = new IO::File($file_in);
  die "Can't open $file_in:$!\n" unless defined $fh;

  while (<$fh>) {
    chomp;
    my @S = split(/\s+/, $_);
    push @$aref, \@S;
  }
}

sub w_count {

  my $aref = shift;
  my $n = 0;
  foreach (@{$aref}) {
    $n+=scalar @{$_};
  }
  return $n;
}
