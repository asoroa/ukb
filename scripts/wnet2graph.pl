#!/usr/bin/perl

# convert wordnet data.[verb|noun|adj|adv] files to ukb input format
#

use strict;
use File::Basename;
use Getopt::Std;

my %opts;

getopt('s', \%opts);

my %poschar = ("noun" => 'n',
	       "verb" => 'v',
	       "adj" => 'a',
	       "adv" => 'r');

# pointer_symbols for nouns

my %ptr_sym_noun = (
		    '!' =>  ["antonym", "Antonym"],
		    '@' =>  ["hypernym", "Hypernym"],
		    '@i' => ["instance_hypernym", "Instance Hypernym"],
		    '~'  => ["hyponym", "Hyponym"],
		    '~i' => ["instance_hyponym", "Instance Hyponym"],
		    '#m' => ["member_holonym", "Member holonym"],
		    '#s' => ["substance_holonym", "Substance holonym"],
		    '#p' => ["part_holonym", "Part holonym"],
		    '%m' => ["member_meronym", "Member meronym"],
		    '%s' => ["substance_meronym", "Substance meronym"],
		    '%p' => ["part_meronym", "Part meronym"],
		    '='  => ["attribute", "Attribute"],
		    '+'  => ["derivationally_related_form", "Derivationally related form"],
		    ';c' => ["domain_of_synset_topic", "Domain of synset - TOPIC"],
		    '-c' => ["member_of_this_domain_topic", "Member of this domain - TOPIC"],
		    ';r' => ["domain_of_synset_region", "Domain of synset - REGION"],
		    '-r' => ["member_of_this_domain_region", "Member of this domain - REGION"],
		    ';u' => ["domain_of_synset_usage", "Domain of synset - USAGE"],
		    '-u' => ["member_of_this_domain_usage", "Member of this domain - USAGE"]
		   );

# pointer_symbols for verbs

my %ptr_sym_verb = (
		    '!' =>  ["antonym", "Antonym"],
		    '@' =>  ["hypernym", "Hypernym"],
		    '~' =>  ["hyponym", "Hyponym"],
		    '*' =>  ["entailment", "Entailment"],
		    '>' =>  ["cause", "Cause"],
		    '^' =>  ["also_see", "Also see"],
		    '$' =>  ["verb_group", "Verb Group"],
		    '+' =>  ["derivationally_related_form", "Derivationally related form"],
		    ';c' => ["domain_of_synset_topic", "Domain of synset - TOPIC"],
		    ';r' => ["domain_of_synset_region", "Domain of synset - REGION"],
		    ';u' => ["domain_of_synset_usage", "Domain of synset - USAGE"]
		   );

# pointer_symbols for adjectives
my %ptr_sym_adj = (
		   '!' =>  ["antonym", "Antonym"],
		   '&' =>  ["similar_to", "Similar to"],
		   '<' =>  ["participle_of_verb", "Participle of verb"],
		   '\\' => ["pertainym", "Pertainym (pertains to noun)"],
		   '=' =>  ["attribute", "Attribute"],
		   '^' =>  ["also_see", "Also see"],
		   ';c' => ["domain_of_synset_topic", "Domain of synset - TOPIC"],
		   ';r' => ["domain_of_synset_region", "Domain of synset - REGION"],
		   ';u' => ["domain_of_synset_usage", "Domain of synset - USAGE"]
		  );
# The pointer_symbols for adverbs are:

my %ptr_sym_adv = (
		   '!' =>  ["antonym", "Antonym"],
		   '\\' => ["derived_from_adjective", "Derived from adjective"],
		   ';c' => ["domain_of_synset_topic", "Domain of synset - TOPIC"],
		   ';r' => ["domain_of_synset_region", "Domain of synset - REGION"],
		   ';u' => ["domain_of_synset_usage", "Domain of synset - USAGE"]
		  );

&usage() unless @ARGV;

my $opt_d = $opts{'d'};
my $opt_v = $opts{'v'};
my $opt_srcname = $opts{'s'};

foreach my $fullname (@ARGV) {
  my $fname = basename($fullname);
  next unless $fname =~ /^data\.([^\.]+)$/;
  my $fpos = $1;
  my $pos = $poschar{$fpos};
  die $fullname." has unknown pos.\n" unless $pos;
  open(my $fh, $fullname) or die "Can't open $fullname:$!\n";
  print STDERR "Processing $fullname " if $opt_v;
  my $l_n = 0;
  while (my $line = <$fh>) {
    next if $line =~ /^\s/;
    chomp($line);
    &parse_line($line, $pos);
    $l_n++;
    print STDERR " $l_n" if ($opt_v and !($l_n % 100));
  }
  print STDERR "\n" if $opt_v;
}


# synset_offset lex_filenum ss_type w_cnt word lex_id [word lex_id...] p_cnt [ptr...] [frames...] | gloss

sub parse_line {

  my ($str, $src_pos) = @_;

  my ($line) = split(/\|/, $str); # skip the glosses

  my ($src_offset, undef, undef, $w_cnt_hex, @rest) = split(/\s+/, $line);
  #my $w_cnt = sprintf("%x", $w_cnt_hex);
  my $w_cnt = hex($w_cnt_hex);
  die "$.\n" unless $w_cnt;
  while ($w_cnt) {
    shift @rest;
    shift @rest;
    $w_cnt--;
  }

  my $src_cid = "$src_offset-$src_pos";
  my $p_cnt = shift @rest;

  return if $p_cnt == 0;

  while ($p_cnt) {
    my $ptr_symbol = shift @rest;
    my $tgt_offset = shift @rest;
    my $tgt_pos = shift @rest;
    shift @rest;
    my $rel = &rel_sym_name($ptr_symbol, $src_pos);
    print "u:$src_cid v:$tgt_offset-$tgt_pos s:$rel";
    print " s:$opt_srcname" if $opt_srcname;
    print " d:1" if $opt_d;
    print "\n";
    $p_cnt--;
  }

}

sub rel_sym_name {

  my ($ptr_symbol, $pos) = @_;
  my $tr;

  if ($pos eq 'n') {
    $tr = $ptr_sym_noun{$ptr_symbol}->[0];
  } elsif ($pos eq 'v') {
    $tr = $ptr_sym_verb{$ptr_symbol}->[0];
  } elsif ($pos eq 'a') {
    $tr = $ptr_sym_adj{$ptr_symbol}->[0];
  } elsif ($pos eq 'r') {
    $tr = $ptr_sym_adv{$ptr_symbol}->[0];
  }
  $tr = "unknown" unless $tr;
  return $tr;
}

sub usage {

die<<"USG"
Usage: $0 [-s src_name] [-d] [-v] wn_data_file1 ...
\t-s source_name\tName of the source
\t-d\t\tCreate a directed graph
\t-v\t\tBe verbose

Example: perl $0 path_to_wordnet/dict/* > graph_ukb.txt\n
USG

}
