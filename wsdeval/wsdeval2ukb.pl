#!/usr/bin/perl

use strict;
use XML::LibXML;

die "$0 src/ALL.data.xml > wsd_evalf_ctx.txt\n" unless @ARGV;
my $fname = shift @ARGV;
# PARSER
my $parser = XML::LibXML->new();
$parser->keep_blanks(0);

my $doc = $parser->parse_file($fname);
my $doc_elem = $doc->getDocumentElement;

foreach my $s_elem ($doc_elem->findnodes('/corpus/text/sentence')) {
	my $s_id = $s_elem->getAttribute('id');
	my @C = ();
	foreach my $i_elem ($s_elem->getChildrenByTagName('instance')) {
		my $id = $i_elem->getAttribute('id');
		my $lemma = $i_elem->getAttribute('lemma');
		my $pos = &xr_pos($i_elem->getAttribute('pos'));
		next unless defined ($id) and defined ($lemma) and defined ($pos);
		push @C, "$lemma#$pos#$id";
	}
	next unless @C;
	#print "$s_id\n";
	print join(" ", @C)."\n";
}


# perl -ne 'next unless /\<instance /; next unless /pos=\"([^\"]+?)\"/; print $1."\n";' src/ALL.data.xml | sort | uniq -c | sort -nr
#    4300 NOUN
#    1652 VERB
#     955 ADJ
#     346 ADV

sub xr_pos {
	my $str = shift;

	my %hpos = ( 'NOUN' => 'n',
				 'VERB' => 'v',
				 'ADJ' => 'a',
				 'ADV' => 'r');
	return $hpos{$str};
}
