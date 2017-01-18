#!/usr/bin/perl

# Create ILI links between gNet and Wnet 3.0.

use strict;
use File::Basename;
use XML::LibXML;

binmode STDOUT, ":utf8";

die "usage: $0 GN_V110_XML\n" unless @ARGV;
my $iDir = shift @ARGV;
die "$iDir is not a dir\n" unless -d $iDir;

my $parser = XML::LibXML->new();
$parser->keep_blanks(0);
my $H = &read_lu2s($iDir);

my $fname = "$iDir/interLingualIndex_DE-EN.xml";

&dump_xlinks($fname, $H);

sub read_lu2s {
	my $iDir = shift;

	my %hpos = ('adj' => 'a',
				'nomen' => 'n',
				'verben' => 'v');

	my $H = {}; # { lu => sid-pos }
	opendir my $dh, $iDir or die "Can't open $iDir:$!\n";
	foreach my $fname (grep { !/^\.\.?$/ } readdir $dh) {
		my ($pos, undef, $ext) = split(/\./, $fname);
		next unless $ext eq 'xml';
		next unless $hpos{$pos};
		&populate_lu2s("$iDir/$fname", $hpos{$pos}, $H);
	}
	closedir $dh;
	return $H;
}

sub populate_lu2s {
	my ($fname, $postag, $H) = @_;

	my $doc = $parser->parse_file($fname);
	my $doc_elem = $doc->getDocumentElement;

	foreach my $synset_elem ($doc_elem->findnodes('/synsets/synset')) {
		my $sid = $synset_elem->getAttribute('id');
		foreach my $lu_elem ($synset_elem->getChildrenByTagName('lexUnit')) {
			my $luid = $lu_elem->getAttribute('id');
			$H->{$luid} = "$sid-$postag";
		}
	}
}


sub dump_xlinks {
	my ($fname, $H) = @_;
	my $doc = $parser->parse_file($fname);
	my $doc_elem = $doc->getDocumentElement;

	foreach my $ili_elem ($doc_elem->findnodes('/interLingualIndex/iliRecord')) {
		my $lid = $ili_elem->getAttribute('lexUnitId');
		my $from = $H->{$lid};
		next unless defined $from;
		my $to = $ili_elem->getAttribute('pwn30Id');
		next unless defined $to;
		$to =~ s/^ENG30/eng-30/;
		print "u:$from v:$to s:gnet2wn30\n";
	}
}
