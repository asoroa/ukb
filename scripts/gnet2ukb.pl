#!/usr/bin/perl

# Create ukb files (dictionary, relations) from GermaNet. Needs one
# parameter, the directory where GermaNet lies.
#
# notes:
#
#   - the script creates two files: "gnet.lex" (dictionary) and
#   "gnet_rels.txt" (relations file). Use -D and -R switches to change
#   those.
#
#   - requires XML::LibXML library
#
#   - in the dictionary, spaces in multiword expressions are replaced with
#     underscores.
#
#   - By default, all headwords are lowercased in the dictionary. Use switch
#   'c' to preserve case in headwords.

use strict;
use File::Basename;
use Getopt::Std;
use XML::LibXML;

sub usage {

	my $cmd = basename $0;
	my $str = <<"USG";
Usage: $cmd [-hcd] [-s src_name] [-D dict_file] [-R relations file] GNet_XML_DIR
\t-h\t\t\t\t help
\t-c\t\t\t\t Preserve case in headwords
\t-d\t\t\t\t Do not create an undirected graph (respect relation type)
\t-D dict_file\t Name of the dictionary file (default 'gnet.lex')
\t-R rels_file\t Name of the relations file (default 'gnet_rels.txt')
\t-s source_name\t Name of the source (by default, 'gnet')

Example: perl $cmd GN_V110_XML\n
USG
	$str .= join("\n", @_) if @_;
	chomp($str);
	die $str."\n";
}

my %opts;
getopts('hcdD:R:s:', \%opts);
&usage() if $opts{'h'};
&usage("Error: Missing directory") unless @ARGV;

my $iDir = $ARGV[0];
&usage("Error: $iDir must be a directory") unless -d $iDir;
my $gnrfname = "$iDir/gn_relations.xml";
&usage("Error: $gnrfname not found\n") unless -f $gnrfname;

my $dfname = "gnet.lex";
$dfname = $opts{'D'} if defined $opts{'D'};
my $rfname = "gnet_rels.txt";
$rfname = $opts{'R'} if defined $opts{'R'};
my $srcname = "gnet";
$srcname = $opts{'s'} if defined $opts{'s'};

my $parser = XML::LibXML->new();
$parser->keep_blanks(0);

# dict

my ($D, $DP) = &read_dict($iDir);
open(my $dfo, ">$dfname") or die "$dfname: $!\n";
binmode $dfo, ":utf8";
if ($opts{'c'}) {
	&dump_dict_tc($D, $dfo);
} else {
	&dump_dict($D, $dfo);
}

# rels
open(my $rfo, ">$rfname") or die "$rfname: $!\n";
binmode $rfo, ":utf8";
&gn_rels($gnrfname, $DP, $rfo);

sub gn_rels {
	my ($fname, $DP, $fo) = @_;

	my $doc = $parser->parse_file($fname);
	my $doc_elem = $doc->getDocumentElement;

	foreach my $cr_elem ($doc_elem->findnodes('/relations/con_rel')) {
		my $f = $cr_elem->getAttribute('from');
		my $pf = $DP->{$f};
		next unless $pf;
		my $from = $f."-".$pf;
		my $t = $cr_elem->getAttribute('to');
		my $pt = $DP->{$t};
		next unless $pt;
		my $to = $t."-".$pt;
		my $name = $cr_elem->getAttribute('name');
		next unless $name;
		my $dir = $cr_elem->getAttribute('dir');
		my $d_src = "d:0";
		$d_src = "d:1" if defined $opts{'d'} and $dir ne "both";
		print $fo "u:$from v:$to t:$name $d_src";
		print $fo " s:".$srcname if $srcname;
		print $fo "\n";
	}
}

sub read_dict {
	my $iDir = shift;

	my %hpos = ('adj' => 'a',
				'nomen' => 'n',
				'verben' => 'v');

	my $D = {}; # { hw => {hw_TC => { pos } => [sid, sid]} }
	my $DP = {}; # { sid => pos }
	opendir my $dh, $iDir or die "Can't open $iDir:$!\n";
	foreach my $fname (grep { !/^\.\.?$/ } readdir $dh) {
		my ($pos, undef, $ext) = split(/\./, $fname);
		next unless $ext eq 'xml';
		next unless $hpos{$pos};
		&populate_dict("$iDir/$fname", $hpos{$pos}, $D, $DP);
	}
	closedir $dh;
	return ($D, $DP);
}

sub dump_dict {

	my ($D, $fo) = @_;
	foreach my $hw (sort keys %{ $D }) {
		# write order: 'n', 'v', 'a'
		my $entry = $D->{$hw};
		my @S;
		foreach my $hw_tc (sort keys %{ $entry }) {
			my $entry_tc = $entry->{$hw_tc};
			foreach my $postag ('n', 'v', 'a') {
				next unless $entry_tc->{$postag};
				foreach my $sid (@{ $entry_tc->{$postag} }) {
					next unless defined $sid;
					push @S, $sid;
				}
			}
		}
		next unless @S;
		my @SS = &remove_dupli(@S);
		print $fo $hw." ".join(" ", @SS)."\n";
	}
}

sub dump_dict_tc {

	my ($D, $fo) = @_;
	foreach my $hw (sort keys %{ $D }) {
		# write order: 'n', 'v', 'a'
		my $entry = $D->{$hw};
		foreach my $hw_tc (sort keys %{ $entry }) {
			my $entry_tc = $entry->{$hw_tc};
			my @S;
			foreach my $postag ('n', 'v', 'a') {
				next unless $entry_tc->{$postag};
				foreach my $sid (@{ $entry_tc->{$postag} }) {
					next unless defined $sid;
					push @S, $sid;
				}
			}
			next unless @S;
			print $fo $hw_tc." ".join(" ", @S)."\n";
		}
	}
}

sub populate_dict {
	my ($fname, $postag, $D, $DP) = @_;

	my $doc = $parser->parse_file($fname);
	my $doc_elem = $doc->getDocumentElement;

	foreach my $synset_elem ($doc_elem->findnodes('/synsets/synset')) {
		my $sid = $synset_elem->getAttribute('id');
		die unless defined $sid;
		foreach my $lu_elem ($synset_elem->getChildrenByTagName('lexUnit')) {
			# single orthForm elem inside lexUnit
			my $of_elems = $lu_elem->getChildrenByTagName('orthForm');
			die unless $of_elems->size() == 1;
			my $hw_tc = $of_elems->item(0)->textContent;
			$hw_tc =~ s/\s+/_/go;
			my $hw = lc($hw_tc);
			my $idx = $lu_elem->getAttribute('sense');
			my $entry = $D->{$hw};
			if (not defined $entry) {
				$D->{$hw} = {};
				$entry = $D->{$hw};
			}
			my $entry_tc = $entry->{$hw_tc};
			if (not defined $entry_tc) {
				$entry->{$hw_tc} = {};
				$entry_tc = $entry->{$hw_tc};
			}
			my $rhs = $entry_tc->{$postag};
			if (not defined $rhs) {
				$entry_tc->{$postag} = [];
				$rhs = $entry_tc->{$postag};
			}
			die "Duplicate entry $hw_tc, pos $postag, sense $idx (file $fname)\n" if defined $rhs->[$idx - 1];
			$rhs->[$idx - 1] = "$sid-$postag";
			$DP->{$sid} = $postag;
		}
	}
}

sub remove_dupli {
	my @A = ();
	my %h;
	foreach my $s (@_) {
		push @A, $s unless $h{$s};
		$h{$s} = 1;
	}
	return @A;
}
