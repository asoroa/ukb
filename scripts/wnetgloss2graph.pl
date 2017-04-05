#!/usr/bin/perl

# wnetgloss2graph.pl
# generate gloss relations

use strict;
use File::Basename;
use XML::LibXML;
use Getopt::Std;

my %opts;

getopt('vds:', \%opts);

my $opt_v = $opts{'v'};
my $dir_str = $opts{'d'} ? " d:1" : "";
my $src_str = $opts{'s'} ? " s:".$opts{'s'} : "";

my $suff_str = $dir_str.$src_str;

sub usage {
	my $cmd = basename $0;
	my $msg = <<"USG";
Usage: $0 [-v] index.sense glosstag_dir
\t-d\t\tDirect links
\t-s name\tSource name
\t-v\t\tBe verbose
example: perl $0 path_to_wordnet/dict/index.sense path_to_glosstag/merged > glossrels.txt

USG
	$msg .= join("\n", @_) if @_;
	chomp($msg);
	die $msg."\n";
}

&usage("Too few parameters") unless @ARGV == 2;

&usage("$ARGV[0] does not exist (or is not a file)") unless -f $ARGV[0];
&usage("$ARGV[1] is not a directory") unless -d $ARGV[1];

my $parser = XML::LibXML->new();
$parser->keep_blanks(0);

my $gdir = $ARGV[1];
my @F = &read_glossfiles($gdir);
&usage($gdir." has no xml files.") unless @F;
my $idx_sense_file = $ARGV[0];
my $D = &read_indexsense($idx_sense_file);

foreach my $f (@F) { &proc_file($f, $D); }

sub proc_file {
	my ($fname, $D) = @_;
	my $tree = $parser->parse_file($fname);
	foreach my $synset ($tree->getDocumentElement->getElementsByTagName('synset')) {
		my $src_syn = $synset->getAttribute("ofs")."-".$synset->getAttribute("pos");
		foreach my $gloss ($synset->getElementsByTagName('gloss')) {
			my @gloss_tgt;
			next unless $gloss->getAttribute('desc') eq 'wsd';
			foreach my $tgt ($gloss->getElementsByTagName('id')) {
				my $tgt_sk = $tgt->getAttribute('sk');
				my $tgt_syn = $D->{$tgt_sk};
				next unless defined $tgt_syn;
				next if $tgt_syn eq $src_syn;
				push @gloss_tgt, $tgt_syn;
			}
			foreach my $tgt_syn (@gloss_tgt) {
				print "u:$src_syn v:$tgt_syn t:gloss".$suff_str."\n";
			}
		}
	}
}

sub read_indexsense {
	my ($idx_sense_file) = @_;
	open(my $fh, $idx_sense_file) || die "Can't open $idx_sense_file:$!\n";
	my $H;
	while (<$fh>) {
		my ($kl, $offset, $sense_n, $f) = split(/\s+/, $_);
		die "Bad index in line $.\n" unless $sense_n;
		$sense_n--;	 # start at 0
		my ($word, $pos) = &split_kl($kl);
		my $cid_str = "$offset-$pos";
		$H->{$kl} = "$offset-$pos";
	}
	return $H;
}

sub read_glossfiles {
	my $gdir = shift;
	my @F;
	opendir my $dh, $gdir or die "Can't open $gdir:$!\n";
	foreach my $fname (grep { /\.xml$/ } readdir $dh) {
		push @F, "$gdir/$fname";
	}
	return @F;
}

sub split_kl {
	my $kl=shift;

	my @Aklpos = ('n', 'v', 'a', 'r', 'a');

	my ($w, $aux) = split(/\%/, $kl);
	my $pos = (split(/\:/, $aux))[0];
	return ($w, $Aklpos[$pos-1]);
}
