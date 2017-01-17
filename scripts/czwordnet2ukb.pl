
use strict;
use XML::LibXML;

die "usage: $0 CzechWordNet1.9.xml\n" unless @ARGV == 1;

my $fname = shift @ARGV;

my $gfname = "graph.txt";
my $dfname = "dict.txt";

open (my $gfh, ">$gfname") or die "Can't create $gfname:$!\n";
open (my $dfh, ">$dfname") or die "Can't create $dfname:$!\n";
binmode $dfh, ":utf8";

my $parser = XML::LibXML->new();
$parser->keep_blanks(0);

my $doc = $parser->parse_file($fname);
my $doc_elem = $doc->getDocumentElement;

my %Dict; # { word => { pos => [sense, sense] } }

foreach my $synset_elem ($doc->getDocumentElement->findnodes('/SYNSETS/SYNSET')) {
    my ($pos) = map { $_->textContent } $synset_elem->getChildrenByTagName('POS');
    next unless defined $pos;
    my ($cid) = map { $_->textContent } $synset_elem->getChildrenByTagName('ID');
    next unless defined $cid;
    # Dict
    my ($syn_elem) =  $synset_elem->getChildrenByTagName('SYNONYM');
    foreach  my $lit_elem ($syn_elem->getChildrenByTagName('LITERAL')) {
        my $word = $lit_elem->firstChild->textContent;
        next unless defined $word;
        my ($sense) = map { $_->textContent } $lit_elem->getChildrenByTagName('SENSE');
        next unless defined $sense;
        $Dict{$word}->{$pos}->[$sense] = $cid;
    }
    # print graph relations
    foreach my $ilr_elem ($synset_elem->findnodes('./ILR')) {
        my $v = $ilr_elem->firstChild->textContent;
        next unless defined $v;
        print $gfh "u:$cid v:$v\n";
    }
}

# prind dictionary

foreach my $hw (sort { $Dict{$a} cmp $Dict{$b} } keys %Dict) {
    my $h = $Dict{$hw};
    my @A;
    foreach my $pos (sort { $a cmp $b } keys %{ $h }) {
        foreach (@{ $h->{$pos} }) {
            push @A, $_ if defined $_;
        }
    }
    next unless @A;
    print $dfh $hw." ".join(" ", @A)."\n";
}
