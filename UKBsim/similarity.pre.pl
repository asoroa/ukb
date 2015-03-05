#!/usr/bin/perl -w
#

=head1 $0

=head1 SYNOPSIS

 similarity.pre.pl [options]

 Options:
  --help
  --sim=cos|dot   similarity used: dot product (DEFAULT), cosine.
  --trunc=N       number of dimensions to use (default is to use all)
  --pos=[nvar]+   which pos to check for each lemma (default n)
  --ppvdir=dir    directory where ppv files live

 This script takes two nouns from the standard input and outputs their
 similarity according to the precomputed random walks stored in "dir".

 Examples:

   Given precomputed ppv files (as downloaded from UKB website)

   $ cat >input
   house farm
   ^D
   $ ./similarity.pre.pl --sim dot --ppvdir wn30g.trunc1000/en \\
     <input >output

 Authors: Eneko Agirre and Aitor Soroa

=head1 CHANGES

 20/10/2010

   --trunc option added
   --pos   option added

=cut

use Getopt::Long qw(:config auto_help);
use Pod::Usage ;
use FindBin qw($Bin);
use lib $Bin;
use Similarity ;
use strict ;

my $TRUNC ;
my $SIM = 'dot';
my $PPVDIR = "./" ;
my $POS = "n" ;
GetOptions(
	   "sim=s" => \$SIM,
	   "pos=s" => \$POS,
	   "trunc=s" => \$TRUNC,
	   "ppvdir=s" => \$PPVDIR)
    or
    pod2usage() ;

if ($SIM !~ /^(cos|dot)$/) { pod2usage() } ;
if ($POS !~ /^[nvra]+$/) { pod2usage() } ;
$Similarity::SIM = $SIM ; # Similarity package parameter

if (! $PPVDIR)    { pod2usage() } ;
if (! -d $PPVDIR) { warn "Non existing directory: $PPVDIR\n\n" ; pod2usage(1) } ;


# read input pairs
# convert to lowercase
my %vector ;
while (<STDIN>) {
    chop ;
    next if /^#/ ;
    my($n1,$n2) = split ;
    die "I need word pairs in input\n" if ! $n1 or ! $n2 ;

    $vector{"0"} = makevectorlemma(lc $n1) ;
    $vector{"1"} = makevectorlemma(lc $n2) ;

    if (%{$vector{"0"}} and %{$vector{"1"}}) {
	printf "$n1 $n2 %.10f\n", similarity($vector{"0"},$vector{"1"}) ;
    }
    else {
	print "$n1 $n2 -0\n" ;
    }
}


sub makevectorlemma {
    my ($lemma) = @_ ;
    my $vector = {} ;
    my $lemma0 = substr($lemma,0,1);
    my $lemma1 = substr($lemma,1,1);
    my @files =  <$PPVDIR/$lemma0/$lemma0$lemma1/$lemma.[$POS]*ppv*> ;

    # If $POS is n then only nominal lemmas are checked, otherwise
    # all matching PoS are used, and the vector is a combination

    if (@files ) {
	# Initialize $vector with first file
	$vector = makevectorlemmapos(shift(@files)) ;
	foreach my $f (@files) {
	    # Add vectors of rest of files
	    addvector($vector,
		      makevectorlemmapos($f));
	}
    }
    else {
	warn "No ppv file for $lemma\n" ;
    }
    return $vector ;
}

sub makevectorlemmapos {
    my ($f) = @_ ;
    my $vector = {} ;

    if (! -e $f) {
	warn "$f does not exist\n";
    } elsif ($f =~ /.bz2$/) {
	open(I,"bzcat $f |") or die $!;
    } else {
	open(I,"$f") or die $!;
    }
    my ($i) ;
    while (<I>) {
	chomp ;
	my @F = split(/\t/,$_) ;
	$vector->{$F[0]} = $F[1] ;
	$i++ ;
	last if $TRUNC and $TRUNC < $i ;
    }
    return $vector ;
}


sub addvector {
    my ($v1,$v2) = @_ ;
    while (my($k,$w) = each %$v2) {
	$v1->{$k} = $w + ($v1->{$k} or 0);
    }
}
