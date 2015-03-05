#!/usr/bin/perl -w
#


=head1 $0

=head1 SYNOPSIS

 similarity.pl [options]

 Options:
  --help
  --sim=cos|dot    similarity used: dot product (default), cosine.
  --x=program      path and name of executable ukb_ppv (default ../bin/ukb_
  --ukbargs="..."  string with the arguments to be passed verbatim to ukb_ppv
  --tmpdir="..."   string with the path to a temporal directory
  --debug=0        0 or 1

  Check ukb_ppv for details on arguments.

 This script takes two nouns from the standard input and outputs their
 similarity as follows:

 1) for noun w1, compute the probability distribution of a random walk
    over the given graph where the walk is initialized with w1.
 2) for noun w2, compute the probability distribution of a random walk
    over the given graph where the walk is initialized with w2.
 3) return the similarity of the two probability distributions.

 As the initialization of the graph is costly, we first create a file
 with all nouns, then compute walks for all nouns, and finally compute
 the similarity.

 NOTE:
 - if you need to select other parts of speech or specific synsets
   you can use the # notations, e.g. house#v#w1#1 (a verb) or
   08249817-n##c1#2 (concept) or 00980806-v##c2#2#0.6 (concept with
   weight).


 - if you want to compare similarity across all parts of speech, include
   --nopos among ukbargs, e.g.

   --ukbargs "--nopos ..."


 Examples:

   Given wnet30_dict.txt and wnet30_rels.txt (as provided and built
   using the ukb release):

   $ cd bin
   $ ./compile_kb -o wnet30_rels.bin ../lkb_sources/30/wnet30_rels.txt
   $ cat >input
   house farm
   house#n#1#1 house#v#1#1
   ^D
   $ ./similarity.pl --sim cos \\
        --ukbargs "--dict_file ../lkb_sources/30/wnet30_dict.txt --kb_binfile wnet30_rels.bin" \\
   <input >output


 Author: Eneko Agirre and Aitor Soroa


=cut

use Getopt::Long qw(:config auto_help);
use Pod::Usage;
use File::Temp qw(tempdir) ;
use FindBin qw($Bin);
use lib $Bin;
use Similarity ;
use strict ;

# default parameter
my $SIM = 'cos';
my $RUN = '../bin/ukb_ppv' ;
my $UKBARGS = "" ;
my $TMPDIR = "." ;
my $DEBUG = 0 ;

GetOptions("x=s" => \$RUN,
	   "debug=i" => \$DEBUG,
	   "sim=s" => \$SIM,
	   "ukbargs=s" => \$UKBARGS,
           "tmpdir=s" => \$TMPDIR)
    or
    pod2usage() ;

$Similarity::SIM  = $SIM ;               # Similarity package parameter
pod2usage() if not availablesim($SIM) ;  # check if $SIM is available

# Dictionary in UKBARGS needed for some similarity computations
($Similarity::DICT) = ($UKBARGS =~ /--dict_file (\S*)/);
($Similarity::DICT) = ($UKBARGS =~ /-D (\S*)/) if not defined $Similarity::DICT ;


my $tmpdir =  tempdir(DIR => $TMPDIR, CLEANUP => not $DEBUG) ;

# read input pairs
my @pairs ;
while (<STDIN>) {
    chop ;
    next if /^#/ ;
    my($n1,$n2) = split ;
    die "I need word pairs in input\n" if ! $n1 or ! $n2 ;
    push @pairs, [$n1, $n2] ;
}


# write temporary file in input format of ukb program
# convert to lowercase
my ($n) ;
open(O,">$tmpdir/00-sim.txt") or die $! ;
foreach my $pair (@pairs)  {
    $n++ ;
    foreach my $x ("0","1",) {
	print O "$n-$pair->[$x]\n" ;
	print O $pair->[$x] ;
	print O "#n#1#1" unless $pair->[$x] =~ /#/ ;
	print O "\n" ;
    }
}
close O ;

# call to ukb
my $status = system("$RUN $UKBARGS -O  $tmpdir $tmpdir/00-sim.txt") ;
die "Program $RUN $UKBARGS -O  $tmpdir $tmpdir/00-sim.txt exited unexpectedly: $?" unless $status == 0 ;

# collect ppv vectors for each word
my %vector ;
$n=0;
foreach my $pair (@pairs)  {
    $n++ ;
    if ((-e "$tmpdir/$n-$pair->[0].ppv" or -e "$tmpdir/$n-$pair->[0].ppv.bz2")
	and (-e "$tmpdir/$n-$pair->[1].ppv" or -e "$tmpdir/$n-$pair->[1].ppv.bz2")){

	$vector{"0"} = makevector("$tmpdir/$n-$pair->[0].ppv") ;
	$vector{"1"} = makevector("$tmpdir/$n-$pair->[1].ppv") ;

        printf "$pair->[0] $pair->[1] %.10f\n", similarity($vector{"0"},$vector{"1"}) ;
    }
    else {
	print "$pair->[0] $pair->[1] -0\n" ;
    }
}


print STDERR "Please see .ppv files in $tmpdir (and remove afterwards)\n" if $DEBUG ;

sub makevector {
    my ($f) = @_ ;
    my $vector = {} ;
    if (-e $f ) {
	open(I,$f) or die $!;
    } elsif (-e "$f.bz2") {
	open(I,"bzcat $f.bz2 |") or die $!;
    } else {
	die "$f does not exist\n";
    }
    while (<I>) {
	chomp ;
	my @F = split(/\t/,$_) ;
	$vector->{$F[0]} = $F[1] ;
    }
    return $vector ;
}
