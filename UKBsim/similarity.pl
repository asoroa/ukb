#!/usr/bin/perl -w
#


=head1 $0

=head1 SYNOPSIS

 similarity.pl [options]

 Options:
  --help
  --sim=cos|dot    similarity used: dot product (default), cosine.
  --x=path         path to executable ukb_ppv (default ../bin/)
  --ukbargs="..."  string with the arguments to be passed verbatim to ukb_ppv

  Check ukb_ppv for details on arguments.

 This script takes two nouns from the standard input and outputs their
 similarity as follows:

 1) for word w1, compute the probability distribution of a random walk
    over the given graph where the walk is initialized with w1.
 2) for word w2, compute the probability distribution of a random walk
    over the given graph where the walk is initialized with w2.
 3) return the similarity of the two probability distributions.

 As the initialization of the graph is costly, we first create a file
 with all words, then compute walks for all words, and finally compute
 the similarity.

 Examples:

   Given wnet30_dict.txt and wnet30_rels.txt (as provided and built
   using the ukb release):

   $ cd bin
   $ ./compile_kb -o wnet30_rels.bin ../lkb_sources/30/wnet30_rels.txt
   $ ./similarity.pl --sim dot \\
        --ukbargs "--only_ctx_words --only_synsets --dict_file ../lkb_sources/30/wnet30_dict.txt --kb_binfile wnet30_rels.bin"

 Authors: Eneko Agirre and Aitor Soroa

=cut

use Getopt::Long qw(:config auto_help);
use Pod::Usage;
use File::Temp qw(tempdir) ;
use Similarity ;
use strict ;

my $SIM = 'dot';
my $RUNDIR = '../bin/' ;
my $UKBARGS = "" ;
GetOptions("x=s" => \$RUNDIR,
	   "sim=s" => \$SIM,
	   "ukbargs=s" => \$UKBARGS)
    or
    pod2usage() ;

if ($SIM !~ /^(cos|dot)$/) { pod2usage() } ;
$Similarity::SIM = $SIM ; # Similarity package parameter


my $tmpdir =  tempdir(CLEANUP => 1) ;

# read input pairs
# convert to lowercase
my @pairs ;
while (<STDIN>) {
    chop ;
    next if /^#/ ;
    my($n1,$n2) = split ;
    die "I need word pairs in input\n" if ! $n1 or ! $n2 ;
    push @pairs, [lc $n1, lc $n2] ;
}


# write temporary file in input format of ukb program
my ($n) ;
open(O,">$tmpdir/sim-$$.txt") or die $! ;
foreach my $pair (@pairs)  {
    $n++ ;
    foreach my $x ("0","1",) {
	print O "sim-$$-$n-$pair->[$x]\n" ;
	print O $pair->[$x] . "#n#1#1" ;
	print O "\n" ;
    }
}
close O ;

# call to ukb
system("$RUNDIR/ukb_ppv  $UKBARGS -O  $tmpdir $tmpdir/sim-$$.txt") ;


# collect ppv vectors for each word
my %vector ;
$n=0;
foreach my $pair (@pairs)  {
    $n++ ;
    if ((-e "$tmpdir/sim-$$-$n-$pair->[0].ppv" or -e "$tmpdir/sim-$$-$n-$pair->[0].ppv.bz2")
	and (-e "$tmpdir/sim-$$-$n-$pair->[1].ppv" or -e "$tmpdir/sim-$$-$n-$pair->[1].ppv.bz2")){

	$vector{"0"} = makevector("$tmpdir/sim-$$-$n-$pair->[0].ppv") ;
	$vector{"1"} = makevector("$tmpdir/sim-$$-$n-$pair->[1].ppv") ;

        printf "$pair->[0] $pair->[1] %.10f\n", similarity($vector{"0"},$vector{"1"}) ;
    }
    else {
	print "$pair->[0] $pair->[1] -0\n" ;
    }
}

system "rm $tmpdir/*-$$*" ;

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
