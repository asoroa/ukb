

=head1 $0

=head1 SYNOPSIS

 similarity.pm

 This module exports four objects:

 - SIM: a global variable to know which is the similarity function to
   be used

 - DICT: a global variable which holds the name of the file with the
   word-concept information (only necessary for 1st ave and 1st max
   similarity functions

 - similarity($v1,$v2): a function that applies the corresponding
   similarity function to vectors $v1 and $v2

 - availablesim($string): a function that returns whether the
   similarity function for $string is available

Examples:

   use lib "." ;
   use Similarity ;

   $Similarity::SIM = 'dot' ;
   die if not availablesim($Similarity::SIM) ;
   return similarity($v1,$v2) ;

   $Similarity::SIM = '1stave' ;
   $Similarity::DICT = <~/ukb-0.1.5r1/lkb_sources/30/wnet30_dict.txt> ;
   return similarity($v1,$v2,'cat','dog') ; # 1stave and 1stmax need to know the two NOUNS
                                            # that correspond to $v1 and $v2
                                            # AND the location of the dictionary file

 Author: Eneko Agirre and Aitor Soroa

=cut

package Similarity ;

use Exporter () ;
@ISA = qw(Exporter) ;
@EXPORT = qw(availablesim similarity SIM DICT) ;

sub availablesim {
    my ($s) = @_ ;
    if ($s =~ /^(cos|dot|l1|l2|proj|kl|klsim|1stave|1stmax)$/) {
	return 1}
    else {
	return 0}
}

sub similarity {
    my ($a,$b,$worda,$wordb) = @_ ;
    $worda = "" if not defined $worda ;
    $wordb = "" if not defined $wordb ;
    if ($SIM eq "cos") {
	return  mycos($a,$b) }
    elsif ($SIM eq "dot") {
	return  biderkaduraEsk($a,$b) }
    elsif ($SIM eq "l1") {
	return  -1 * l1norm($a,$b) }
    elsif ($SIM eq "l2") {
	return  -1 * l2norm($a,$b) }
    elsif ($SIM eq "proj") {
	return  averageprojection($a,$b) }
    elsif ($SIM eq "kl") {
	return  -1 * kl($a,$b) }
    elsif ($SIM eq "klsim") {
	return  -1 * klsim($a,$b) }
    elsif ($SIM eq "1stmax") {
	return  firstmax($a,$b,$worda,$wordb) }
    elsif ($SIM eq "1stave") {
	return  firstave($a,$b,$worda,$wordb) }
    else {
	die "wrong similarity $SIM\n" ;
    }
}

sub mycos {
    my ($v1,$v2) = @_ ;
    my ($biderk) ;

    $biderk = biderkaduraEsk ($v1,$v2);
    if ($biderk == 0){
	return 0;
    }
    else{
	return $biderk/(modulua($v1)*modulua($v2));
    }
}

sub modulua {
  my ($v) = @_;
  my ($weight,$mod,$word) ;

  $mod = 0;

  foreach $word (keys %{ $v } ){
      $weight = $v->{$word};
      $mod += ($weight*$weight) ;
  }
  return sqrt($mod);

}

sub biderkaduraEsk {
  my ($v1,$v2) = @_;
  my ($mul,$word) ;
  $mul = 0 ;
  foreach $word (keys %$v1) {
      $mul += $v1->{$word} * $v2->{$word} if defined($v2->{$word}) ;
  }
  return $mul ;
}

sub l1norm {
  my ($v1,$v2) = @_;
  my ($mul) ;
  $mul = 0 ;
  while (my($k,$w) = each %$v1) {
      next if not defined $v2->{$k};
      $mul += abs($w - $v2->{$k});
  }
  return $mul ;
}

sub l2norm {
  my ($v1,$v2) = @_;
  my ($mul) ;
  $mul = 0 ;
  while (my($k,$w) = each %$v1) {
      next if not defined $v2->{$k};
      $mul += ($w - $v2->{$k})**2;
  }
  return sqrt($mul) ;
}

sub averageprojection {
# scalar projection of A onto B is |A|cos(alpha) = A.B/|B|
  my ($v1,$v2) = @_;
  my ($proj1,$proj2) ;
  my $mod1 = modulua($v1) ;
  my $mod2 = modulua($v2) ;
  my $bidesk = biderkaduraEsk($v1,$v2) ;
  $proj1 = $bidesk/$mod2 if ($mod2 > 0) ;
  $proj2 = $bidesk/$mod1 if ($mod1 > 0) ;
  return ($proj1+$proj2)/2 ;
}


# kullback-leibler
#
# D_{KL}(P\|Q) =  -\sum_x p(x) \log q(x) +  \sum_x p(x) \log p(x) \\
#                       =  H(P,Q) - H(P)


# MATLAB IMPLEMENTATION
#
# function KL = kldiv(varValue,pVect1,pVect2,varargin)
# %KLDIV Kullback-Leibler or Jensen-Shannon divergence between two distributions.
# %   KLDIV(X,P1,P2) returns the Kullback-Leibler divergence between two
# %   distributions specified over the M variable values in vector X.  P1 is a
# %   length-M vector of probabilities representing distribution 1, and P2 is a
# %   length-M vector of probabilities representing distribution 2.  Thus, the
# %   probability of value X(i) is P1(i) for distribution 1 and P2(i) for
# %   distribution 2.  The Kullback-Leibler divergence is given by:
# %
# %       KL(P1(x),P2(x)) = sum[P1(x).log(P1(x)/P2(x))]
# %
# %   If X contains duplicate values, there will be an warning message, and these
# %   values will be treated as distinct values.  (I.e., the actual values do
# %   not enter into the computation, but the probabilities for the two
# %   duplicate values will be considered as probabilities corresponding to
# %   two unique values.)  The elements of probability vectors P1 and P2 must
# %   each sum to 1 +/- .00001.
# %
# %   A "log of zero" warning will be thrown for zero-valued probabilities.
# %   Handle this however you wish.  Adding 'eps' or some other small value
# %   to all probabilities seems reasonable.  (Renormalize if necessary.)
# %
# %   KLDIV(X,P1,P2,'sym') returns a symmetric variant of the Kullback-Leibler
# %   divergence, given by [KL(P1,P2)+KL(P2,P1)]/2.  See Johnson and Sinanovic
# %   (2001).
# %
# %   KLDIV(X,P1,P2,'js') returns the Jensen-Shannon divergence, given by
# %   [KL(P1,Q)+KL(P2,Q)]/2, where Q = (P1+P2)/2.  See the Wikipedia article
# %   for "Kullback

# Function KL = kldiv(varValue,pVect1,pVect2,varargin)
# %KLDIV Kullback-Leibler or Jensen-Shannon divergence between two distributions.
# %   KLDIV(X,P1,P2) returns the Kullback-Leibler divergence%   "Jeffrey divergence."  See Rubner et al. (2000).
# %
# %   EXAMPLE:  Let the event set and probability sets be as follow:
# %                X = [1 2 3 3 4]';
# %                P1 = ones(5,1)/5;
# %                P2 = [0 0 .5 .2 .3]' + eps;
# %
# %             Note that the event set here has duplicate values (two 3's). These
# %             will be treated as DISTINCT events by KLDIV. If you want these to
# %             be treated as the SAME event, you will need to collapse their
# %             probabilities together before running KLDIV. One way to do this
# %             is to use UNIQUE to find the set of unique events, and then
# %             iterate over that set, summing probabilities for each instance of
# %             each unique event.  Here, we just leave the duplicate values to be
# %             treated independently (the default):
# %                 KL = kldiv(X,P1,P2);
# %                 KL =
# %                      19.4899
# %
# %             Note also that we avoided the log-of-zero warning by adding 'eps'
# %             to all probability values in P2.  We didn't need to renormalize
# %             because we're still within the sum-to-one tolerance.
# %
# %   REFERENCES:
# %   1) Cover, T.M. and J.A. Thomas. "Elements of Information Theory," Wiley,
# %      1991.
# %   2) Johnson, D.H. and S. Sinanovic. "Symmetrizing the Kullback-Leibler
# %      distance." IEEE Transactions on Information Theory (Submitted).
# %   3) Rubner, Y., Tomasi, C., and Guibas, L. J., 2000. "The Earth Mover's
# %      distance as a metric for image retrieval." International Journal of
# %      Computer Vision, 40(2): 99-121.
# %   4) <a href="matlab:web('http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence','-browser')">Kullback Between two
# %   distributions specified over the M variable values in vector X.  P1 is a
# %   length-M vector of probabilities representing distribution 1, and P2 is a
# %   length-M vector of probabilities representing distribution 2.  Thus, the
# %   probability of value X(i) is P1(i) for distribution 1 and P2(i) for
# %   distribution 2.  The Kullback-Leibler divergence%
# %   See also: MUTUALINFO, ENTROPY

# if ~isequal(unique(varValue),sort(varValue)),
#     warning('KLDIV:duplicates','X contains duplicate values. Treated as distinct values.')
# end
# if ~isequal(size(varValue),size(pVect1)) || ~isequal(size(varValue),size(pVect2)),
#     error('All inputs must have same dimension.')
# end
# % Check probabilities sum to 1:
# if (abs(sum(pVect1) - 1) > .00001) || (abs(sum(pVect2) - 1) > .00001),
#     error('Probablities don''t sum to 1.')
# end

# if ~isempty(varargin),
#     switch varargin{1},
#         case 'js',
#             logQvect = log2((pVect2+pVect1)/2);
#             KL = .5 * (sum(pVect1.*(log2(pVect1)-logQvect)) + ...
#                 sum(pVect2.*(log2(pVect2)-logQvect)));

#         case 'sym',
#             KL1 = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
#             KL2 = sum(pVect2 .* (log2(pVect2)-log2(pVect1)));
#             KL = (KL1+KL2)/2;

#         otherwise
#             error(['Last argument' ' "' varargin{1} '" ' 'not recognized.'])
#     end
# else
#     KL = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
# end


sub kl {
  my ($v1,$v2) = @_;
  my ($mul) ;
  $mul = 0 ;
  # we should go through all keys
  # but:
  # a) if key is in     v1 and in     v2, no problem
  # b)  if key is in     v1 and not in v2, assign 0.00000001
  # c)  if key is not in v1 and in     v2, assign 0.00000001
  # d)  if key is not in v1 and not in v2, mul += 0, so don't do anything
  #
  # for fast computation, we go through keys in v1 (a,b)
  # then go through keys in v1 (if not already in v1) (c)
  my $min = 0.0000000001 ;
  while (my($k,$w) = each %$v1) {
      my $value = $min ;
      if (defined($v2->{$k})) {$value = $v2->{$k} } ;
      $mul += $w * ( log($w)/log(2) - log($value)/log(2)) ;
  }
  while (my($k,$w) = each %$v2) {
      next if $v1->{$k} ;
      $mul += $w * ( log($min)/log(2) - log($w)/log(2)) ;
  }
  return $mul ;
#   while (my($k,$w) = each %$v1) {
#       next if not (defined $v2->{$k}) or ($v2->{$k} == 0);
#       $mul += $w * ( log($w)/log(2) - log($v2->{$k})/log(2)) ;
#   }
}

sub klsim {
  my ($v1,$v2) = @_;
  return (kl($v1,$v2)+kl($v2,$v1))/2;
}


sub kljs { #bukatu gabe
  my ($v1,$v2) = @_;
  my ($mul) ;
  $mul = 0 ;
  while (my($k,$w) = each %$v1) {
      next if not (defined $v2->{$k}) or ($v2->{$k} == 0);
      $mul += $w * ( log($w)/log(2) - log($v2->{$k})/log(2)) ;
  }
  return $mul ;
}



#
# First order formula: look directly nouns in the vectors
#

my %dict ;

sub firstmax {
    my ($v1,$v2,$n1,$n2) = @_ ;
    my ($tot1,$tot2) ;

    die "Similarity.pm::firstave n1 not defined\n" if not $n1;
    die "Similarity.pm::firstave n2 not defined\n" if not $n2;

    loaddict() if not keys %dict;

    $tot1 = $tot2 = 0;
    # look for $n2 in $v1
    foreach my $syn (split(/\s/,$dict{$n2})) {
	next if $syn !~ /-n:/ ; # only nouns!
	$syn =~ s/:\d+$// ;
	$tot1 = $v1->{$syn} if $v1->{$syn} > $tot1 ;
    }
    # look for $n1 in $v2
    foreach my $syn (split(/\s/,$dict{$n1})) {
	next if $syn !~ /-n/ ;
	$syn =~ s/:\d+$// ;
	$tot2 = $v2->{$syn} if $v2->{$syn} > $tot2 ;
    }
    return ($tot1+$tot2)/2 ;
}

sub firstave {
    my ($v1,$v2,$n1,$n2) = @_ ;
    my ($tot,$N) ;

    die "Similarity.pm::firstave n1 not defined\n" if not $n1;
    die "Similarity.pm::firstave n2 not defined\n" if not $n2;

    loaddict() if not keys %dict;

    $tot = 0;
    $N = 0;
    # look for $n2 in $v1
    foreach my $syn (split(/\s/,$dict{$n2})) {
	next if $syn !~ /-n/ ;
	$syn =~ s/:\d+$// ;
	$tot += $v1->{$syn} ;
	$N++ ;
    }
    # look for $n1 in $v2
    foreach my $syn (split(/\s/,$dict{$n1})) {
	next if $syn !~ /-n/ ;
	$syn =~ s/:\d+$// ;
	$tot += $v2->{$syn} ;
	$N++ ;
    }
    if ($N>0) {
	return $tot/$N ;
    } else {
	return 0 ;
    }
}


sub loaddict {

    die "Similarity::DICT not defined\n" if (not defined $DICT) or (not $DICT);
    open(I,$DICT) or die "Similarity.pm::loaddict can't open file $DICT\n"  ;
    while (<I>) {
	chomp ;
	my($lemma,@synsets) = split ;
	$dict{$lemma} = join(" ",@synsets) ;
    }
    close(I) ;
}


(1) ;
