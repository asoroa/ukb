# Package to compute similarity of two vectors
#
# Authors: Eneko Agirre and Aitor Soroa
#

package Similarity ;

use Exporter () ;
@ISA = qw(Exporter) ;
@EXPORT = qw(similarity SIM) ;

sub similarity {
    my ($a,$b) = @_ ;
    if ($SIM eq "cos") {
	return  mycos($a,$b) }
    elsif ($SIM eq "dot") {
	return  biderkaduraEsk($a,$b) }
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


