#!/usr/bin/perl

use strict ;
use FindBin qw($Bin);
use lib "$Bin";

use ukbProtocol;

my $sess = new ukbProtocol(6969);

my $str;
my $ok;

open(my $fh, $ARGV[0]) or die;

my $l_n = 0;
while(my $line = <$fh>) {
	$l_n++;
	$sess->send($l_n);
	$str = $sess->receive();
	$sess->send($line);
	if ($line ne "stop") {
		$str = $sess->receive();
		print STDOUT "$str";
	}
}
