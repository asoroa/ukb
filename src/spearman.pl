#!/usr/bin/perl

use File::Temp;

die "Usage: $0 file1 file2\n" unless @ARGV == 2;

my $h1 = &rfile($ARGV[0]);
my $h2 = &rfile($ARGV[1]);

my $hi = &intersect($h1, $h2); # { key => index }
print STDERR "intersect: ". scalar ( keys %{ $hi } ). " elements\n";

my @A;
my @B;
foreach my $k ( sort { $hi->{$b} <=> $hi->{$a} } keys %{ $hi } ) {
        push @A, $h1->{$k};
        push @B, $h2->{$k};
}

#my $ftmp = File::Temp->new();
#die "Can not create temporal file:$!\n" unless $ftmp;
open(my $ftmp, ">/tmp/R.$$.in");
print $ftmp "a <- c(" . join(",",@A) . ")\n" ;
print $ftmp "b <- c(" . join(",",@B) . ")\n" ;
print $ftmp "cor.test(a, b, method=\"pearson\")\n" ;
print $ftmp "cor.test(a, b, method=\"spearman\")\n" ;
close($ftmp);

system "(cat /tmp/R.$$.in | R --vanilla > /tmp/R.$$.out) > /dev/null 2>&1" ;

my ($pearson, $conf1, $conf2);
my $spearman;
my $done = 0;
open(I,"/tmp/R.$$.out") or die $! ;
while(<I>) {
    if (/Spearman\'s rank correlation/) {
	while (<I>) {
	    next if ! /^\s*(-?[01]\.?\d*)\s*$/ ;
	    $spearman = $1 ;
		$done += 1;
	    last ;
	}
    }
    if (/Pearson\'s product-moment/) {
	while (<I>) {
	    if (/^\s*(0\.\d+)\s+(0\.\d+)\s*$/) {
		$conf1 = $1 ;
		$conf2 = $2 ;
	    }
	    if  (/^\s*(-?[01]\.?\d*)\s*$/ ) {
		$pearson = $1 ;
		$done += 1;
		last ;
	    }
	}
    }
    last if $done >= 2;
}

system "rm /tmp/R.$$.*" ;

#print "Pearson: $pearson ($conf1, $conf2)\n" ;
print "Spearman: $spearman\n" ;

sub intersect {

        my ($ha, $hb) = @_;

        my $hi = {};
        my ($h1, $h2) = ($ha, $hb);
        ($h1, $h2) = ($hb, $ha)        if scalar keys %{ $hb } < scalar keys %{ $ha };
        my $i = 0;
        while( my ($k, $v) = each %{ $h1 } ) {
                next unless $h2->{$k};
                $hi->{$k} = $i;
                $i++;
        }
        return $hi;
}

sub rfile {

        my ($fname) = @_;
        my $h = {};
        open(my $fh, $fname) or die "Can not open $fname:$!\n";
        while(<$fh>) {
                next if /^\s*$/;
                next if /^!!/;
                chomp;
                my ($s, $v) = split(/\s+/, $_);
                $h->{$s} = $v;
        }
        return $h;
}
