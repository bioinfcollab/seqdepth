use strict;
use warnings;
use Test::More tests => 1;

my $result = `perl t/hypergeom.pl`;
is($result,'0.87',"Count hypergeometric probability density function");
