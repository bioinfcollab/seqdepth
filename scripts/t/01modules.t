use feature qw(say);
use strict;
use warnings;
use Test::More tests => 10;

require_ok ('Getopt::Long');
require_ok ('Parallel::ForkManager');
require_ok ('Digest::MD5');
require_ok ('MCE::Loop');
require_ok ('List::Compare');
require_ok ('File::Temp');
require_ok ('Text::CSV');
require_ok ('Sys::CPU');
require_ok ('Statistics::Basic');
require_ok ('Math::GSL::Randist');
