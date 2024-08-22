use feature qw(say);
use strict;
use warnings;
use Test::More tests => 2;
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use testing qw(t_readfiles t_openfd);

=pod

=head1 Test Logic

 - generate random number of bed files in a temp directory
 - run fast code building table over the files
 - run median and hypergeometric code
 - fail tests if non zero return values
 - cleanup tmpdir

=cut

my $table_code = './samples_table.pl';
my $median_code= './median.pl';
my $hyper_code= './hypergeometric.pl';
my $gendata = 't/gendata.pl';
my $datadir = 't/data';
my $tmp = tempdir( DIR => $datadir, CLEANUP => 1 );
system "$gendata | gzip -c > ${tmp}/sample_${_}.bed.gz" for (rand (3) .. (rand(3) +10));
my $fnumber = scalar(t_readfiles($tmp));
my $result = `$table_code -s $tmp -d $tmp`;
#say "$median_code -i $tmp/result_table.txt.gz -N $fnumber -n 4 -d $tmp";
my @median_result = `$median_code -i "$tmp/result_table.txt.gz" -N $fnumber -n 4 -d $tmp`;
#say "$hyper_code -s \"$tmp/result_median.txt.gz\" -C 1 -n 4 -c 20 --ctd 30e6 -m 30e6 --mindepth 15e6 -r \"$tmp/hypergeometric.txt\"";
my @hypergeometric_result = `$hyper_code -s "$tmp/result_median.txt.gz" -C 1 -n 4 -c 20 --ctd 30e6 -M 30e6 --mindepth 15e6 -r "$tmp/hypergeometric.txt"`;
is (@median_result, 0, "Median code exited cleanly");
is (@hypergeometric_result, 0, "Hypergeometric code exited cleanly");
