#!/usr/bin/perl
use strict;
use warnings;
use MCE;
use MCE::Loop;
use List::Util qw(max);
use Getopt::Long qw(GetOptions);
use Statistics::Basic qw(:all);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use samples qw(nproc read_f write_f $compression $compressor $decompressor);
use Data::Dumper;
use feature qw(say);
$Data::Dumper::Indent   = 1;
$Data::Dumper::Sortkeys = 1;

sub help_and_exit() {
    print <<"HELP";
    $0 - calculate median over table of vectors
requires --input, --dstdir and --fnumber options

Flags           Description
---------------------
 -h, --help             Display this help
 -i, --input            Input table file
 -d, --dstdir           Destination directory
 -n, --nproc            The number of parallel processes
 -N, --fnumber          Number of samples
HELP
    exit;
}

my $start = time;
my @oldARGS = @ARGV;
our ($fnumber,$nproc);
my ($inputf,$destdir,$help,$filter);
sub read_command_line_args() {
    my $help;
    Getopt::Long::Configure("bundling");
    GetOptions
    'input|i=s'     => \$inputf,
    'destdir|d=s'   => \$destdir,
    'fnumber|N=i'   => \$fnumber,
    'nproc|n=i'     => \$nproc,
    'help|h'        => \$help;
    help_and_exit if ! $inputf || ! $destdir || ! $fnumber;
    help_and_exit if $help;
    return;
}
read_command_line_args;

die "Input file (-i, --input) $inputf does not exist\n"  if ! -f $inputf;
die "Destination directory (-d, --destdir) $destdir does not exist\n"  if ! -d $destdir;
die "-n, --nproc $nproc must be integer\n" if $nproc && ($nproc < '0' || $nproc !~ /\d+/);

my $maxproc = nproc($nproc);
our $c = "$compressor -p $maxproc";
our $dc = "$decompressor -fdc";

my $result_file = "${destdir}/result_median.txt";
my $rd = write_f($result_file,$compression);

MCE::Loop->init(
    max_workers => nproc(), chunk_size => 1,
    user_output => sub {
        print $rd $_[0];
    },
);

my $fd = read_f($inputf);

mce_loop_f {
    my ($mce, $chunk_ref, $chunk_id) = @_;
    my @values = split();
    my $pos = shift @values;
    return if ! @values;
    my @svalues = sort(@values);
    my ($mypos,$result) = compute(\@svalues,$fnumber,$pos);
    #MCE->say("$mypos\t".sprintf("%.2f",$result)."\t".sprintf("%.2f",$pval));
    MCE->say("$mypos\t".sprintf("%u",$result));
} $fd;

MCE::Loop->finish;

sub round {
    my ($val)  = @_;
    return int ($val + 0.5);
}

# todo, implement cache for shapiro_wilk test?
sub compute{
    my ($svalues,$fnumber,$pos) = @_;
    my $vector = vector(@$svalues);
    $vector->set_size($fnumber);
    # median returns an object, query() method suppose to return the number
    my $result = round(median($vector->{v})->query());
    return ($pos,$result);
    return "grrr" if ! $result;
}
