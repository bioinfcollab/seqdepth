#!/usr/bin/perl
use 5.10.0;
use strict;
use warnings;
use open qw( :encoding(UTF-8) :std );
use Parallel::ForkManager;
use File::Copy qw(copy move);
use Data::Dumper;
use List::Compare;
use List::Util qw(max);
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempdir);
use Statistics::Basic qw(:all);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
#use Statistics::KernelEstimation qw(pdf extended_tange);
use samples qw(readfiles run_pairmerge nproc $compression $compressor $decompressor);
#use Statistics::Basic::Vector;
#use Statistics::Basic::Median;

$Data::Dumper::Indent   = 1;
$Data::Dumper::Sortkeys = 1;

our $mode = 'table';

sub help_and_exit() {
    print <<"HELP";
    $0 - calculate average reads over multiple bed files
requires --srcdir and --dstdir options

Flags           Description
---------------------
 -h, --help             Display this help
 -s, --srcdir           Directory with bed files
 -d, --dstdir           Result directory, by default used for temporary files
 -t, --tmpdir           Temp directory is used for temporary files if specified
 -n, --nproc            The number of parallel processes
 -v, --verbose          Verbose mode;
HELP
    exit;
}

my $start = time;
my @oldARGS = @ARGV;
my $result_file = 'result_table.txt';
$result_file = 'result_table.txt.gz' if $compression && $compression eq 'gzip';
our ($fnumber, $verbose, $nproc);
my ($sdir,$destdir,$tmpdir,$help);
sub read_command_line_args() {
    my $help;
    Getopt::Long::Configure("bundling");
    GetOptions
    'srcdir|s=s'    => \$sdir,
    'destdir|d=s'   => \$destdir,
    'tmpdir|t=s'    => \$tmpdir,
    'nproc|n=i'     => \$nproc,
    'verbose|v'     => \$verbose,
    'help|h'        => \$help;
    help_and_exit if ! $sdir || ! $destdir;
    help_and_exit if $help;
    return;
}
read_command_line_args;

die "Source directory (-s, --srcdir) $sdir does not exist\n"  if ! -d $sdir;
die "Destination directory (-d, --destdir) $destdir does not exist\n"  if ! -d $destdir;
die "Temporary directory (-t, --tmpdir) $tmpdir does not exist\n"  if $tmpdir && ! -d $tmpdir;
die "-n, --nproc $nproc must be integer\n" if $nproc && ($nproc < '0' || $nproc !~ /\d+/);
#die "The destination dir $destdir is not emtpy\n" if readfiles($destdir);

my $maxproc = nproc($nproc);
our $c = "$compressor -p $maxproc";
our $dc = "$decompressor -fdc";

# parralel core, fork for every pair of files
my $pm = Parallel::ForkManager->new($maxproc);
sub prepare {
    my ($bedpairs,$destdir,$iter) = @_;
    print "Iteration: $iter, combined files:".scalar (keys %$bedpairs).", ";
    my $start = time;
    # data structure retrieval and handling
    PAIRS:
    for my $id (keys %$bedpairs) {
        my $pid = $pm->start and next PAIRS;
        my @files = @{$bedpairs->{$id}};
        my $error = run_pairmerge($id,\@files,$destdir,$iter);
        $pm->finish(0,\$error);
    }
    $pm->wait_all_children;
    my $duration = time - $start;
    say "time: $duration sec";
}

### main
my $h = readfiles($sdir,'0');
foreach my $key(keys %$h) {
    my $f = $h->{$key};
    $fnumber += scalar(@$f);
}
die "Can't find files to work on, exitting\n" if ! $fnumber || $fnumber == 1;
say "Will fork up to $maxproc children";
say "Got $fnumber files, up to ".int(log($fnumber)/log(2)+1.5).' iterations:';

$tmpdir = $tmpdir ? tempdir ( DIR => $tmpdir, CLEANUP => 1) : tempdir (DIR => $destdir, CLEANUP => 1);
say "Using tempdir $tmpdir" if $verbose;
prepare($h,$tmpdir,"1");

#merge samples until only one is left in the destination dir;
my $iter = 1;
while (my $result = readfiles($tmpdir,$iter)){
    $iter ++;
    my @keys = keys %$result;
    # exit if only one "bed" file left in the dstdir directory
    last if scalar @keys == '1' && scalar @{($result->{$keys[0]})} == '1';
    prepare($result,$tmpdir,$iter);
}
my $result = (values %{readfiles($tmpdir,$iter)})[0]->[0];
copy "$result", "${destdir}/${result_file}";
say "Final file: ${destdir}/${result_file}";
