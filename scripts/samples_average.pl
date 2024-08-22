#!/usr/bin/perl
# calculate average for position values without creating intermediate table file
use strict;
use warnings;
use feature qw(say);
use Parallel::ForkManager;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempdir);
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use samples qw(readfiles run_pairmerge nproc read_f write_f $compression $compressor $decompressor);

$Data::Dumper::Indent   = 1;
$Data::Dumper::Sortkeys = 1;
our $mode = 'average';

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
 -f, --filter           Filter values smaller than defined
 -v, --verbose          Verbose mode;
HELP
    exit;
}

my $start = time;
my @oldARGS = @ARGV;
my $result_file = 'result_averages.txt';
# this vars are visible in the lib/samples.pm module
our ($fnumber, $verbose, $nproc);
my ($sdir,$destdir,$tmpdir,$help,$filter);
sub read_command_line_args() {
    my $help;
    Getopt::Long::Configure("bundling");
    GetOptions
    'srcdir|s=s'    => \$sdir,
    'destdir|d=s'   => \$destdir,
    'tmpdir|t=s'    => \$tmpdir,
    'nproc|n=i'     => \$nproc,
    'filter|f=f'    => \$filter,
    'verbose|v'     => \$verbose,
    'help|h'        => \$help;
    help_and_exit if ! $sdir || ! $destdir;
    help_and_exit if $help;
    return;
}
read_command_line_args;

die "Source directory (-s, --srcdir) $sdir does not exist\n"  if ! -e $sdir;
die "Destination directory (-d, --destdir) $destdir does not exist\n"  if ! -e $destdir;
die "Temporary directory (-t, --tmpdir) $tmpdir does not exist\n"  if $tmpdir && ! -e $tmpdir;
die "-n, --nproc $nproc must be integer\n" if $nproc && ($nproc < '0' || $nproc !~ /\d+/);


my $maxproc = nproc($nproc);
our $c = "$compressor -p $maxproc";
our $dc = "$decompressor -fdc";

# parralel core, fork for every pair of files
my $pm = Parallel::ForkManager->new($maxproc);
sub prepare {
    my ($bedpairs,$destdir,$iter) = @_;
    print "Iteration: $iter, combined files:".scalar (keys %$bedpairs).", ";
    my $start = time;

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
say "Will fork up to $maxproc children";
say "Got $fnumber files, ".int(log($fnumber)/log(2)+0.5).' iterations:';

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
my $fd = read_f($result);
say "Got compression = $compression";
my $rd = write_f("${destdir}/${result_file}",$compression);

while (<$fd>) {
    my ($pos,$value) = split();
    my $average = int($value/$fnumber + 0.5);
    next if $filter && $average  <= $filter;
    say $rd "$pos\t$average";
}
unlink $result if -e $result;
close $rd;
say "Result file: ${destdir}/${result_file}";
my $duration = time - $start;
say "run took: ~ $duration sec";
