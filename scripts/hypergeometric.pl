#!/usr/bin/perl
use strict; use warnings; use feature qw (say);
use Math::GSL::Randist qw(gsl_ran_hypergeometric_pdf);
use List::Util qw(max sum first);
use MCE::Loop;
use Data::Dumper;
use Getopt::Long qw(GetOptions);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use samples qw (readfiles nproc write_f read_f $compression $compressor $decompressor);

sub help_and_exit() {
    print <<"HELP";
    $0 - calculate probabilities and generate result file
    requires: --srcfile, --rfile, --mindepth, --tocheck, --chr, --ctd
    header: chrID position median_value expected_depth probability prob_number_tock

Flags           Description
---------------------
 -h, --help             Display this help
 -s, --srcfile          Source file with i n format "pos\taverage"
 -r, --resultfile       Result file
 -n, --nproc            The number of parallel processes
 -C, --chr              Chromosome
 -m, --mindepth         minDepth
 -c, --tocheck          Depth to check
 -T, --ctd              currentTotalDepth
 -M, --minap            minAcceptedProb
 -v, --verbose          Verbose mode
HELP
    exit;
}

my $start = time;
our ($nproc,$verbose);
my ($sfile,$rfile,$help,$chr,$currentTotalDepth,$minDepth,$numberReadtoCk,$steps,$minAcceptedProb);
$steps = int(10e6);

sub read_command_line_args() {
    my $help;
    Getopt::Long::Configure("bundling");
    GetOptions
    'srcfile|s=s'    => \$sfile,
    'rfile|r=s'      => \$rfile,
    'mindepth|m=s'   => \$minDepth,
    'chr|C=s'        => \$chr,
    'tocheck|c=i'    => \$numberReadtoCk,
    'nproc|n=i'      => \$nproc,
    'verbose|v'      => \$verbose,
    'ctd|T=s'        => \$currentTotalDepth,
    'minap|M=s'      => \$minAcceptedProb,
    'help|h'         => \$help;
    help_and_exit if ! $sfile|| ! $rfile || ! $chr || !$minDepth || !$minAcceptedProb || !$currentTotalDepth || !$numberReadtoCk;
    help_and_exit if $help;
    return;
}
read_command_line_args;

die "-n, --nproc $nproc must be integer\n" if $nproc && ($nproc < '0' || $nproc !~ /\d+/);
die "Can't read the source file (-s, --srcfile): $sfile " if ! -f $sfile;

my $maxproc = nproc($nproc);
our $c = "$compressor -p $maxproc";
our $dc = "$decompressor -fdc";


say "Forking up to $maxproc" if $verbose;

sub up_to_min_accepted_prob {
    my ($n) = @_;
    my ($n1,$n2,$t) = ($n,$currentTotalDepth-$n,$minDepth);
    my $cache;
    my $index;
    my $accum;
    for my $x ( reverse 0 .. $n) {
        my $prob  = gsl_ran_hypergeometric_pdf($x,$n1,$n2,$t);
        $cache->{$x} = $prob;
        $accum += $prob;
        $index = $x;
        last if $accum >= $minAcceptedProb;
    }
    return ($index,$accum,$cache);
}

sub prob_for_number_read_toCk {
    my ($n,$s,$cache) = @_;
    my ($n1,$n2,$t) = ($n,$currentTotalDepth-$n,$minDepth);
    my $accum;
    my $prob;
    return 0 if $s > $n;
    for my $x ($s .. $n) {
        $prob = $cache->{$x} ? $cache->{$x} : gsl_ran_hypergeometric_pdf($x,$n1,$n2,$t);
        $accum += $prob;
    }
    return $accum;
}

sub calc {
    my ($n,$s) = @_;
    # ensure $n is int, required for gsl_ran_hypergeomtric function
    # kind of ceiling
    $n = int($n+0.5);
    my ($n1,$n2,$t) = ($n,$currentTotalDepth-$n,$minDepth);
    my $result->{'input'} = $n;
    my ($expect_depth,$prob_expect_depth,$cache) = up_to_min_accepted_prob($n);
    $result->{'expect_depth'} = $expect_depth;
    $result->{'prob_expect_depth'} = sprintf "%.2f", $prob_expect_depth;
    if ($s) {
        my $prob_number_tock = prob_for_number_read_toCk($n,$s,$cache);
        $result->{'prob_number_tock'} = sprintf "%.2f", $prob_number_tock;
    }
    return $result;
}

sub collect_uniq_vals_from_{
    say "Collecting the uniq values from $sfile" if $verbose;
    my $file = shift;
    my %vals;
    my $fd = read_f($file);
    while (<$fd>) {
        my ($pos,$value) = split;
        $vals{$value} = 1;
    }
    my @uniqvals = keys %vals;
    close $fd;
    return @uniqvals if @uniqvals;
    return;
}

my $fd = read_f($sfile);
my @uniqvals = collect_uniq_vals_from_($sfile);
say "Size of \@uniqvals is " .@uniqvals if $verbose;

MCE::Loop->init(
    max_workers => $maxproc, chunk_size => 1,
);
my %h = mce_loop {
    my ($mce, $chunk_ref, $chunk_id) = @_;
    my %ret;
    $ret{$_} = calc($_, $numberReadtoCk);
    MCE->gather(%ret);
} \@uniqvals;


MCE::Loop->finish;

my $rfd = write_f($rfile,$compression);

say "Final step, generating results file: $rfile" if $verbose;
while  (<$fd>) {
    my ($pos,$value) = split;
    my $str = $chr."\t$pos"."\t".$value."\t".$h{$value}->{'expect_depth'}."\t".$h{$value}->{'prob_expect_depth'};
    $str = $numberReadtoCk ? $str . "\t".$h{$value}->{'prob_number_tock'} : $str;
    say $rfd $str;
#    if ($numberReadtoCk) {
#        say  $rfd $chr."\t$pos"."\t".$value."\t".$h{$value}->{'expect_depth'}."\t".$h{$value}->{'prob_expect_depth'}."\t".$h{$value}->{'prob_number_tock'};
#    } else {
#        say  $rfd $chr."\t$pos"."\t".$value."\t".$h{$value}->{'expect_depth'}."\t".$h{$value}->{'prob_expect_depth'};
#    }
}
close $rfd;
