use feature qw(say);
use strict;
use warnings;
use Test::More tests => 1;
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use testing qw(t_readfiles t_openfd);

=pod

=head1 Test Logic

 - generate random number of bed files in a temp directory
 - run fast code counting average over the files
 - run slow code counting average over the files
 - compare the results, fail test if the results are different
 - cleanup tmpdir

=cut

sub readdata{
    my $file = shift;
    my $fd = t_openfd($file);
    my @result;
    while (<$fd>) {
        my @data = split;
        shift @data if scalar @data == '3';
        push @result,"$data[0]\t$data[1]";
    }
    return @result if @result;
    return;
}

sub average {
    my ($nr,$tmp) = @_;
    my %h;
    for my $file (t_readfiles($tmp)) {
        my @data = readdata("${tmp}/${file}");
        map{ my ($pos,$read) = split; $h{$pos}+= $read } @data;
    }
    map{ $h{$_} = int ($h{$_}/${nr} +0.5)} keys (%h);
    return \%h;
}

my $tool = './samples_average.pl';
my $gendata = 't/gendata.pl';
my $datadir = 't/data';
my $tmp = tempdir( DIR => $datadir, CLEANUP =>1 );
system "$gendata > ${tmp}/sample_${_}.bed" for (rand (3) .. (rand(3) +3));

my $result = `$tool -s $tmp -d $tmp`;

my @files = t_readfiles($tmp);
my $nr = scalar(@files);

my $h = average($nr,$tmp);
my @result = (readdata("${tmp}/result_averages.txt.gz"));

my $h2;
map {my ($pos,$read) = split; $h2->{$pos} = $read} @result;

my $fail;
# compare both hashes
for my $key (keys %$h) {
    for my $k2 ( keys %$h2) {
        $fail = '1' if $h->{$k2} != $h2->{$k2};
    }
}

is ($fail, undef, "Test that average values calculated with different methods are the same across all positions");
