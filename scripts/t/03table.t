use feature qw(say);
use strict;
use warnings;
use Test::More tests => 1;
use File::Temp qw(tempdir);
use experimental 'smartmatch';
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use testing qw(t_readfiles t_openfd);

=pod

=head1 Test Logic

 - generate random number of bed files in a temp directory
 - run fast code building table over the files
 - run slow code collecting values into a hash from the files
 - compare the results, fail test if the results are different
 - cleanup tmpdir

=cut

sub readdata{
    my ($file,$type) = @_;
    my $fd = t_openfd($file);
    my @result;
    while (<$fd>) {
        my @data = split;
        shift @data if $type && $type eq 'bedfiles';
        my $newdata = join("\t",@data);
        push @result,$newdata;
    }
    return @result if @result;
    return;
}

sub table{
    my ($tmp) = @_;
    my %h;
    for my $file (t_readfiles($tmp)) {
        my @data = readdata("${tmp}/${file}",'bedfiles');
        map{ my ($pos,$read) = split; push @{$h{$pos}},$read} @data;
    }
    return \%h;
}

my $tool = './samples_table_pp.pl';
my $gendata = 't/gendata.pl';
my $datadir = 't/data';
my $tmp = tempdir( DIR => $datadir, CLEANUP => 1 );
system "$gendata | gzip -c > ${tmp}/sample_${_}.bed.gz" for (rand (3) .. (rand(3) +10));

my $result = `$tool -s $tmp -d $tmp`;

my @files = t_readfiles($tmp);

my $h = table($tmp);
# compress files by default, hence .gz
my @result = (readdata("${tmp}/result_table.txt.gz"));

my $h2;
map {my ($pos,@vals) = split; push @{$h2->{$pos}},@vals} @result;

my $fail;
for my $key (keys %$h) {
    for my $k2 ( keys %$h2) {
        my @table_data = sort (grep {!/^0$/} @{$h->{$k2}});
        my @result_arr = sort (@{$h2->{$k2}});
        $fail = '1' if ! (@table_data ~~ @result_arr);
    }
}

is ($fail, undef, "Test that table build with different methods have same values across all positions");
