#!/usr/bin/perl
use strict;
no strict "refs";
use warnings;
use List::Util qw(max);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use samples qw(nproc read_f write_f $compression $compressor $decompressor);
use Data::Dumper;
use feature qw(say);
$Data::Dumper::Indent   = 1;
$Data::Dumper::Sortkeys = 1;


my $inputf = $ARGV[0] // die "need gtf file";

my $maxproc = nproc(1);
say "maxproc = $maxproc";
our $c = "$compressor -p $maxproc";
our $dc = "$decompressor -fdc";
my $fd = read_f($inputf);

my $chr_regex = '^chr([\d]{1,2}|[xyXY])$';

while (my $line = <$fd>) {
    my @list = split("\t",$line);
    next if (! $list[2] || $list[2] ne 'gene');
    my $chr_id = $list[0];
    next if $chr_id =~ m/_/;
    next if $chr_id !~ m/$chr_regex/;
        my $result_file ='chr_'. ${chr_id} . '.txt';
        open my $afd, ">>", $result_file
            or die "Error writing file $result_file: $!\n";
        print $afd $line;
}

