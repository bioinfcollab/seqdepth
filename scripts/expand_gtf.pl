#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(max);
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use samples qw(nproc read_f write_f $compression $compressor $decompressor);
use Data::Dumper;
use feature qw(say);
$Data::Dumper::Indent   = 1;
$Data::Dumper::Sortkeys = 1;


my $inputf = $ARGV[0] // die "need a single chromosome from a gtf file";
my $result_file = ${inputf};
$result_file =~ s/.*?chr_(chr)?([\d]+|[XxYy])\.txt$/$2.gz/g;
#my $result_file = "${inputf}_expanded.gz";
unlink $result_file if -e $result_file;

my $maxproc = nproc(1);
say "maxproc = $maxproc";
our $c = "$compressor -p $maxproc";
our $dc = "$decompressor -fdc";
my $fd = read_f($inputf);

while (<$fd>) {
    my @list = split("\t",$_);
    next if (! $list[2] || $list[2] ne 'gene');
    my($chr,$p_start,$p_end,$geneid) = @list[0,3,4,8];
    my ($gene_id) = $geneid =~ m/ID=gene-(.*?);/;
        #open my $afd, ">>", $result_file or die "Error writing gziped file $result_file: $!\n";
        open my $afd, ("| $main::c -c >> $result_file") or die "Error writing gziped file $result_file: $!\n";
    #my $dstfd = write_f($result_file,$main::compression);
    while ( $p_start++ <= $p_end) {
            #my $string = join("\t", ($chr_id,$p_start,$gene_id));
            my $string = join("\t", ($p_start,$gene_id));
            say $afd $string;
    }
}

