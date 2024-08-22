package testing;
use strict;
use warnings;
use PerlIO::gzip;

our ( @ISA, @EXPORT, @EXPORT_OK );

BEGIN {
    require Exporter;
    @ISA    = qw(Exporter);
    @EXPORT = qw();
    @EXPORT_OK = qw(t_readfiles t_openfd);
}

sub t_openfd {
    my $file = shift;
    my $type = $file =~ m/.*?\.(bed|txt)\.gz/ ? 'gzip' : 'plain';
    my $filehandle;
    if ($type  eq 'gzip') {
        open $filehandle, '-|', '/bin/gunzip', '-c', $file or die "Can't open file $file: $!\n";
    }
    else {
        open $filehandle, "<", $file or die "Openfd, can't open file $file: $!";
    }
    return $filehandle if $filehandle;
    return;
}

sub t_readfiles {
    my $dir = shift;
    opendir my $dh, $dir || die "Can't open dir $dir, $!";
    my @bedfiles = grep {/.*\.bed(.gz)?$/} readdir($dh);
    closedir ($dh);
    @bedfiles = sort @bedfiles;
    return @bedfiles if @bedfiles;
    return;
}

1;
