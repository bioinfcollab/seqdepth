package samples;
# apt-get install liblist-compare-perl libparallel-forkmanager-perl libtext-csv-perl libsys-cpu-perl
use strict;
use warnings;
use open qw( :encoding(UTF-8) :std );
use feature qw(say);
use File::Copy qw(copy move);
use List::Compare;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
$Data::Dumper::Indent   = 1;
$Data::Dumper::Sortkeys = 1;

our ( @ISA, @EXPORT, @EXPORT_OK );

BEGIN {
    require Exporter;
    @ISA    = qw(Exporter);
    @EXPORT = qw($compression $compressor $decompressor);
    @EXPORT_OK = qw(readfiles run_pairmerge run_pairmerge_coreutils nproc read_f write_f);
}

#defaults
our $compression  = 'gzip'; # can be gzip or plain
our $compressor   = 'pigz';
our $decompressor = 'gzip';
our $bash  = '/bin/bash -euo pipefail -c';

sub nproc {
    my $argnproc = shift;
    our $hardcoded_MAX_PROCESSES = '10';
    my $rc = eval {
        require Sys::CPU;
        Sys::CPU->import();
        1;
    };
    # try to guess CPU on the host and fallback to hardcoded value, or use provided arguments
    my  $MAX_PROCESSES = $rc ? (Sys::CPU::cpu_count() - 2) : $hardcoded_MAX_PROCESSES;
    $MAX_PROCESSES = $ENV{'NPROC'} if $ENV{'NPROC'};
    $MAX_PROCESSES = $argnproc if $argnproc;
    return $MAX_PROCESSES;
}

sub parsecsv {
    my $rc = eval {
        require Text::CSV;
        Text::CSV->import();
        1;
    };
    die "Can't parse CSV files, Text::CSV is missing" if ! $rc;
    my ($csvfile,$csvopts) =  @_;
    my ($samplename,$csvfilter,$value) = split(/=/,$csvopts);
    my @samples;

    my $csv = Text::CSV->new ({ binary => 1, auto_diag => 1, sep => "\t" });
    open my $fh, "<:encoding(utf8)", "$csvfile" or die "Can't open file: $csvfile, $!";
    $csv->header ($fh);
    while (my $row = $csv->getline_hr ($fh)) {
        die "Can't find in header $samplename" if ! $row->{$samplename};
        push @samples, $row->{$samplename} if $row->{$csvfilter} eq $value;
    }
    close $fh;
    return @samples if @samples;
    return;
}

#read directory with samples
sub readfiles {
    my ($sdir,$iter,$csvfile, $csvopts) = @_;
    opendir my $dh, $sdir || die "Can't open dir $sdir, $!";
    my @bedfiles = grep {/.*\.bed(.gz)?$/} readdir($dh);
    closedir ($dh);

    @bedfiles = sort @bedfiles;
    my @selected_samples;
    ## quick and dirty workaround
    ## read csvfile and select only samples which exist in the source dir and in the list
    ## run only for the iteration 0, when reading the source files
    if ($csvfile && $csvopts && $iter == '0') {
        my @preselected = parsecsv($csvfile,$csvopts);
        for my $sample (@preselected) {
            for my $bed (@bedfiles) {
                push @selected_samples, $bed if $bed =~ m/\A($sample)\.bed\Z/;
            }
        }
        die "Couldn't find any sample from CSV $csvfile in the source directory\n" if ! @selected_samples;
        if (scalar @preselected != scalar @selected_samples ) {
            my $lc = List::Compare->new(\@preselected, \@selected_samples);
            my @complement = $lc->get_complement;
            warn "Can't find the following samples in the dst directory:\n";
            print "$_ " for @complement;
            say "\n\t";
        }
        @bedfiles = @selected_samples;
    }

    my $h;
    while (@bedfiles) {
        my @pair;
        push @pair, $sdir .'/'. shift @bedfiles;
        push @pair, $sdir .'/'. shift @bedfiles if @bedfiles;
        my $id = join('_', getsamples(@pair));
        $h->{$id} = \@pair;
    }
    return $h if $h;
    return;
}

# initial and secondary processing of samples
# more universal sub, but single files can be copied when calculating averages
# instead of stripping the position

sub run_pairmerge_coreutils{
    my ($id,$files,$destdir,$iter) = @_;
    my $flist = join(" ",@$files);
    say "PAIRS: $flist" if $main::verbose;
    my @processed;
    # process every pair and treat separately single ones
    if (scalar @$files == '2') {
        pairmerge_coreutils ($files,$destdir,$id,$iter);
        push @processed, @$files;
    }
    if (scalar (@$files) == '1') {
        my ($sample) = getsamples(@$files);
        say "SINGLE:\t$flist > ${destdir}/${sample}_iter_${iter}.bed" if $main::verbose;
        my $dstfile = "${destdir}/${sample}_iter_$iter.bed";
        local $ENV{F1} = $iter == '1' ? "$main::dc $flist | cut -f2- | sort -k 1b,1 | $main::c > ${dstfile}.gz " : "$main::dc $flist | $main::c > ${dstfile}.gz";
        say "Single cmd: $ENV{F1}" if $main::verbose;
        system($bash .' "$F1"');
        push @processed,$flist;
    }
    # don't delete the source files!
    if ($iter != '1') {
        say "REMOVE: @processed" if @processed && $main::verbose;
        unlink @processed if @processed;
    }
}

sub run_pairmerge{
    my ($id,$files,$destdir,$iter) = @_;
    my $flist = join(" ",@$files);
    say "PAIRS: $flist" if $main::verbose;
    my @processed;
    # process every pair and treat separately single ones
    if (scalar @$files == '2') {
        pairmerge ($files,$destdir,$id,$iter);
        push @processed, @$files;
    }
    if (scalar (@$files) == '1') {
        my ($sample) = getsamples(@$files);
        say "SINGLE:\t$flist > ${destdir}/${sample}_iter_${iter}.bed" if $main::verbose;
        my $dstfile = "${destdir}/${sample}_iter_$iter.bed";
        my $dstfd = write_f($dstfile,$main::compression);
        #open my $dstfd, ">", $dstfile  or die "Can't open file $dstfile, $!\n";
        my $f1 = read_f($flist);
        my $line = read_file_line($f1,$iter);
        while ($line) {
            my ($pos,$value) = @$line;
            say $dstfd "$pos\t".join("\t",@$value);
            $line = read_file_line($f1,$iter);
        }
        close $dstfd;
        push @processed,$flist;
    }
    # don't delete the source files!
    if ($iter != '1') {
        say "REMOVE: @processed" if @processed && $main::verbose;
        unlink @processed if @processed;
    }
}

# initial and secondary processing of samples
sub run_pairmerge_averages{
    my ($id,$files,$destdir,$iter) = @_;
    my $flist = join(" ",@$files);
    say "PAIRS: $flist" if $main::verbose;
    my @processed;
    # process every pair and treat separately single ones
    if (scalar @$files == '2') {
        pairmerge ($files,$destdir,$id,$iter);
        push @processed, @$files;
    }
    if (scalar (@$files) == '1') {
        my ($sample) = getsamples(@$files);
        my $dstfile = "${destdir}/${sample}_iter_$iter.bed";
        say "SINGLE: $dstfile" if $main::verbose;
        copy $flist, $dstfile;
        push @processed,$flist;
    }

    # don't delete the source files!
    if ($iter != '1') {
        say "REMOVE: @processed" if @processed && $main::verbose;
        unlink @processed if @processed;
    }
}

# convert filenames to samples by stripping path and extension
sub getsamples {
    my (@files) = @_;
    my @samples;
    for my $file (@files) {
        $file = (split(/\//,$file))[-1];
        $file =~ s/(.*?)(?:_iter_(?:\d+))?\.bed(\.gz)?$/$1/g;
        push @samples,$file;
    }
    return @samples;
    return;
}

# join 2 bed files by using coreutils join
# massive performace gain
sub pairmerge_coreutils {
    my ($files,$destdir,$id,$iter) = @_;
    $id = md5_hex($id);
    my $dstfile = "${destdir}/${id}_${iter}.bed";
    # some coreutils magick here:
    # since data is huge and perl code appears to be significantly slower it got replaced
    # by hackish coreutils
    # from man 1 join:
    # Important:  FILE1  and  FILE2 must be sorted on the join fields.
    # E.g., use "sort -k 1b,1" if 'join' has no options
    # Note, comparisons honor the rules specified by 'LC_COLLATE'.
    # so we read gzipped files directly in bash process substitution <()
    # cutting the first column - chromosome name on the fly with cut -f2-
    # and sorting by LC_COLLATE as described in man 1 join
    # ENV{} is a cool way to avoid shell escaping hell
    # see here: https://github.com/tpf/perldotcom/issues/202#issue-481945998
    local $ENV{DSTF} = " | $main::c > ${dstfile}.gz";
    local $ENV{JOIN_DEL} = "-t\$'\t'";
    local $ENV{SUBP1} = $iter == '1' ? "<( $main::dc $files->[0] | cut -d\$'\t' -f2- |sort -k 1b,1 )" : "<( $main::dc $files->[0])";
    local $ENV{SUBP2} = $iter == '1' ? "<( $main::dc $files->[1] | cut -d\$'\t' -f2- |sort -k 1b,1 )" : "<( $main::dc $files->[1])";
    say "Pairmerge join cmd: $bash \"join -j 1 $ENV{JOIN_DEL} -a 1 -a 2 $ENV{SUBP1} $ENV{SUBP2} $ENV{DSTF}" if $main::verbose;
    system($bash.' "join -j 1 $JOIN_DEL -a 1 -a 2 $SUBP1 $SUBP2 $DSTF"');
}

# join 2 bed files
sub pairmerge {
    my ($files,$destdir,$id,$iter) = @_;
    $id = md5_hex($id);
    my $dstfile = "${destdir}/${id}_${iter}.bed";
    my $dstfd = write_f($dstfile,$main::compression);
    #open my $dstfd, ">", $dstfile  or die "Can't open file $dstfile, $!\n";
    my $f1 = read_f($files->[0]);
    my $f2 = read_f($files->[1]);

    my $pair1 = read_file_line($f1,$iter);
    my $pair2 = read_file_line($f2,$iter);

    while ($pair1 or $pair2) {
        if ($pair1 && $pair2) {
            if ($pair1->[0] < $pair2->[0]) {
                compute($dstfd,$pair1->[0],$pair1->[1],0);
                $pair1 = read_file_line($f1,$iter);
            }
            elsif ($pair2->[0] < $pair1->[0]) {
                compute($dstfd,$pair2->[0],0,$pair2->[1]);
                $pair2 = read_file_line($f2,$iter);
            }
            else {
                compute($dstfd,$pair1->[0],$pair1->[1],$pair2->[1]);
                $pair1 = read_file_line($f1,$iter);
                $pair2 = read_file_line($f2,$iter);
            }
        }
        elsif ($pair1 && ! $pair2) {
            compute($dstfd,$pair1->[0],$pair1->[1],0);
            $pair1 = read_file_line($f1,$iter);
        }
        elsif (! $pair1 && $pair2) {
            compute($dstfd,$pair2->[0],0,$pair2->[1]);
            $pair2 = read_file_line($f2,$iter);
        }
    }
    close $f1,$f2,$dstfd;
}

# not needed afaik
sub read_f_iter1 {
    my $file = shift;
    die "Can't open file $file: $!\n" if ! -e $file;
    open my $fh, "$main::dc $file| cut -f2-| sort -k 1b,1|" or die "Can't open file $file: $!\n";
    return $fh if $fh;
}

# open plain or gzipped/bgzipped files
sub read_f {
    my $file = shift;
    die "Can't open file $file: $!\n" if ! -e $file;
    # gzip happily parses plain files as well
    open my $fh, "$main::dc $file|" or die "Can't open file $file: $!\n";
    return $fh if $fh;
}

sub write_f {
    my ($file,$mode) = @_;
    open my $fh, ($mode eq 'gzip' ?  "| $main::c -c > ${file}.gz" : "> $file" ) or die "Error writing gziped file $file: $!\n";
    return $fh;
}

sub read_file_line {
    my ($fh,$iter) = @_;
    if ($fh and my $line = readline $fh) {
        chomp $line;
        my @data = split (/\t/, $line);
        shift @data if $iter == '1';
        my $pos = shift(@data);
        return [ ($pos,\@data) ];
    }
    return;
}

# accepts both 2 and 3 values, and ignores the chromosome in 3 value input
# read bed files and return splitted data
sub read_file_line_averages {
    my ($fh) = @_;
    if ($fh and my $line = readline $fh) {
        chomp $line;
        my @data = split (/\t/, $line);
        shift @data if scalar @data == '3';
        return [ @data ];
    }
    return;
}

sub sum{
    my $total;
    $total += $_ for @_;
    return $total;
}

sub compute {
    my ($dstfd,$pos,$val1,$val2) = @_;
    my @result;
    my $result_str;
    if ( ref $val1 eq 'ARRAY') {
        push @result,@$val1;
    }
    else {
        push @result,$val1;
    }
    if ( ref $val2 eq 'ARRAY') {
        push @result,@$val2;
    }
    else {
        push @result,$val2;
    }
    # table is needed
    if ($main::mode  eq 'table') {
        @result = grep {!/^0$/} @result;
        $result_str = join("\t",@result);
    }
    # average case
    else {
        $result_str = sum(@result);
    }
    say $dstfd "$pos\t$result_str";
}

1;
