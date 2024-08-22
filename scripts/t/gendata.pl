#!/usr/bin/perl
use strict;
use warnings;
use feature qw (say);

my $base = '10';
my %data = map {$_ => int (sin(rand(10))+0.5+int(rand 100))} ((rand($base) + $base) .. (rand($base) + 2*$base));
my %empty = map {$_ => ''} ((rand($base) + $base) .. (rand($base) + 2*$base));
say "chr1\t$_\t$data{$_}" for sort keys %data;
#gen some emtpy values
#say "chr1\t$_\t$empty{$_}" for sort keys %empty;
