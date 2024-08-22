#!/usr/bin/perl
use Math::GSL::Randist qw(gsl_ran_hypergeometric_pdf);
use List::Util qw(max sum);
use feature qw (say);

#my ($M,$n,$N) = (20,7,12);
#my ($M,$n,$N) = (60e6,30000,30e6);
my ($M,$n,$N) = (60e6,3,30e6);
my ($n1,$n2,$t)  =  ($n,$M-$n,$N);

#aka pmf? https://en.wikipedia.org/wiki/Probability_mass_function
my @p2 = map { $_= gsl_ran_hypergeometric_pdf($_, $n1, $n2, $t)} (0 .. $n);
my @list_prob = map {$_ = sum @p2[$_ .. $#p2]} (0 .. $#p2);
# since the array is sorted, in order to find the min, it is sufficient
# to pickup the last value where condition is met
my $depth;
map { $depth++ if $_ >= '0.8'} @list_prob;
my $prob_expect_depth = sprintf "%.2f", $list_prob[$depth-1];
print $prob_expect_depth;
