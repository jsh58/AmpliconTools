#!/usr/bin/perl

# John M. Gaspar
# June 2015

# Reformat primers for input to ipcress.

use strict;
use warnings;

sub usage {
  print q(Usage: perl formatPrimers.pl  <infile>  <outfile>
  Required:
    <infile>   File listing primer and target sequences (produced by getPrimers.pl)
    <outfile>  Output file listing primers, reformatted for use by ipcress
  Optional:
    <minLen>   Minimum amplicon length (def. 50)
    <maxLen>   Maximum amplicon length (def. 300)
);
  exit;
}

usage() if (scalar @ARGV < 2 || $ARGV[0] eq "-h");

# open files
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(OUT, ">$ARGV[1]");

# load min/max lengths
my $min = 50;
$min = $ARGV[2] if (scalar @ARGV > 2);
my $max = 300;
$max = $ARGV[3] if (scalar @ARGV > 3);

# reformat lines
while (my $line = <IN>) {
  chomp $line;
  my @spl = split(',', $line);
  next if (scalar @spl < 3);
  my $rev = reverse $spl[2];
  $rev =~ tr/ACGT/TGCA/;
  print OUT "$spl[0]\t$spl[1]\t$rev\t$min\t$max\n";
}
close IN;
close OUT;
