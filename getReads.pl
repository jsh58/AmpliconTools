#!/usr/bin/perl

# John M. Gaspar
# July 2014

# Retrieve a set of reads from a fastq file.

use strict;
use warnings;

sub usage {
  print q(Usage: perl getReads.pl  <infile1>  <infile2>  <outfile>
  Required:
    <infile1>  FASTQ file listing reads to be retrieved (used only
                 to get read headers)
    <infile2>  FASTQ file listing original reads
    <outfile>  Output file for reads
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(FQ, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]") || die "Cannot open $ARGV[2] for writing\n";

# load read headers from first file
my %un;
my $count = 0;
while (my $line = <IN>) {
  chomp $line;
  next if (substr($line, 0, 1) ne "@");

  my @spl = split(" ", $line);
  $un{$spl[0]} = 1;

  for (my $x = 0; $x < 3; $x++) {
    my $waste = <IN>;
  }
  $count++;
}
close IN;

# produce output file from second fastq
my $print = 0;
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne "@");
  my @spl = split(" ", $line);
  my $q = <FQ>;
  my $w = <FQ>;
  my $e = <FQ>;
  if (exists $un{$spl[0]}) {
    print OUT "$line$q$w$e";
    $print++;
  }
}
close FQ;
close OUT;
