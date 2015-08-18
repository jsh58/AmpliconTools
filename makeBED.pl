#!/usr/bin/perl

# John M. Gaspar
# April 2015

# Convert a tabular Fluidigm file describing designed
#   amplicons to a BED file listing primer locations.

use strict;
use warnings;

sub usage {
  print q(Usage: perl makeBED.pl  <infile>  <outfile>
  Required:
    <infile>   Input file describing designed amplicons --
                 Header line should contain fields Assay_Name,
                 F-sp, R-sp, Chr, From, and To.
    <outfile>  Output BED file
);
  exit;
}

usage() if (scalar @ARGV < 2 || $ARGV[0] eq "-h");

# open files
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(OUT, ">$ARGV[1]");

# get fields from header line
my @col;
for (my $x = 0; $x < 6; $x++) {
  $col[$x] = -1;
}
my $line = <IN>;
chomp $line;
my @spl = split("\t", $line);
for (my $x = 0; $x < scalar @spl; $x++) {
  if ($spl[$x] eq "Assay_Name") {
    $col[0] = $x;
  } elsif ($spl[$x] eq "F-sp") {
    $col[1] = $x;
  } elsif ($spl[$x] eq "R-sp") {
    $col[2] = $x;
  } elsif ($spl[$x] eq "Chr") {
    $col[3] = $x;
  } elsif ($spl[$x] eq "From") {
    $col[4] = $x;
  } elsif ($spl[$x] eq "To") {
    $col[5] = $x;
  }
}
for (my $x = 0; $x < 6; $x++) {
  if ($col[$x] == -1) {
    die "Error! Missing info from header in $ARGV[0]\n";
  }
}

# read file
my $total = 0;
my $count = 0;
while (my $line = <IN>) {
  chomp $line;
  $total++;
  my @spl = split("\t", $line);
  if (scalar @spl < 15) {
    #print "Warning! Skipping $line\n";
    next;
  }
  print OUT "$spl[$col[3]]\t", $spl[$col[4]]-1,
    "\t", $spl[$col[4]] - 1 + length $spl[$col[1]],
    "\t$spl[$col[0]]\n",
    "$spl[$col[3]]\t", $spl[$col[5]] - length $spl[$col[2]],
    "\t$spl[$col[5]]\t$spl[$col[0]]\n";
  $count++;
}
close IN;
close OUT;

print "Records in $ARGV[0]: $total",
  "\nValid amplicons: $count",
  "\n";
