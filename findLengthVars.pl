#!/usr/bin/perl

# John M. Gaspar
# Nov. 2014

# Find putative small structural variants (based on
#   length variants) in a fastq file.

use strict;
use warnings;

sub usage {
  print q(Usage: perl findLengthVars.pl  <infile1>  <infile2>  <outfile>  <percent>  <distance>
  Required:
    <infile1>  File containing input reads in fastq format, with primers removed
                 and amplicon identification in header (produced by removePrimer;
                 may be gzip compressed, with ".gz" extension)
    <infile2>  BED file listing locations of primers
    <outfile>  Output file for length variants
  Optional:
    <percent>  Fraction of reads to consider, e.g. 0.01 means that 1% of reads from
                 a particular amplicon, such as STK11_t1_2, must be significantly
                 different in length from 79bp, in order to list the variants in
                 <outfile>  (def. 0.01)
    <distance> Minimum length difference that is considered "significant", e.g.
                 a value of 5 means all reads 75-83bp are NOT counted as different
                 from the expected 79bp for STK11_t1_2  (def. 5bp)
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

if (substr($ARGV[0], -3) eq ".gz") {
  die "Cannot open $ARGV[0]\n" if (! -f $ARGV[0]);
  open(IN, "zcat $ARGV[0] |");
} else {
  open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
}
open(BED, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]") || die "Cannot open $ARGV[2] for writing\n";

# load length parameters
my $pct = 0.01;  # fraction of reads to consider
my $dist = 5;    # minimum distance from expected to consider
$pct = $ARGV[3] if (scalar @ARGV > 3);
$dist = $ARGV[4] if (scalar @ARGV > 4);

# load expected lengths
my %pos;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    warn "Warning! Improperly formatted line in $ARGV[1]: $line\n  ",
      "Need chromName, start, end, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $pos{$spl[3]}) {
    my @div = split("\t", $pos{$spl[3]});
    if ($div[0] ne $spl[0]) {
      warn "Warning: skipping amplicon $spl[3] --\n",
        "  located at chromosomes $spl[0] and $div[0]!?\n";
      delete $pos{$spl[3]};
    }
    if ($spl[1] < $div[1]) {
      $pos{$spl[3]} = "$spl[0]\t$spl[1]\t".($spl[2]-$spl[1]).
        "\t".($div[1]-$spl[2])."\t$div[2]";
    } else {
      $pos{$spl[3]} .= "\t".($spl[1]-$div[1]-$div[2]).
        "\t".($spl[2]-$spl[1]);
    }
  } else {
    $pos{$spl[3]} = "$spl[0]\t$spl[1]\t".($spl[2]-$spl[1]);
  }
}
close BED;

# check and reconfigure expected lengths hash
foreach my $k (keys %pos) {
  my @spl = split("\t", $pos{$k});
  if (scalar @spl < 4) {
    warn "Warning! Insufficient information in $ARGV[1] for amplicon $k\n";
    delete $pos{$k};
  } else {
    $pos{$k} = $spl[3];
  }
}

# parse fastq file
my $count = 0;
my %tot;  # for lengths of each amplicon
my %len;  # for length variants
while (my $line = <IN>) {
  next if (substr($line, 0, 1) ne "@");
  chomp $line;

  # determine amplicon ID
  my @spl = split(" ", $line);
  my $id = "";
  for (my $x = 1; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $id = $spl[$x-1];
      last;
    }
  }
  die "Error! $ARGV[0] is improperly formatted:\n  ",
    "no amplicon ID in $line\n" if (!$id);

  if (! exists $pos{$id}) {
    warn "Warning! Skipping read $spl[0] with unknown amplicon $id\n";
    for (my $x = 0; $x < 3; $x++) {
      $line = <IN>;
    }
    next;
  }
  $line = <IN>;
  chomp $line;
  $count++;
  my $ln = length $line;
  $tot{$id}++;
  # record length if it is variant
  $len{$id}{$ln}++ if (abs($ln - $pos{$id}) >= $dist);
  for (my $x = 0; $x < 2; $x++) {
    $line = <IN>;
  }
}
close IN;

#print "Reads analyzed: $count\n";

# print output -- only variants
print OUT "#Amplicon\tExpected\tVarLength\tPercent\n";
foreach my $x (sort keys %len) {
  foreach my $y (sort keys %{$len{$x}}) {
    # only if variant > specified percentage $pct
    my $ratio = $len{$x}{$y} / $tot{$x};
    if ($ratio >= $pct) {
      printf OUT "$x\t$pos{$x}\t$y\t%.1f\n", $ratio*100;
    }
  }
}
close OUT;
