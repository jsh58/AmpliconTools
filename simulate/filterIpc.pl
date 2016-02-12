#!/usr/bin/perl

# John M. Gaspar
# Mar. 2015

# Filter output from ipcress.

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterIpc.pl  <infile1>  <infile2>  <genome>  <outfile>  <score>
  Required:
    <infile1>  Output from ipcress
    <infile2>  File listing primer sequences (input to ipcress)
    <genome>   Fasta file of reference genome
    <outfile>  Output file
  Optional:
    <score>    Minimum primer matching score (scale 0-1; def. 0.75)
);
  exit;
}

usage() if (scalar @ARGV < 4 || $ARGV[0] eq "-h");

open(IP, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(PR, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(FA, $ARGV[2]) || die "Cannot open $ARGV[2]\n";
open(OUT, ">$ARGV[3]");
my $pct = 0.75;
if (scalar @ARGV > 4) {
  die "Minimum primer matching score must be in [0,1]\n"
    if ($ARGV[4] < 0 || $ARGV[4] > 1);
  $pct = $ARGV[4];
}

# load primer sequences
my %fwd; my %rev;
while (my $line = <PR>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split("\t", $line);
  next if (scalar @spl < 3);
  $fwd{$spl[0]} = $spl[1];
  $rev{$spl[0]} = $spl[2];
}
close PR;

# load genome
my %chr;
$/ = '>';
my $waste = <FA>;
while (my $chunk = <FA>) {
  chomp $chunk;
  my @spl = split("\n", $chunk);
  my $ch = shift @spl;
  $chr{$ch} = join("", @spl);
}
close FA;

# parse ipcress output
$/ = "\n";
my $count = 0;
while (my $line = <IP>) {
  chomp $line;
  my @spl = split(" ", $line);
  next if ($spl[0] ne "ipcress:");
  next if ($spl[10] eq "single_A" || $spl[10] eq "single_B");

  my @div = split(':', $spl[1]);
  # skip if no genomic segment loaded
  if (! exists $chr{$div[0]}) {
    print "Warning! No sequence loaded for reference $div[0]\n";
    next;
  }

  my $pr5 = ($spl[4] eq 'B' ? $rev{$spl[2]} : $fwd{$spl[2]});
  my $pr3 = ($spl[4] eq 'B' ? $fwd{$spl[2]} : $rev{$spl[2]});
  my $fwdP = substr($chr{$div[0]}, $spl[5], length $pr5);
  my $revP = revComp(substr($chr{$div[0]}, $spl[8], length $pr3));
  $fwdP =~ tr/a-z/A-Z/;
  $revP =~ tr/a-z/A-Z/;

  # score primer-genome alignment
  my $scoreF = scoreAlign($fwdP, $pr5);
  my $scoreR = scoreAlign($revP, $pr3);

  if ($scoreF >= $pct && $scoreR >= $pct) {
    printf OUT "$line %.3f %.3f\n", $scoreF, $scoreR;
  }
}
close IP;
close OUT;

# reverse-complement a sequence
sub revComp {
  my $seq = $_[0];
  $seq = reverse $seq;
  $seq =~ tr/acgtACGT/tgcaTGCA/;
  return $seq;
}

# score alignment -- no in/dels allowed
# weighting function (from 5' end):
#   bases before last 20: 1
#   bases 1-10:           2*pos
#   bases 11-19:          3*pos
#   base 20:              5*pos
sub scoreAlign {
  my $que = $_[0];
  my $ref = $_[1];
  my $len = length $ref;

  my $best = -1;
  my $pos = -1;
  my $score = 0; my $total = 0;
  for (my $x = 0; $x < $len; $x++) {
    my $val = 21 - $len + $x;  # value of match at this position
    if ($val > 19) {
      $val *= 5;
    } elsif ($val > 10) {
      $val *= 3;
    } elsif ($val > 0) {
      $val *= 2;
    } else {
      $val = 1;
    }
    $score += $val if (substr($ref, $x, 1) eq substr($que, $x, 1));
    $total += $val;
  }
  return $score/$total;
}
