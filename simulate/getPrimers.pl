#!/usr/bin/perl

# John M. Gaspar
# Nov. 2014

# Retrieve primer and target sequences from a BED file
#   that lists primer positions.

use strict;
use warnings;

sub usage {
  print q(Usage: perl getPrimers.pl  <infile>  <genome>  <outfile>
  Required:
    <infile>   BED file listing locations of primers, tab-delimited.
                 For example:
                   chr7   127413326   127413347   amplicon1
                   chr7   127413430   127413448   amplicon1
                   chr9   89010943    89010965    amplicon2
                   chr9   89011062    89011085    amplicon2
                 Two primers are required for each amplicon.
    <genome>   Fasta file of reference genome
    <outfile>  Output file containing primer and target sequences
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

open(BED, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(GEN, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]");

# load primer locations: chr# - 5'Loc - 5'PLen - targLen - 3'PLen
my %pos;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    print "Warning! Improperly formatted line in $ARGV[0]: $line\n  ",
      "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $pos{$spl[3]}) {
    my @div = split("\t", $pos{$spl[3]});
    if ($div[0] ne $spl[0]) {
      print "Warning: skipping amplicon $spl[3] -- ",
        "located at chromosomes $spl[0] and $div[0]!?\n";
      delete $pos{$spl[3]};
    } elsif ($spl[1] < $div[1]) {
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

# load amplicon locations -- sort by chromosome for efficient retrieval
my %loc;
foreach my $amp (keys %pos) {
  my @spl = split("\t", $pos{$amp});
  if (scalar @spl < 5) {
    print "Warning: skipping amplicon $amp -- ",
      "not enough info in BED file\n";
    next;
  }
  $loc{$spl[0]}{$amp} = "$spl[1]\t$spl[2]\t$spl[3]\t$spl[4]";
}

# retrieve primers from genome
$/ = '>';
my $waste = <GEN>;
while (my $chunk = <GEN>) {
  chomp $chunk;
  my @spl = split("\n", $chunk);
  my $ch = shift @spl;
  my $chr = join("", @spl);
  foreach my $amp (sort keys %{$loc{$ch}}) {
    my @div = split("\t", $loc{$ch}{$amp});
    my $seg = substr($chr, $div[0], $div[1]+$div[2]+$div[3]);
    $seg =~ tr/a-z/A-Z/;
    print OUT "$amp,", substr($seg, 0, $div[1]),
      ",", substr($seg, $div[1]+$div[2], $div[3]),
      ",", substr($seg, $div[1], $div[2]),
      "\n";
  }
}
close GEN;
close OUT;
