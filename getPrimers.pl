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
                 The chromosome names (first column) must match the
                   first space-delimited token in the headers of the
                   fasta reference genome.
    <genome>   Fasta file of reference genome
    <outfile>  Output file containing primer and target sequences
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

open(BED, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(GEN, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]") || die "Cannot open $ARGV[2] for writing\n";

# load primer locations: chr# - 5'Loc - 5'PLen - targLen - 3'PLen
my %pos;
my $count = 0;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    print STDERR "Warning! Improperly formatted line in $ARGV[0]: $line\n  ",
      "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $pos{$spl[3]}) {
    my @div = split("\t", $pos{$spl[3]});
    if ($div[0] ne $spl[0]) {
      print STDERR "Warning! Skipping amplicon $spl[3] --\n",
        "  located at chromosomes $spl[0] and $div[0]!?\n";
      delete $pos{$spl[3]};
    } elsif ($spl[1] < $div[1]) {
      $pos{$spl[3]} = "$spl[0]\t$spl[1]\t".($spl[2]-$spl[1]).
        "\t".($div[1]-$spl[2])."\t$div[2]";
      $count++;
    } else {
      $pos{$spl[3]} .= "\t".($spl[1]-$div[1]-$div[2]).
        "\t".($spl[2]-$spl[1]);
      $count++;
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
    print STDERR "Warning! Skipping amplicon $amp --\n",
      "  Not enough info in BED file\n";
    next;
  }
  $loc{$spl[0]}{$amp} = "$spl[1]\t$spl[2]\t$spl[3]\t$spl[4]";
}

# retrieve primers from genome
$/ = '>';
my %rep_ch;  # to check for repeated chromosome names
my $total = 0;
my $waste = <GEN>;
while (my $chunk = <GEN>) {
  chomp $chunk;
  my @spl = split("\n", $chunk);
  my @head = split(" ", shift @spl);
  my $ch = $head[0];
  if (exists $rep_ch{$ch}) {
    die "Error! In reference genome $ARGV[1]:\n" .
      "  Chromosome name $ch repeated\n";
  }
  $rep_ch{$ch} = 1;
  my $chr = join("", @spl);
  foreach my $amp (sort keys %{$loc{$ch}}) {
    my @div = split("\t", $loc{$ch}{$amp});
    if ($div[0] + $div[1] + $div[2] + $div[3] > length $chr) {
      print STDERR "Warning! Skipping amplicon $amp --\n",
        "  Outside bounds of chromosome $ch\n";
      next;
    }
    my $seg = substr($chr, $div[0], $div[1]+$div[2]+$div[3]);
    $seg =~ tr/a-z/A-Z/;
    print OUT "$amp,", substr($seg, 0, $div[1]),
      ",", substr($seg, $div[1]+$div[2], $div[3]),
      ",", substr($seg, $div[1], $div[2]),
      "\n";
    $total++;
  }
}
close GEN;
close OUT;

if ($total < $count) {
  print STDERR "Warning! Could not find all amplicons\n",
    "  Amplicons specified: $count\n",
    "  Sequences found: $total\n";
}
