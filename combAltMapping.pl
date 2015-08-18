#!/usr/bin/perl

# John M. Gaspar
# Dec. 2014

# Combine outputs from checkAltMapping.pl run on
#   multiple samples.

use strict;
use warnings;

sub usage {
  print q(Usage: perl combAltMapping.pl  <outfile>  <infile1>  <infile2>  [...]
  Required:
    <outfile>  Output file listing combined match judgments
    <infile1>  Output from checkAltMapping.pl for sample #1
    <infile2>    "     "         "             "  sample #2
                 (etc. for multiple samples)
);
  exit;
}

usage() if (scalar @ARGV < 2 || $ARGV[0] eq "-h");

open(OUT, ">$ARGV[0]");

# loop through folders
my %tot;
my $samp = 0;
for (my $x = 1; $x < scalar @ARGV; $x++) {
  if ( ! open(IN, $ARGV[$x]) ) {
    print "Error! Cannot open $ARGV[$x]\n";
    next;
  }
  $samp++;

  # load checkAltMapping results -- 
  #   disagreements favor 1 (match)
  while (my $line = <IN>) {
    next if (substr($line, 0, 1) eq '#');
    chomp $line;
    my @spl = split("\t", $line);
    my @div = split('-', $spl[4]);
    for (my $x = $div[0]; $x < $div[$#div]+1; $x++) {
      if (exists $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x}) {
        if ($spl[5] != $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x}) {
          #print "Warning: $dr has $spl[5] for $spl[0],$spl[1],$spl[2],$spl[3],$x,",
          #  " previously $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x}\n";
          $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x} = 1;
        }
      } else {
        $tot{$spl[0]}{$spl[1]}{$spl[2]}{$spl[3]}{$x} = $spl[5];
      }
    }
  }
  close IN;
}

# print output
print OUT "#Amplicon\tPrimersRemoved\tStrand\tChrom\tPosition(s)\tMatch?\n";
foreach my $am (sort {$a cmp $b} keys %tot) {
  foreach my $bo (sort {$b cmp $a} keys %{$tot{$am}}) {
    foreach my $st (sort {$a cmp $b} keys %{$tot{$am}{$bo}}) {
      foreach my $ch (sort {$a cmp $b} keys %{$tot{$am}{$bo}{$st}}) {

        my $res; my $min;
        my $prev = -10;  # previous position analyzed
        my @loc = sort {$a <=> $b} keys %{$tot{$am}{$bo}{$st}{$ch}};
        for (my $x = 0; $x < scalar @loc; $x++) {
          # combine results for neighboring positions --
          #   consider a match if any is a match
          if ($loc[$x] < $prev + 4) {
            $res = 1 if ($tot{$am}{$bo}{$st}{$ch}{$loc[$x]});
          } else {
            if ($x) {
              printf OUT "$am\t$bo\t$st\t$ch\t%s\t$res\n",
                ($min != $prev ? "$min-$prev" : $min);
            }
            $res = $tot{$am}{$bo}{$st}{$ch}{$loc[$x]};
            $min = $loc[$x];
          }
          $prev = $loc[$x];
        }
        if (@loc) {
          printf OUT "$am\t$bo\t$st\t$ch\t%s\t$res\n",
            ($min != $prev ? "$min-$prev" : $min);
        }
      }
    }
  }
}
close OUT;
#print "Samples analyzed: $samp\n";
