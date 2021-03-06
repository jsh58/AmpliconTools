#!/usr/bin/perl

# John M. Gaspar
# Dec. 2014

# Classify alignments of simulated reads based on their actual origins
#   into one of four categories:
#     - correct location
#     - incorrect location
#     - correct and at least one other location (equivalency defined
#         as larger edit distance and worse alignment score)
#     - unmapped


use strict;
use warnings;

sub usage {
  print q(Usage: perl countAltMappingSim.pl  <infile1>  <infile2>  <outfile>
  Required:
    <infile1>  File containing input reads in fastq format, with primers removed
                 and amplicon identification in header (produced by removePrimer)
    <infile2>  SAM file containing mapping information
    <outfile>  Output file listing counts of read alignments
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

open(FQ, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(SAM, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]") || die "Cannot open $ARGV[2] for writing\n";

# load amplicon info for reads from fastq
my %loc;
my $count = 0;
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne '@');
  chomp $line;
  $count++;
  my @spl = split(" ", $line);
  my $que = substr($spl[0], 1);
  die "$ARGV[0] is not properly formatted\n" if (scalar @spl < 4);
  my $ref = "";
  my $pl = 0;
  my $am = "";
  for (my $x = 1; $x < scalar @spl; $x++) {
    my @div = split('=', $spl[$x]);
    $ref = $div[1] if ($div[0] eq "ref");
    $pl = $div[1] if ($div[0] eq "pos");
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $am = $spl[$x-1];
      last;
    }
  }
  die "Unknown amplicon in read header $line\n" if ((! $pl) || (! $ref));
  $loc{$que} = "$am\t$ref\t$pl";
  for (my $x = 0; $x < 3; $x++) {
    $line = <FQ>;
  }
}
close FQ;
print "Reads: $count\n";

# parse SAM, analyze alternative mapping locations
$count = 0;
my $mapped = 0;
my $skip = 0;
my $un = 0;
my %tot;  # number of reads for an amplicon
my %yes;  # mapped to correct origin
my %no;   # mapped to incorrect origin
my %both; # mapped to correct + incorrect
my %unm;  # unmapped
my %res;  # saves previous result with this read
my $line = <SAM>;
while ($line) {
  if (substr($line, 0, 1) eq '@') {
    $line = <SAM>;
    next;
  }
  my @spl = split("\t", $line);
  $count++;

  # skip if no amplicon info for read
  if (! exists $loc{$spl[0]}) {
    $skip++;
    #print "Warning: skipping read $spl[0] (no amplicon info)\n";
    while ($line = <SAM>) {
      my @div = split("\t", $line);
      last if ($div[0] ne $spl[0]);
    }
    next;
  }
  # load amplicon info for read
  my $ramp = $loc{$spl[0]};
  my @cut = split("\t", $loc{$spl[0]});  # expected location: chr $cut[0],
                                         #  min pos $cut[1], max pos $cut[2]
  shift @cut;
  $cut[1] -= 99;
  $cut[2] = $cut[1] + 199;

  my $prev = 0;  # previous result for this read
  if (! exists $res{$spl[0]}) {
    $tot{$ramp}++;
  } else {
    die "Too many results for $spl[0]\n" if ($res{$spl[0]} > 4);
    $prev = $res{$spl[0]};  # save previous result
    $res{$spl[0]} = 5;  # will throw error if analyzed a third time
  }

  # skip if unmapped
  if ($spl[1] & 0x4) {
    my $val = 1;
    if ($prev) {
      $val = 0.5;
      # decrease previous result by 0.5
      if ($prev == 1) {
        $yes{$ramp} -= 0.5;
      } elsif ($prev == 2) {
        $no{$ramp} -= 0.5;
      } elsif ($prev == 3) {
        $both{$ramp} -= 0.5;
      } elsif ($prev == 4) {
        $un -= 0.5;
        $unm{$ramp} -= 0.5;
      } else {
        die "Cannot handle previous result $prev with read $spl[0]\n";
      }
    } else {
      $res{$spl[0]} = 4;
    }
    $un += $val;
    $unm{$ramp} += $val;

    while ($line = <SAM>) {
      my @div = split("\t", $line);
      last if ($div[0] ne $spl[0]);
    }
    next;
  }

  # get data from primary map
  my $nm = getTag("NM", @spl);
  die "Cannot find NM in $line\n" if ($nm == 1000);
  my $as = getTag("AS", @spl);
  die "Cannot find AS in $line\n" if ($as == 1000);
  my $op = getTag("OP", @spl);

  # check read mapping(s) -- only those equiv. to primary
  my $hit = ($op == 1000 ? 0 : 1);  # if realigned, label as expected
  while ($line) {
    chomp $line;
    my @div = split("\t", $line);
    last if ($div[0] ne $spl[0]);
    die "Supp. alignment: $line\n" if ($div[1] & 0x800);
    my $nm2 = getTag("NM", @div);
    die "Cannot find NM in $line\n" if ($nm2 == 1000);
    my $as2 = getTag("AS", @div);
    die "Cannot find AS in $line\n" if ($as2 == 1000);

    # skip if edit distance (NM) is greater than primary map
    #   and alignment score (AS) is worse
    #   -or- primary map was a realignment (OP)
    if (($nm2 > $nm && $as2 < $as) || ($op != 1000)) {
      $line = <SAM>;
      next;
    }

    # determine if read maps to expected location
    my $pos = $div[3];
    if ($div[2] ne $cut[0] || $pos < $cut[1] || $pos > $cut[2]) {
      $hit = ($hit == 1 || $hit == 3 ? 3 : 2);
    } else {
      $hit = ($hit == 2 || $hit == 3 ? 3 : 1);
    }

    $line = <SAM>;
  }

  # record result
  my $val = 1;  # value of a match (decreased to 0.5 in case of duplicate)
  if ($prev) {
    $val = 0.5;
    # decrease previous result by 0.5
    if ($prev == 1) {
      $yes{$ramp} -= 0.5;
    } elsif ($prev == 2) {
      $no{$ramp} -= 0.5;
    } elsif ($prev == 3) {
      $both{$ramp} -= 0.5;
    } elsif ($prev == 4) {
      $un -= 0.5;
      $unm{$ramp} -= 0.5;
    } else {
      die "Cannot handle previous result $prev with read $spl[0]\n";
    }
  } else {
    $res{$spl[0]} = $hit;
  }

  if ($hit == 1) {
    $yes{$ramp} += $val;
  } elsif ($hit == 2) {
    $no{$ramp} += $val;
  } elsif ($hit == 3) {
    $both{$ramp} += $val;
  } else {
    die "Error: No hits for read $spl[0]\n";
  }

}
close SAM;
print "Total reads: $count",
  "\nUnmapped: $un",
  "\nSkipped: $skip",
  "\n";

# print output
my @total;
for (my $x = 0; $x < 5; $x++) {
  $total[$x] = 0;
}
print OUT "#Amplicon\tCorrect\tBoth\t",
  "Incorrect\tUnmapped\tTotal\n";
foreach my $am (sort {$a cmp $b} keys %tot) {
  print OUT "$am";
  if (exists $yes{$am}) {
    print OUT "\t$yes{$am}";
    $total[0] += $yes{$am};
  } else {
    print OUT "\t0";
  } 
  if (exists $both{$am}) {
    print OUT "\t$both{$am}";
    $total[1] += $both{$am};
  } else {
    print OUT "\t0";
  }
  if (exists $no{$am}) {
    print OUT "\t$no{$am}";
    $total[2] += $no{$am};
  } else {
    print OUT "\t0";
  }
  if (exists $unm{$am}) {
    print OUT "\t$unm{$am}";
    $total[3] += $unm{$am};
  } else {
    print OUT "\t0";
  }
  if (exists $tot{$am}) {
    print OUT "\t$tot{$am}";
    $total[4] += $tot{$am};
  } else {
    print OUT "\t0";
  }
  print OUT "\n";
}
print OUT "Total";
for (my $x = 0; $x < 5; $x++) {
  print OUT "\t$total[$x]";
}
print OUT "\n";
close OUT;

# retrieve value from optional field of SAM
sub getTag {
  my $tag = shift @_;
  my @spl = @_;
  my $ret = 1000;
  my $x;
  for ($x = 11; $x < scalar @spl; $x++) {
    my @cut = split(':', $spl[$x]);
    if ($cut[0] eq $tag) {
      $ret = $cut[$#cut];
      last;
    }
  }
  return $ret;
}
