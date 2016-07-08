#!/usr/bin/perl

# John M. Gaspar
# Dec. 2014

# Classify read alignments based on their target locations
#   into one of four categories:
#     - target location
#     - alternative location
#     - target and at least one other location (equivalency defined
#         as larger edit distance and worse alignment score)
#     - unmapped
# Version PE: for paired-end read alignments

use strict;
use warnings;

sub usage {
  print q(Usage: perl countAltMappingPE.pl  <infile1>  <infile2>  <infile3>  <outfile>
  Required:
    <infile1>  File containing input reads in fastq format, with primers removed
                 and amplicon identification in header (produced by removePrimer)
    <infile2>  SAM file containing mapping information
    <infile3>  BED file listing locations of primers
    <outfile>  Output file listing counts of read alignments
);
  exit;
}

usage() if (scalar @ARGV < 4 || $ARGV[0] eq "-h");

open(FQ, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(SAM, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(BED, $ARGV[2]) || die "Cannot open $ARGV[2]\n";
open(OUT, ">$ARGV[3]") || die "Cannot open $ARGV[3] for writing\n";

# load expected amplicon locations (convert to 1-based;
#   allow 100bp "wiggle room" to both sides of amplicon)
my %pos;
my %loc;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    warn "Warning! Improperly formatted line in BED file: $line\n",
      "Need chromName, start, end, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $pos{$spl[3]}) {
    my @div = split("\t", $pos{$spl[3]});
    if ($div[0] ne $spl[0]) {
      warn "Warning! Skipping amplicon $spl[3] -- ",
        "located at chromosomes $spl[0] and $div[0]!?\n";
    } else {
      # 100bp wiggle room on either side
      $loc{$spl[3]} = ($spl[1] < $div[1] ?
        "$spl[0]\t".($spl[1]-99)."\t".($div[2]+100) :
        "$spl[0]\t".($div[1]-99)."\t".($spl[2]+100));
    }
    delete $pos{$spl[3]};
  } else {
    $pos{$spl[3]} = "$spl[0]\t$spl[1]\t$spl[2]";
  }
}
close BED;

# load amplicon info for reads from fastq
my %amp;  # saves amplicon
my $count = 0;
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne '@');
  chomp $line;
  $count++;
  my @spl = split(" ", $line);
  my $que = substr($spl[0], 1);
  die "$ARGV[0] is not properly formatted\n" if (scalar @spl < 3);
  die "Error with read $que: amplicons $amp{$que} and $spl[2] don't match\n"
    if ((exists $amp{$que}) && ($amp{$que} ne $spl[2]));
  my $id = "";
  for (my $x = 1; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $id = $spl[$x-1];
      last;
    }
  }
  die "Unknown amplicon in read header $line\n" if (! $id);
  $amp{$que} = $id;
  for (my $x = 0; $x < 3; $x++) {
    $line = <FQ>;
  }
}
close FQ;
print "Reads: $count\n";
print "Unique: ", scalar keys %amp,
  "\n";

# parse SAM, analyze alternative mapping locations
$count = 0;
my $mapped = 0;
my $un = 0;
my $skip = 0;
my %tot;  # number of reads for an amplicon
my %yes;  # mapped to target location
my %no;   # mapped to alt. location
my %both; # mapped to target + alt.
my %unm;  # unmapped
my $line = <SAM>;
while ($line) {
  if (substr($line, 0, 1) eq '@') {
    $line = <SAM>;
    next;
  }
  chomp $line;
  my @spl = split("\t", $line);
  die "Supp. alignment: $line\n" if ($spl[1] & 0x800);
  $count++;

  # skip if no amplicon info for read
  if (! exists $amp{$spl[0]}) {
    $skip++;
    #print "Warning: skipping read $spl[0] (no amplicon info)\n";
    while ($line = <SAM>) {
      my @div = split("\t", $line);
      last if ($div[0] ne $spl[0]);
    }
    next;
  }
  # load amplicon info for read
  my $ramp = $amp{$spl[0]};
  my @cut = split("\t", $loc{$ramp});  # expected location: chr $cut[0],
                                       #  min pos $cut[1], max pos $cut[2]
  $tot{$ramp}++;

  # handle "not properly paired" as separate reads worth 0.5 each
  if (!($spl[1] & 0x2)) {

    for (my $x = 0; $x < 2; $x++) {
      chomp $line;
      my @div = split("\t", $line);
      die "Supp. alignment: $line\n" if ($div[1] & 0x800);

      if ($div[1] & 0x4) {
        $un += 0.5;
        $unm{$ramp} += 0.5;
        $line = <SAM>;
        next;
      }

      # get data from primary map
      my $bit = $div[1] & 0xC0;  # 0x40 if first, 0x80 if second
      my $nm = getTag("NM", @div);
      die "Cannot find NM in $line\n" if ($nm == 1000);
      my $as = getTag("AS", @div);
      die "Cannot find AS in $line\n" if ($as == 1000);

      # check read mapping(s) -- only those equiv. to primary
      my $hit = 0;
      while ($line) {
        chomp $line;
        my @brk = split("\t", $line);
        die "Supp. alignment: $line\n" if ($brk[1] & 0x800);
        last if ($brk[0] ne $spl[0]);
        last if (!($brk[1] & $bit));  # diff. PE read
        my $nm2 = getTag("NM", @brk);
        die "Cannot find NM in $line\n" if ($nm2 == 1000);
        my $as2 = getTag("AS", @brk);
        die "Cannot find AS in $line\n" if ($as2 == 1000);

        # skip if edit distance (NM) is greater than primary map
        #   and alignment score (AS) is worse
        if ($nm2 > $nm && $as2 < $as) {
          $line = <SAM>;
          next;
        }

        # determine if mapping is to expected location
        my $pos = $brk[3];
        if ($brk[2] ne $cut[0] || $pos < $cut[1] || $pos > $cut[2]) {
          $hit = ($hit == 1 || $hit == 3 ? 3 : 2);
        } else {
          $hit = ($hit == 2 || $hit == 3 ? 3 : 1);
        }
        $line = <SAM>;
      }

      # record result
      if ($hit == 1) {
        $yes{$ramp} += 0.5;
      } elsif ($hit == 2) {
        $no{$ramp} += 0.5;
      } elsif ($hit == 3) {
        $both{$ramp} += 0.5;
      } else {
        die "Error: No hits for read $spl[0]\n";
      }
    }
    next;
  }

  # get data from primary map
  my $nm = getTag("NM", @spl);
  die "Cannot find NM in $line\n" if ($nm == 1000);
  my $as = getTag("AS", @spl);
  die "Cannot find AS in $line\n" if ($as == 1000);

  # check read mapping(s) -- only those equiv. to primary
  my $hit = 0;
  my $x = 0;
  while ($line) {
    chomp $line;
    my @div = split("\t", $line);
    die "Supp. alignment: $line\n" if ($div[1] & 0x800);
    last if ($div[0] ne $spl[0]);

    # load PE map
    my $line2 = <SAM>;
    chomp $line2;
    my @d2 = split("\t", $line2);

    # check for problems
    die "Unmapped with $div[0]!?\n"
      if (($div[1] & 0x4) || ($d2[1] & 0x4));
    die "Not properly paired failure with $div[0]\n"
      if ((!($div[1] & 0x2)) || (!($d2[1] & 0x2)));
    die "Paired-end failure with $div[0] and $d2[0]\n"
      if ($div[0] ne $d2[0]);
    die "Supp. alignment: $line\n"
      if ($div[1] & 0x800 || $d2[1] & 0x800);

    # load edit distance, alignment score data
    if ($x == 0) {
      # adjust data for primary map
      $nm += getTag("NM", @d2);
      $as += getTag("AS", @d2);
      $x++;
    }
    my $n1 = getTag("NM", @div);
    my $n2 = getTag("NM", @d2);
    die "Cannot find NM in $line\n" if ($n1 == 1000 || $n2 == 1000);
    my $a1 = getTag("AS", @div);
    my $a2 = getTag("AS", @d2);
    die "Cannot find AS in $line\n" if ($a1 == 1000 || $a2 == 1000);
    my $nm2 = $n1 + $n2;
    my $as2 = $a1 + $a2;

    # skip if edit distance (NM) is greater than primary map
    #   and alignment score (AS) is worse
    if ($nm2 > $nm && $as2 < $as) {
      $line = <SAM>;
      next;
    }

    # check if mapping is to expected location
    my $pos = $div[3];
    my $p2 = $d2[3];
    if (($div[2] ne $cut[0] || $pos < $cut[1] || $pos > $cut[2]) ||
        ($d2[2] ne $cut[0] || $p2 < $cut[1] || $p2 > $cut[2])) {
      $hit = ($hit == 1 || $hit == 3 ? 3 : 2);
    } else {
      $hit = ($hit == 2 || $hit == 3 ? 3 : 1);
    }

    $line = <SAM>;
  }

  # record result
  if ($hit == 1) {
    $yes{$ramp}++;
  } elsif ($hit == 2) {
    $no{$ramp}++;
  } elsif ($hit == 3) {
    $both{$ramp}++;
  } else {
    die "Error: No hits for read $spl[0]\n";
  }
}
close SAM;
print "Reads analyzed: $count",
  "\nUnmapped: $un",
  "\nSkipped: $skip",
  "\n";

# print output
my @total;
for (my $x = 0; $x < 5; $x++) {
  $total[$x] = 0;
}
print OUT "#Amplicon\tTarget\tBoth\t",
  "Alternative\tUnmapped\tTotal\n";
foreach my $am (sort {$a cmp $b} keys %loc) {
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
