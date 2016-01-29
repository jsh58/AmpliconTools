#!/usr/bin/perl

# John M. Gaspar
# Dec. 2014

# Check read mapping to locations other than expected.

use strict;
use warnings;

sub usage {
  print q(Usage: perl checkAltMapping.pl  <infile1>  <infile2>  <infile3>  <infile4>  \
                      <genome>  <outfile>  <score>  <logfile>
  Required:
    <infile1>  File containing input reads in fastq format, with primers removed
                 and amplicon identification in header (produced by removePrimer)
    <infile2>  SAM file containing mapping information
    <infile3>  File listing primer and target sequences (produced by getPrimers.pl)
    <infile4>  BED file listing locations of primers
    <genome>   Fasta file of reference genome
    <outfile>  Output file listing match judgments
  Optional:
    <score>    Minimum primer matching score (scale 0-1; def. 0.75)
    <logfile>  Output file that lists primer alignments and scores
);
  exit;
}

usage() if (scalar @ARGV < 6 || $ARGV[0] eq "-h");

open(FQ, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(SAM, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(PR, $ARGV[2]) || die "Cannot open $ARGV[2]\n";
open(BED, $ARGV[3]) || die "Cannot open $ARGV[3]\n";
open(FA, $ARGV[4]) || die "Cannot open $ARGV[4]\n";
open(OUT, ">$ARGV[5]");
my $pct = 0.75;
if (scalar @ARGV > 6) {
  die "Minimum primer matching score must be in [0,1]\n"
    if ($ARGV[6] < 0 || $ARGV[6] > 1);
  $pct = $ARGV[6];
  if (scalar @ARGV > 7) {
    open(LOG, ">$ARGV[7]");
    print LOG "Amplicon\tPrimersRemoved\tAltChr\tAltPos\tStrand\n";
  }
}

# load primer sequences
my %fwd; my %rev;
while (my $line = <PR>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split(',', $line);
  next if (scalar @spl < 3);
  $fwd{$spl[0]} = $spl[1];
  $rev{$spl[0]} = $spl[2];
}
close PR;

# load expected amplicon locations (convert to 1-based;
#   allow 20bp "wiggle room" to both sides of amplicon)
my %pos;
my %loc;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    print "Warning! Improperly formatted line in $ARGV[3]: $line\n  ",
      "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $pos{$spl[3]}) {
    my @div = split("\t", $pos{$spl[3]});
    if ($div[0] ne $spl[0]) {
      print "Warning: skipping amplicon $spl[3] -- ",
        "located at chromosomes $spl[0] and $div[0]!?\n";
    } else {
      # 20bp wiggle room on either side
      $loc{$spl[3]} = ($spl[1] < $div[1] ?
        "$spl[0]\t".($spl[1]-19)."\t".($div[2]+20) :
        "$spl[0]\t".($div[1]-19)."\t".($spl[2]+20));
    }
    delete $pos{$spl[3]};
  } else {
    $pos{$spl[3]} = "$spl[0]\t$spl[1]\t$spl[2]";
  }
}
close BED;

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

# load amplicon info for reads from fastq
my %amp;  # saves amplicon
my %pri;  # saves removed-primer information: (fwd|rev)(\tboth)?
my %dup;  # saves removed-primer information for reads listed
          #   multiple times (singletons, unjoined)
my %seq;  # saves sequences for reads listed multiple times --
          #   used as 2nd key for %dup
my $count = 0; my $cdup = 0;
$/ = "\n";
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne '@');
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
  die "Error! $ARGV[0] is improperly formatted:\n",
    "  no amplicon ID in $line\n" if (!$id);
  my $que = substr($spl[0], 1);

  # check amplicon, and for duplicates
  foreach my $am (keys %loc) {
    if ($id eq $am) {
      if (exists $amp{$que}) {
        if ($id ne $amp{$que}) {
          print "Warning! Primers do not match for read $que\n";
          delete $amp{$que};
          last;
        }
        if (exists $pri{$que}) {
          # save previous to %dup
          $dup{$que}{$seq{$que}} = $pri{$que};
          delete $pri{$que};
        }
        chomp ($line = <FQ>);
        $dup{$que}{$line} = ($spl[$#spl] eq "both" ?
          "$spl[$#spl-1]\t$spl[$#spl]" : $spl[$#spl]);
        $cdup++;
      } else {
        $amp{$que} = $am;
        $pri{$que} = ($spl[$#spl] eq "both" ?
          "$spl[$#spl-1]\t$spl[$#spl]" : $spl[$#spl]);
        chomp ($line = <FQ>);
        $seq{$que} = $line;  # save sequence in case of duplicate
      }
      last;
    }
  }

  if (! exists $amp{$que}) {
    print "Warning! Skipping read $que -- no amplicon found\n";
    $line = <FQ>;
  } 
  $count++;
  for (my $x = 0; $x < 2; $x++) {
    $line = <FQ>;
  }
}
close FQ;
#print "Reads: $count\n";
#print "Unique: ", scalar keys %amp,
#  "\nDuplicates: $cdup\n";
%seq = ();

# parse SAM, analyze alternative mapping locations
$count = 0; my $mapped = 0;
my $ct = 0; my $un = 0;
my %alt;  # 1 if a match, else 0
my %tot;  # number of reads at an alt. location
my %ln;   # lengths of reads analyzed at an alt. location
my $line = <SAM>;
while ($line) {
  if (substr($line, 0, 1) eq '@') {
    $line = <SAM>;
    next;
  }
  my @spl = split("\t", $line);
  die "Error! $ARGV[1] is improperly formatted\n" if (scalar @spl < 11);
  die "Error! $ARGV[1] contains a supplementary alignment:\n$line\n"
    if ($spl[1] & 0x800);
  die "Error! $ARGV[1] is improperly formatted ",
    "(possibly coordinate sorted?)\n$line\n" if ($spl[1] & 0x100);
  $count++;
  if ($spl[1] & 0x4) {
    $un++;
    while ($line = <SAM>) {
      my @div = split("\t", $line);
      last if (scalar @div < 11 || $div[0] ne $spl[0] ||
        ($div[9] ne $spl[9] && $div[9] ne revComp($spl[9])
        && $div[9] ne '*'));
    }
    next;
  }

  # load amplicon info for read
  if (! exists $amp{$spl[0]}) {
    print "Warning: skipping read $spl[0] (no amplicon info)\n";
    while ($line = <SAM>) {
      my @div = split("\t", $line);
      last if (scalar @div < 11 || $div[0] ne $spl[0] ||
        ($div[9] ne $spl[9] && $div[9] ne revComp($spl[9]) 
        && $div[9] ne '*'));
    }
    next;
  }
  my $ramp = $amp{$spl[0]};
  my @cut = split("\t", $loc{$ramp});  # expected location: chr $cut[0],
                                       #  min pos $cut[1], max pos $cut[2]

  # check read mapping
  while ($line) {
    chomp $line;
    my @div = split("\t", $line);
    last if (scalar @div < 11 || $div[0] ne $spl[0] ||
      ($div[9] ne $spl[9] && $div[9] ne revComp($spl[9]) 
      && $div[9] ne '*'));
    die "Error! $ARGV[1] contains a supplementary alignment:\n$line\n"
      if ($div[1] & 0x800);
    my $rc = ($div[1] & 0x10 ? 1 : 0); # read maps to minus strand

    # load removed-primer info
    my @pr = ();
    if (exists $pri{$div[0]}) {
      @pr = split("\t", $pri{$div[0]});
    } elsif (exists $dup{$div[0]}) {
      my $test = ($rc ? revComp($div[9]) : $div[9]);
      foreach my $seq (keys %{$dup{$div[0]}}) {
        if ($test eq $seq) {
          @pr = split("\t", $dup{$div[0]}{$seq});
          last;
        }
      }
    }
    if (!@pr) {
      print "Warning! No removed-primer info for read $div[0]\n";
      $line = <SAM>;
      next;
    }

    # get offset (D/I)
    my $off = 0;
    while ($div[5] =~ m/(\d+)([ID])/g) {
      if ($2 eq 'D') {
        $off += $1;
      } elsif ($2 eq 'I') {
        $off -= $1;
      }
      # can also throw an error for 'S' or 'H' here
    }

    # save 1st position after fwd primer:
    my $pos = ($rc ? $div[3] + $off + length $div[9] : $div[3]);

    # save length if both primers removed, else 0
    my $len = ($pr[$#pr] eq "both" ? $off + length $div[9] : 0);
    my $both = ($len ? 1 : 0);  # boolean if both primers removed

    # examine only reads that do not map to expected location
    if ($div[2] ne $cut[0] || $pos < $cut[1] || $pos > $cut[2]) {

      # skip if no genomic segment loaded
      if (! exists $chr{$div[2]}) {
        print "Warning! No sequence loaded for reference $div[2]\n";
        $line = <SAM>;
        next;
      }

      # skip if it has been analyzed before
      my $flag = 0;    # skipping flag
      if (exists $alt{$ramp}{$both}{$rc}{$div[2]}{$pos}) {
        if ($alt{$ramp}{$both}{$rc}{$div[2]}{$pos}) {
          # already found to be a match
          $flag = 1;
        } else {
          my @lns = split("\t", $ln{$ramp}{$both}{$rc}{$div[2]}{$pos});
          for (my $x = 0; $x < scalar @lns; $x++) {
            if ($len == $lns[$x]) {
              # length has been examined before
              $flag = 1;
              last;
            }
          }
        }
      }
      # skipping subroutine:
      if ($flag) {
        $tot{$ramp}{$both}{$rc}{$div[2]}{$pos}++;
        $line = <SAM>;
        next;
      }

      # get primer sequences (based on which primers were removed)
      my $fiveP; my $threeP = "";
      if ($pr[0] eq "fwd") {
        $fiveP = $fwd{$ramp};
        $threeP = revComp($rev{$ramp}) if ($both);
      } else {
        $fiveP = revComp($rev{$ramp});
        $threeP = $fwd{$ramp} if ($both);
      }

      # get upstream and downstream sequences
      #   -- allow 1bp wiggle room
      my $fwdP; my $revP = "";
      if ($rc) {
        $fwdP = substr($chr{$div[2]}, $div[3] - 2 + $off + length $div[9], 2 + length $fiveP);
        $fwdP = revComp($fwdP);
        if ($both) {
          $revP = substr($chr{$div[2]}, $div[3] - 2 - length $threeP, 2 + length $threeP);
        }
      } else {
        $fwdP = substr($chr{$div[2]}, $div[3] - 2 - length $fiveP, 2 + length $fiveP);
        if ($both) {
          $revP = substr($chr{$div[2]}, $div[3] - 2 + $off + length $div[9], 2 + length $threeP);
          $revP = revComp($revP);
        }
      }
      $fwdP =~ tr/a-z/A-Z/;
      $revP =~ tr/a-z/A-Z/;

      # score primer-genome alignment
      my ($scoreF, $fPos) = scoreAlign($fwdP, $fiveP);
      my ($scoreR, $rPos) = scoreAlign($revP, $threeP) if ($both);

      # determine result
      my $res = 1;
      $res = 0 if ($scoreF < $pct || ($both && $scoreR < $pct));

      # record score
      $tot{$ramp}{$both}{$rc}{$div[2]}{$pos}++;
      if ($res) {
        # if match, add results to position AND neighbors (3bp)
        for (my $x = -3; $x < 4; $x++) {
          $alt{$ramp}{$both}{$rc}{$div[2]}{$pos+$x} = $res;
        }
        # if both, copy results to other side and singletons as well
        if ($both) {
          my $newRC = ($rc ? 0 : 1);
          my $newPos = ($rc ? $pos-$len : $pos+$len);
          for (my $x = -3; $x < 4; $x++) {
            $alt{$ramp}{$both-1}{$rc}{$div[2]}{$pos+$x} = $res;
            $alt{$ramp}{$both-1}{$newRC}{$div[2]}{$newPos+$x} = $res;
            $alt{$ramp}{$both}{$newRC}{$div[2]}{$newPos+$x} = $res;
          }
        }
      } else {
        $alt{$ramp}{$both}{$rc}{$div[2]}{$pos} = $res;
      }

      # add length analyzed to %ln
      if (exists $ln{$ramp}{$both}{$rc}{$div[2]}{$pos}) {
        $ln{$ramp}{$both}{$rc}{$div[2]}{$pos} .= "\t$len";
      } else {
        $ln{$ramp}{$both}{$rc}{$div[2]}{$pos} = $len;
      }

      # log alignment
      if (scalar @ARGV > 7) {
        printf LOG "$ramp\t%s\t$div[2]\t$pos\t%s\t$div[0]",
          ($both ? "both" : "one"), ($rc ? '-' : '+');
        printf LOG "\nfwd primer $fiveP\ngenome seg %s\t%.3f\n",
          substr($fwdP, $fPos, length $fiveP), $scoreF;
        printf LOG "rev primer $threeP\ngenome seg %s\t%.3f\n",
          substr($revP, $rPos, length $threeP), $scoreR if ($both);
        print LOG "\n";
      }
    }

    $line = <SAM>;
  }

}
close SAM;
close LOG if (scalar @ARGV > 7);
#print "Reads analyzed: $count",
#  "\nUnmapped: $un",
#  "\n";

# print output
print OUT "#Amplicon\tPrimersRemoved\tStrand\tChrom\tPosition(s)\tMatch?\n";
foreach my $am (sort {$a cmp $b} keys %tot) {
  foreach my $bo (sort {$a <=> $b} keys %{$tot{$am}}) {
    foreach my $st (sort {$a <=> $b} keys %{$tot{$am}{$bo}}) {
      foreach my $ch (sort {$a cmp $b} keys %{$tot{$am}{$bo}{$st}}) {

        my $res; my $min;
        my $prev = -10;  # previous position analyzed
        my @loc = sort {$a <=> $b} keys %{$tot{$am}{$bo}{$st}{$ch}};
        for (my $x = 0; $x < scalar @loc; $x++) {
          # combine results for neighboring positions --
          #   consider a match if any is a match
          if ($loc[$x] < $prev + 4) {
            $res = 1 if ($alt{$am}{$bo}{$st}{$ch}{$loc[$x]});
          } else {
            if ($x) {
              printf OUT "$am\t%s\t%s\t$ch\t%s\t$res\n",
                ($bo ? "both" : "one"), ($st ? '-' : '+'),
                ($min != $prev ? "$min-$prev" : $min);
            }
            $res = $alt{$am}{$bo}{$st}{$ch}{$loc[$x]};
            $min = $loc[$x];
          }
          $prev = $loc[$x];
        }
        if (@loc) {
          printf OUT "$am\t%s\t%s\t$ch\t%s\t$res\n",
            ($bo ? "both" : "one"), ($st ? '-' : '+'),
            ($min != $prev ? "$min-$prev" : $min);
        }
      }
    }
  }
}
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
  for (my $y = 0; $y < 3; $y++) {
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
      $score += $val if (substr($ref, $x, 1) eq substr($que, $x+$y, 1));
      $total += $val;
    }
    if ($score/$total > $best) {
      $best = $score/$total;
      $pos = $y;
    }
  }
  return $best, $pos;
}
