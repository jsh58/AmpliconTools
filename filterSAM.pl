#!/usr/bin/perl

# John M. Gaspar
# Dec. 2014

# Filter a SAM file, taking into account multi-mapping
#   reads and realignments of length-variant reads.

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterSAM.pl  <infile1>  <infile2>  <infile3>  <infile4>  \
                <outfile>  <infile5>  <logfile>  <score>
  Required:
    <infile1>  File containing input reads in fastq format, with primers removed
                 and amplicon identification in header (produced by removePrimer)
    <infile2>  BED file listing locations of primers
    <infile3>  File listing alternative mapping locations and whether putative primers
                 are exact/close (1) or way off (0) (produced by checkAltMapping.pl)
    <infile4>  Input SAM file containing mapping information
    <outfile>  Output SAM file
  Optional:
    <infile5>  File containing realignment info (produced by alignLengthVars.pl)
    <logfile>  Log file for realignments
    <score>    Maximum primer matching score of external in/del to consider
                 (scale 0-1; def. 0.75)
);
  exit;
}

usage() if (scalar @ARGV < 5 || $ARGV[0] eq "-h");

# open files
open(FQ, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(BED, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(ALT, $ARGV[2]) || die "Cannot open $ARGV[2]\n";
open(SAM, $ARGV[3]) || die "Cannot open $ARGV[3]\n";
open(OUT, ">$ARGV[4]");
if (scalar @ARGV > 5) {
  open(ALN, $ARGV[5]) || die "Cannot open $ARGV[5]\n";
  open(LOG, ">$ARGV[6]") if (scalar @ARGV > 6);
}

my $maxExt = 0.75;  # maximum match of external in/del to consider --
                    #   above this value, the primer is a good enough
                    #   match to explain the length variant as an
                    #   artifact
if (scalar @ARGV > 7) {
  die "Maximum primer matching score must be in [0,1]\n"
    if ($ARGV[7] < 0 || $ARGV[7] > 1);
  $maxExt = $ARGV[7];
}

# load expected amplicon locations (convert to 1-based;
#   allow 20bp "wiggle room" to both sides of amplicon)
my %pos;
my %loc;
while (my $line = <BED>) {
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    print "Warning! Improperly formatted line in $ARGV[1]: $line\n  ",
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

# load amplicon info for reads from fastq
my %amp;  # saves amplicon and (one|both)
my %dup;  # saves amplicon and (one|both) for reads listed
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
  die "Error! $ARGV[0] is improperly formatted:\n  ",
    "no amplicon ID in $line\n" if (!$id);
  my $que = substr($spl[0], 1);

  # check amplicon, and for duplicates
  foreach my $am (keys %loc) {
    if ($id eq $am) {
      if (exists $amp{$que}) {
        my @div = split("\t", $amp{$que});
        if ($id ne $div[0]) {
          print "Warning! Primers do not match for read $que\n";
          delete $amp{$que};
          last;
        }
        # save previous to %dup
        $dup{$que}{$seq{$que}} = $amp{$que};
        chomp ($line = <FQ>);
        $dup{$que}{$line} = "$am\t".($spl[$#spl] eq "both" ? "both" : "one");
        $cdup++;
      } else {
        $amp{$que} = "$am\t".($spl[$#spl] eq "both" ? "both" : "one");
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

# load alternative site information
my %alt;
while (my $line = <ALT>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split("\t", $line);
  my $res = pop @spl;
  die "Error! $ARGV[2] is improperly formatted\n"
    if (scalar @spl < 5 || (($res != 0) && ($res != 1))); 
  my @div = split('-', $spl[$#spl]);
  for (my $x = $div[0]; $x < $div[$#div]+1; $x++) {
    $alt{"$spl[0]\t$spl[1]"}{"$spl[2]\t$spl[3]\t$x"} = $res;
  }
}
close ALT;

# load realignments
my %aln;  # for new alignments
my %tot; my %xExt; my %xInv; my %xWorse; my %xSame; my %xFil;
  my %yBetter; my %yUnmap;  # counting variables
if (scalar @ARGV > 5) {
  my $total = 0;
  while (my $line = <ALN>) {
    chomp $line;
    next if (substr($line, 0, 1) eq '#');
    $total++;
    my @spl = split("\t", $line);
    die "Error! $ARGV[5] is improperly formatted\n" if (scalar @spl < 10);
    $tot{"$spl[1] $spl[3]"}++;
    if (($spl[$#spl-1] eq "external") && ($spl[$#spl] > $maxExt)) {
      $xExt{"$spl[1] $spl[3]"}++;  # external artifact
      next;
    }

    # determine if new alignment is valid --
    #   scoring function similar to bowtie2's default,
    #   but considering only subs ($spl[2])
    my $len = length $spl[7];  # length of ref.
    $spl[3] =~ m/^(\d+)([ID])/;
    $len -= $1 if ($2 eq 'I');
    if ($spl[2] > (0.6*$len+0.6) / 5) {
      $xInv{"$spl[1] $spl[3]"}++;  # invalid new alignment
      next;
    }

    my $read = shift @spl;
    $aln{$read} = join("\t", @spl);
  }
  close ALN;
  #print "Total realignments: $total\n",
  #  "Valid realignments: ", scalar keys %aln, "\n";
}

# parse SAM, produce filtered file
$count = 0; my $mapped = 0;
my $ct = 0; my $un = 0;
my $real = 0; my $pral = 0;
my $line = <SAM>;
while ($line) {
  if (substr($line, 0, 1) eq '@') {
    print OUT $line;  # keep header in output
    $line = <SAM>;
    next;
  }
  chomp $line;
  my @spl = split("\t", $line);
  die "Error! $ARGV[3] is improperly formatted\n" if (scalar @spl < 11);
  die "Error! $ARGV[3] contains a supplementary alignment:\n$line\n"
    if ($spl[1] & 0x800);
  die "Error! $ARGV[3] is improperly formatted ",
    "(possibly coordinate sorted?)\n$line\n" if ($spl[1] & 0x100);
  $count++;

  # prepare to check mapping(s)
  my @res = ();                # array for saving good mapping(s)
  my $idx = -1;                # index of correct map
  my $map = 0;                 # primary map to correct location

  if ($spl[1] & 0x4) {
    # unmapped: check only for new alignment
    $un++;
    while ($line = <SAM>) {
      my @div = split("\t", $line);
      last if (scalar @div < 11 || $div[0] ne $spl[0] ||
        ($div[9] ne $spl[9] && $div[9] ne revComp($spl[9])
        && $div[9] ne '*'));
    }
  } elsif (! exists $amp{$spl[0]}) {
    # no amplicon info: save results, check for new alignment
    print "Warning: no amplicon info for read $spl[0]\n";
    push @res, $line;
    while ($line = <SAM>) {
      chomp $line;
      my @div = split("\t", $line);
      last if (scalar @div < 11 || $div[0] ne $spl[0] ||
        ($div[9] ne $spl[9] && $div[9] ne revComp($spl[9])
        && $div[9] ne '*'));
      push @res, $line;
    }
  } else {
    # load amplicon info for read
    my @brk = split("\t", $amp{$spl[0]});
    my @cut = split("\t", $loc{$brk[0]});  # expected location: chr $cut[0],
                                           #   min pos $cut[1], max pos $cut[2]

    # check mapping(s)
    my $as = getTag("AS", @spl); # get alignment score
    while ($line) {
      chomp $line;
      my @div = split("\t", $line);
      last if (scalar @div < 11 || $div[0] ne $spl[0] ||
        ($div[9] ne $spl[9] && $div[9] ne revComp($spl[9])
        && $div[9] ne '*'));
      die "Error! $ARGV[3] contains a supplementary alignment:\n$line\n"
        if ($div[1] & 0x800);
      my $xs = getTag("AS", @div);

      # get location, accounting for in/dels
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
      my $rc = ($div[1] & 0x10 ? 1 : 0);  # read maps to minus strand
      my $pos = ($rc ? $div[3] + $off + length $div[9] : $div[3]);

      # mapped to correct location?
      if (($div[2] eq $cut[0]) && ($pos >= $cut[1]) && ($pos <= $cut[2])) {
        push @res, $line;            # save correct mapping
        if ($idx == -1) {
          $idx = $#res;              # save index (if not already saved)
          $map = 1 if ($xs == $as);  # equivalent to primary
        }
      } else {
        # not mapped to correct location: determine if it is legit

        # load removed-primer info
        my $pr = "";  # for primer name and (both|one)
        if (exists $dup{$div[0]}) {
          my $test = ($rc ? revComp($div[9]) : $div[9]);
          foreach my $seq (keys %{$dup{$div[0]}}) {
            if ($test eq $seq) {
              $pr = $dup{$div[0]}{$seq};
              last;
            }
          }
        } else {
          $pr = $amp{$div[0]};
        }
        if (!$pr) {
          print "Warning! No removed-primer info for read $div[0]\n";
          #push @res, $line;  # save mapping(?)
          $line = <SAM>;
          next;
        }

        # check alternative location
        my $lc = ($rc ? '-' : '+')."\t$div[2]\t$pos"; # 2nd key to %alt
        if ((! exists $alt{$pr}) || (! exists $alt{$pr}{$lc})) {
          print "Warning! No alt. location info for read $div[0], $pr\t$lc\n";
          #push @res, $line;  # save mapping(?)
          $line = <SAM>;
          next;
        }
        if ($alt{$pr}{$lc}) {
          # a match: save the result
          push @res, $line;
        }
      }

      $line = <SAM>;
    }
  }

  # check for new alignment
  if (exists $aln{$spl[0]}) {

    my @div = split("\t", $aln{$spl[0]});
    my $seq = ($spl[1] & 0x10 ? revComp($spl[9]) : $spl[9]);

    # check sequence -- allow 1bp diff in case of external in/del
    if ($div[6] eq $seq || substr($div[6], 1) eq $seq ||
        substr($div[6], 0, -1) eq $seq) {
      delete $aln{$spl[0]};

      # construct new SAM record
      my $mapq = 255;       # mapping quality (255 = "unavailable")
      my $xm = $div[1];     # number of subs
      my $as = -5*$xm - 1;  # for AS, count subs as -5 each, in/del as -1 only
      my $rec = "$spl[0]\t0\t$div[3]\t$div[4]\t$mapq\t$div[5]".
        "\t*\t0\t0\t$div[6]\t$div[7]\tAS:i:$as";
      $div[2] =~ m/^(\d+)([ID])/;
      my $xg = $1;          # length of gap
      my $xn = 0;           # assume no Ns in ref
      my $xo = 1;           # 1 gap open, by definition
      my $nm = $xm + $xg;   # "edit distance"
      my $rec2 = "\tXN:i:$xn\tXM:i:$xm\tXO:i:$xo\tXG:i:$xg".
        "\tNM:i:$nm\tMD:Z:$div[8]";

      # compare new alignment to previous (if any)
      #   (consider only number of subs [XM])
      my $x;  # index to insert new alignment into @res
      for ($x = 0; $x < scalar @res; $x++) {
        my @fun = split("\t", $res[$x]);
        my $oldXM = getTag("XM", @fun);  # number of subs
        if ($xm < $oldXM || (($div[$#div-1] eq "external")
            && ($xm == $oldXM))) {
          last;
        }
      }

      my $skip = 0;
      if ($idx != -1) {
        # compare to previous correct map
        my @fun = split("\t", $res[$idx]);
        if ($x <= $idx) {
          # save previous alignment info
          $rec2 .= "\tOP:i:$fun[3]\tOC:Z:$fun[5]";
          $idx = $x;
          $yBetter{"$div[0] $div[2]"}++; # new alignment is better
        } else {
          # compare CIGARs
          my $d = 0; my $i = 0;
          while ($fun[5] =~ m/(\d+)([ID])/g) {
            $d += $1 if ($2 eq 'D');
            $i += $1 if ($2 eq 'I');
          }
          while ($div[5] =~ m/(\d+)([ID])/g) {
            $d -= $1 if ($2 eq 'D');
            $i -= $1 if ($2 eq 'I');
          }
          if (!$d && !$i) {
            $xSame{"$div[0] $div[2]"}++;  # equiv. I/D
            $skip = 1 if ($fun[5] eq $div[5]);  # skip identical realignment
          } else {
            $xWorse{"$div[0] $div[2]"}++;  # prev. alignment is better (or equal)
          }
        }
      } else {
        $yUnmap{"$div[0] $div[2]"}++;  # previously unmapped to correct loc.
      }

      $rec2 .= "\tYT:Z:UU";
      if (@res && !$x) {
        my $oldAS = getTag("AS", split("\t", $res[$x]));
        $rec .= "\tXS:i:$oldAS";  # save new XS
      }

      # add new SAM record to @res
      if (! $skip) {
        splice(@res, $x, 0, $rec.$rec2);
        $real++;
        $pral++ if (!$x);
      }
    } else {
      print "Warning! Cannot consider realignment for read $spl[0]\n",
        "Sequences do not match:\n$div[6]\n$seq\n";
    }
  }

  # print output
  if (!@res) {
    # no good maps -- print as unmapped
    $spl[1] = ($spl[1] | 0x4);
    print OUT join("\t", @spl), "\n";
    $ct++;
  } else {

    # promote correct map that is equiv. to primary
    if ($map && $idx > 0) {
      my $temp = $res[0];
      $res[0] = $res[$idx];
      $res[$idx] = $temp;
    }

    # print output(s)
    for (my $x = 0; $x < scalar @res; $x++) {
      my @div = split("\t", $res[$x]);
      # make sure first result is primary, others are secondary:
      $div[1] = ($x ? $div[1] | 0x100 : $div[1] & 0xEFF);
      print OUT join("\t", @div), "\n";
    }
    $mapped++;
  }
}
close SAM;
close OUT;

# print realignment log information
if (scalar @ARGV > 6) {
  print LOG "Amplicon\tIn/Del\tBetter\tNew\tExternalInDel\t",
    "Invalid\tWorse\tSame/Equiv\tFiltered\tTotal\n";

  # count realignments not used
  foreach my $x (keys %aln) {
    my @div = split("\t", $aln{$x});
    $xFil{"$div[0] $div[2]"}++;
  }

  foreach my $x (sort keys %tot) {
    my @div = split(" ", $x);
    print LOG "$div[0]\t$div[1]\t",
      (exists $yBetter{$x} ? $yBetter{$x} : 0), "\t",
      (exists $yUnmap{$x} ? $yUnmap{$x} : 0), "\t",
      (exists $xExt{$x} ? $xExt{$x} : 0), "\t",
      (exists $xInv{$x} ? $xInv{$x} : 0), "\t",
      (exists $xWorse{$x} ? $xWorse{$x} : 0), "\t",
      (exists $xSame{$x} ? $xSame{$x} : 0), "\t",
      (exists $xFil{$x} ? $xFil{$x} : 0), "\t",
      (exists $tot{$x} ? $tot{$x} : 0), "\n";
  }
  close LOG;
}

#print "Reads analyzed: $count\n",
#  "  (Previously unmapped: $un)\n",
#  "Post-filtering unmapped: $ct\n",
#  "Post-filtering mapped: $mapped\n";
#print "Total realignments: $real\n",
#  "  Primary maps: $pral\n" if (scalar @ARGV > 5);

# reverse-complement a sequence
sub revComp {
  my $seq = $_[0];
  $seq = reverse $seq;
  $seq =~ tr/ACGT/TGCA/;
  return $seq;
}

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
  if ($ret == 1000) {
    print "Error! Cannot find $tag in SAM record:\n",
      join("\t", @spl), "\n";
    die;
  }
  return $ret;
}
