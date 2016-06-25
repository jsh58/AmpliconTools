#!/usr/bin/perl

# John M. Gaspar
# Nov. 2014

# Align length variants to genomic segments.

use strict;
use warnings;

sub usage {
  print q(Usage: perl alignLengthVars.pl  <infile1>  <infile2>  <infile3>  <outfile>
                        <logfile>  <genome>
  Required:
    <infile1>  File containing length-variant reads without primers attached
                 (produced by getLengthVars.pl)
    <infile2>  File listing primer and target sequences (produced by getPrimers.pl)
    <infile3>  BED file listing locations of primers
    <outfile>  Output file containing SAM mapping info
  Optional:
    <logfile>  Verbose output file listing possible CIGARs and scores for each read
    <genome>   Fasta file of reference genome (to evaluate external insertions)
);
  exit;
}

usage() if (scalar @ARGV < 4 || $ARGV[0] eq "-h");

# open files
open(FQ, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(PR, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(BED, $ARGV[2]) || die "Cannot open $ARGV[2]\n";
open(OUT, ">$ARGV[3]") || die "Cannot open $ARGV[3] for writing\n";
if (scalar @ARGV > 4) {
  open(LOG, ">$ARGV[4]") || die "Cannot open $ARGV[4] for writing\n"
}
if (scalar @ARGV > 5) {
  open(GEN, $ARGV[5]) || die "Cannot open $ARGV[5]\n";
}

# load primer and target sequences
my %fwd; my %rev; my %targ;
while (my $line = <PR>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split(',', $line);
  next if (scalar @spl < 4);
  $fwd{$spl[0]} = $spl[1];
  $rev{$spl[0]} = $spl[2];
  $targ{$spl[0]} = $spl[3];
}
close PR;

# load expected amplicon location
my %loc;
while (my $line = <BED>) {
  next if (substr($line, 0, 1) eq '#');
  chomp $line;
  my @spl = split("\t", $line);
  if (scalar @spl < 4) {
    print STDERR "Warning! Improperly formatted line in $ARGV[2]: $line\n  ",
      "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
    next;
  }
  if (exists $loc{$spl[3]}) {
    my @div = split("\t", $loc{$spl[3]});
    $loc{$spl[3]} = "$spl[0]\t".($spl[2]+1)
      if ($spl[2]+1 < $div[1]);
  } else {
    $loc{$spl[3]} = "$spl[0]\t".($spl[2]+1);
  }
}
close BED;

# load (unique) reads from fastq
my $count = 0; my $uniq = 0;
my %seq;  # amplicons are keys, lengths are keys to %{$seq{$amp}},
          #   seqs are keys to %{$seq{$amp}{$len}},
          # values are lists of reads that share the seq (comma-separated)
my %qual; # for quality scores
while (my $line = <FQ>) {
  next if (substr($line, 0, 1) ne '@');
  chomp $line;
  my @spl = split(" ", $line);  # $spl[$#spl] should have cigar I/D

  # determine amplicon ID
  my $id = "";
  for (my $x = 1; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $id = $spl[$x-1];
      last;
    }
  }
  die "Error! $ARGV[0] is improperly formatted:\n  ",
    "no amplicon ID in $line\n" if (!$id);
  die "Error! $ARGV[0] is improperly formatted:\n  ",
    "no in/del info in $line\n" if ($spl[$#spl] !~ m/[ID]/);

  my $read = substr($spl[0], 1);
  $line = <FQ>;
  chomp $line;
  my $len = length $line;
  if (exists $seq{$id}{$len}{$line}) {
    $seq{$id}{$len}{$line} .= ",$read";
  } else {
    $seq{$id}{$len}{$line} = "$spl[$#spl] $read";
    $uniq++;
  }
  $count++;
  $line = <FQ>;
  # save quality scores
  $line = <FQ>;
  chomp $line;
  $qual{$read} = $line;
}
close FQ;
#print "Reads loaded: $count\n",
#  "Unique reads: $uniq\n";

# align reads to genomic segments
print LOG "Amplicon\tReads\tScore\tCIGAR(s)\tSequence\tReads"
  if (scalar @ARGV > 4);
print OUT "#Read\tAmplicon\tSubs\tIn/del\tChrom\tPos\tCIGAR",
  "\tSeq\tQual\tMD\tExternal\tScore\n";
foreach my $am (sort keys %seq) {

  # loop through length variants
  foreach my $len (sort keys %{$seq{$am}}) {

    if (! exists $loc{$am}) {
      print STDERR "Warning! Skipping $am -- no location info\n";
      last;
    } elsif (! exists $targ{$am}) {
      print STDERR "Warning! Skipping $am -- no sequences loaded\n";
      last;
    }

    print LOG "\n$am len=$len\n" if (scalar @ARGV > 4);
    my %rep = ();  # for I/D and list of reads
    my %cig = ();  # saves counts of reads that match each cigar
    my $tot = 0;   # total number of reads

    # loop through unique reads
    while (my $max = findMax(%{$seq{$am}{$len}})) {
      my @que = split(" ", $max); # $que[0] is sequence, $que[1] is read count
      $tot += $que[1];
      my @spl = split(" ", $seq{$am}{$len}{$que[0]}); # $spl[0] is I/D, $spl[1] lists reads
      $spl[0] =~ m/(\d+)([ID])/;
      if (length $que[0] != ($2 eq 'I' ? $1 : -$1) + length $targ{$am}) {
        print STDERR "Error! $am, $spl[0] does not match sequences:\n",
          "seq  $que[0]\nref  $targ{$am}\n";
        delete $seq{$am}{$len}{$que[0]};
        next;
      }

      # get possible CIGARs
      my $aln = alignSeq($que[0], $targ{$am}, $1, ($2 eq 'I' ? 1 : 0));
      my @cut = split(" ", $aln);
      my $score = shift @cut;
      for (my $x = 0; $x < scalar @cut; $x++) {
        $cig{$cut[$x]} += $que[1];
      }

      # log info
      print LOG "\t$que[1]\t$score\t", join(",", @cut),
        "\t$que[0]\t$spl[1]\n" if (scalar @ARGV > 4);

      $rep{$que[0]} = "$spl[0] $spl[1]"; # save I/D and list of reads
      delete $seq{$am}{$len}{$que[0]};
    }

    # find consensus CIGAR
    my $max = 0;
    my $firstM = 1000;  # first M of cigar -- 0 indicates external in/del
    my $best = "";      # best cigar
    foreach my $key (keys %cig) {
      if ($cig{$key} > $max) {
        $max = $cig{$key};
        $best = $key;
        if ($key =~ m/\D0M$/) {
          $firstM = 0;
        } else {
          $key =~ m/^(\d+)M/;
          $firstM = $1;
        }
      } elsif ($cig{$key} == $max) {
        # if equal, push toward 5' end
        $key =~ m/^(\d+)M/;
        if ($1 < $firstM) {
          $firstM = $1;
          $best = $key;
        }
      }
    }
    # preference for external in/dels:
    if ($firstM) {
      foreach my $key (keys %cig) {
        if ((($key =~ m/^0M/) || ($key =~ m/\D0M$/)) &&
            $cig{$key} > 0.95*$max) {  # 5% less than max is OK
          $firstM = 0;
          $max = $cig{$key};
          $best = $key;
          last;
        }
      }
    }
    printf LOG "consensus\t$max\t%.1f%%\t$best\n", 100*$max/$tot
      if (scalar @ARGV > 4);
    print STDERR "Warning for $am, len=$len: $best not perfect\n"
      if ($max != $tot);

    # evaluate external in/del for artifacts (i.e. if primer is match)
    my $score = -1;
    if (!$firstM) {
      my $three = ($best =~ m/\D0M$/ ? 1 : 0);  # in/del at 3' end
      my $prim = ($three ? reverse $rev{$am} : $fwd{$am});  # primer sequence
      my $base = substr($prim, -1);  # 3' base of primer
      my $gen = "";  # genomic sequence
      if ($best =~ m/D/) {
        if ($three) {
          $gen = reverse(substr($targ{$am}.$rev{$am}, $len, length $prim));
        } else {
          $gen = substr($fwd{$am}.$targ{$am}, - $len - length $prim, length $prim);
        }
      } elsif (scalar @ARGV > 5) {
        # external insertion: must query whole genome for segment
        my @div = split("\t", $loc{$am});
        my $chr = "";
        local $/ = '>';
        my $waste = <GEN>;
        while (my $chunk = <GEN>) {
          chomp $chunk;
          my @spl = split("\n", $chunk);
          my @head = split(" ", shift @spl);
          my $ch = $head[0];
          if ($ch eq $div[0]) {
            $chr = join("", @spl);
            last;
          }
        }
        close GEN;
        open(GEN, $ARGV[5]) || die "Cannot open $ARGV[5]\n";
        if ($chr) {
          if ($three) {
            $gen = reverse(substr($chr, $div[1] - 1 + $len, length $prim));
          } else {
            $best =~ m/^0M(\d+)I/;
            $gen = substr($chr, $div[1] - 1 - $1 - length $prim, length $prim);
          }
          $gen =~ tr/a-z/A-Z/;
        }
      }
      if (!$gen) {
        print STDERR "Warning! Cannot evaluate external in/del of $am, len=$len",
          "  (no reference sequence loaded)\n";
      } else {
        # score match of primer to genomic segment
        $score = scoreAlign($gen, $prim);
        printf LOG "\texternal\t%.3f\n", $score if (scalar @ARGV > 4);
      }

      # adjust $best to include one base of primer
      if ($three) {
        substr($best, -2) = "1M"; # =~ s/0M\$/1M\$/;
      } else {
        substr($best, 0, 2) = "1M"; #$best =~ s/0M/1M/;
      }

      # produce output -- info for a SAM record
      while (my $max = findMax(%rep)) {
        my @que = split(" ", $max); # $que[0] is sequence, $que[1] is read count
        my @cut = split(" ", $rep{$que[0]}); # $cut[0] has in/del length
        my @spl = split(",", $cut[1]); # list of reads
        my $seq1 = ($three ? $que[0].$base : $base.$que[0]);
        my $seq2 = ($three ? $targ{$am}.$base : $base.$targ{$am});
        my ($md0, $md1) = getMD($seq1, $seq2, $best);
        my @lc = split("\t", $loc{$am});
        for (my $x = 0; $x < scalar @spl; $x++) {
          print OUT "$spl[$x]\t$am\t$md1\t$cut[0]\t$lc[0]\t",
            ($three ? $lc[1] : $lc[1]-1),
            "\t$best\t$seq1\t",
            ($three ? "$qual{$spl[$x]}I" : "I$qual{$spl[$x]}"),
            "\t$md0";
          printf OUT "\texternal\t%.3f", $score if ($score != -1);
          print OUT "\n";
        }
        delete $rep{$que[0]};
      }

      next;
    }

    # produce output -- info for a SAM record
    while (my $max = findMax(%rep)) {
      my @que = split(" ", $max); # $que[0] is sequence, $que[1] is read count
      my @cut = split(" ", $rep{$que[0]}); # $cut[0] has in/del length
      my @spl = split(",", $cut[1]); # list of reads
      my ($md0, $md1) = getMD($que[0], $targ{$am}, $best);
      for (my $x = 0; $x < scalar @spl; $x++) {
        print OUT "$spl[$x]\t$am\t$md1\t$cut[0]\t$loc{$am}\t$best\t$que[0]",
          "\t$qual{$spl[$x]}\t$md0";
        printf OUT "\texternal\t%.3f", $score if ($score != -1);
        print OUT "\n";
      }
      delete $rep{$que[0]};
    }
  }
}
close OUT;
close LOG if (scalar @ARGV > 4);
close GEN if (scalar @ARGV > 5);

# produce MD flag
sub getMD {
  my $que = $_[0];
  my $ref = $_[1];
  my $cig = $_[2];

  my $md = "";
  my $qpos = 0; my $rpos = 0;
  my $match = 0;
  my $sub = 0;
  while ($cig =~ m/(\d+)([IMD])/g) {
    my $len = $1;
    if ($2 eq 'M') {
      for (my $x = 0; $x < $len; $x++) {
        if (substr($que, $qpos + $x, 1) ne
            substr($ref, $rpos + $x, 1)) {
          $md .= $match . substr($ref, $rpos + $x, 1);
          $match = 0;
          $sub++;
        } else {
          $match++;
        }
      }
      $rpos += $len;
      $qpos += $len;
    } elsif ($2 eq 'D') {
      if ($match) {
        $md .= $match;
      } elsif (!$md) {
        $md = "0";
      }
      $match = 0;
      $md .= '^' . substr($ref, $rpos, $len);
      $rpos += $len;
    } elsif ($2 eq 'I') {
      $qpos += $len;
    }
  }
  $md .= $match if ($match);
  return $md, $sub;
}

# align sequences -- return cigar(s)
sub alignSeq {
  my $que = $_[0];
  my $ref = $_[1];
  my $gap = $_[2];
  my $ins = $_[3];

  # rearrange so $que is shorter sequence
  if ($ins) {
    my $temp = $que;
    $que = $ref;
    $ref = $temp;
  }

  my $maxscore = -1000;  # best alignment score (actual max. is 0)
  my $cig = "";   # best alignment cigar

  # check each position -- $x is 5' end of gap
  for (my $x = 0; $x < 1 + length $que; $x++) {
    my $score = 0;
    my $rpos = 0;
    my $flag = 0;
    for (my $y = 0; $y < length $que; $y++) {
      $rpos += $gap if ($y == $x);
      if (substr($que, $y, 1) ne substr($ref, $rpos + $y, 1)) {
        if (--$score < $maxscore) {
          $flag = 1;
          last;
        }
      }
    }

    # save cigar
    if ($flag) {
      next;
    } elsif ($score == $maxscore) {
      $cig .= " ${x}M$gap" . ($ins ? "I" : "D") . ((length $que) - $x) . "M";
    } else {
      $maxscore = $score;
      $cig = "${x}M$gap" . ($ins ? "I" : "D") . ((length $que) - $x) . "M";
    }
  }

  return "$maxscore $cig";
}

# find sequence in hash with the most reads
sub findMax {
  my %hash = @_;
  my $max = 0;
  my $seq = "";
  foreach my $key (sort keys %hash) {
    my @spl = split(",", $hash{$key});
    if (scalar @spl > $max) {
      $max = scalar @spl;
      $seq = $key;
    }
  }
  return $max ? "$seq $max" : "";
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
  my $len = length $que;

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
  return $score / $total;
}
