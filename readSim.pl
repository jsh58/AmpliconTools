#!/usr/bin/perl

# John M. Gaspar
# Mar. 2015

# Simulate amplicon-based paired-end reads.

use strict;
use warnings;

sub usage {
  print q(Usage: perl readSim.pl  <infile1>  <infile2>  <genome>  <output1>  <output2>  \
    <number>  <length>  <score>  <lengths>  <infile3>  <infile4>  <percent>
  Required:
    <infile1>  Output from ipcress (or filterIpc.pl)
    <infile2>  File listing primer sequences (input to ipcress)
    <genome>   Fasta file of reference genome
    <output1>  Output file -- PE reads #1
    <output2>  Output file -- PE reads #2
  Optional:
    <number>   Number of PE reads for each perfect primer match score --
                 others will be scaled (using <score> parameter) (def. 1000)
    <length>   Length of PE reads to create (def. 100)
    <score>    Minimum primer matching score (scale 0-1; def. 0.75)
    <lengths>  Allowed amplicon lengths (e.g. 70,160 -> 70-160bp amplicons; def. 0,100000)
    <infile3>  BED file giving expected amplicon locations --
                 alternative overlapping amplicons will be excluded
    <infile4>  File listing variants to make in the reads --
                 should list chromosome, position, reference, and alternate
                 (tab-delimited), with alleles in VCF-primitive style (no MNPs
                 or complex variants).  For example:
                   chr3    109432461   T    C
                   chr17   41343130    GC   G
                   chrX    111077593   A    AACCTCCG
    <percent>  Percent of reads to make variants from <infile4> (def. 10)
);
  exit;
}

usage() if (scalar @ARGV < 5 || $ARGV[0] eq "-h");

# open files
open(IP, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(PR, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(FA, $ARGV[2]) || die "Cannot open $ARGV[2]\n";
open(FQ1, ">$ARGV[3]");
open(FQ2, ">$ARGV[4]");

# load optional parameters:

# number of reads for a perfect primer match --
#   lesser matches will get a proportion thereof
my $num = 1000;
$num = $ARGV[5] if (scalar @ARGV > 5);
die "Number of reads must be greater than 0\n" if ($num <= 0);

# length of PE reads
my $len = 100;
$len = $ARGV[6] if (scalar @ARGV > 6);
die "Length of reads must be greater than 0\n" if ($len <= 0);

# primer match score
my $pct = 0.75;
$pct = $ARGV[7] if (scalar @ARGV > 7);
die "Minimum primer matching score must be in [0,1)\n"
  if ($pct < 0 || $pct >= 1);

# amplicon lengths -- min and max
my $min = 0;
my $max = 100000;
if (scalar @ARGV > 8) {
  my @spl = split(',', $ARGV[8]);
  if (scalar @spl > 1) {
    $min = $spl[0];
    $max = $spl[1];
  } else {
    print "Warning! Cannot load amplicon lengths $ARGV[8]\n",
      "Using defaults of $min,$max\n";
  }
  # make sure $min is less than $max
  if ($min > $max) {
    my $temp = $min;
    $min = $max;
    $max = $temp;
  }
}

# load expected amplicon locations
my %exp;  # exact expected location
if (scalar @ARGV > 9 && open(BED, $ARGV[9])) {
  while (my $line = <BED>) {
    chomp $line;
    my @spl = split("\t", $line);
    if (scalar @spl < 4) {
      print "Warning! Improperly formatted line in $ARGV[9]: $line\n  ",
        "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
      next;
    }
    if (exists $exp{$spl[3]}) {
      my @div = split("\t", $exp{$spl[3]});
      if ($div[1] < $spl[1]) {
        $exp{$spl[3]} .= "\t$spl[1]";
      } else {
        $exp{$spl[3]} = "$div[0]\t$spl[1]\t$div[1]";
      }
    } else {
      $exp{$spl[3]} = "$spl[0]\t$spl[1]";
    }
  }
  close BED;
}

# load variants
my %var;
my $per = 0;
if (scalar @ARGV > 10) {
  my $count = 0;
  open(VAR, $ARGV[10]) || die "Cannot open $ARGV[10]\n";
  while (my $line = <VAR>) {
    chomp $line;
    my @spl = split("\t", $line);
    next if (scalar @spl < 4);
    $var{$spl[0]}{$spl[1]} = "$spl[2]\t$spl[3]";
    $count++;
  }
  close VAR;
  print "Variants loaded: $count\n";
  $per = (scalar @ARGV > 11 ? $ARGV[11] : 10);
}

# adapter sequences: reads from shorter amplicons may contain part of these sequences 
my $il1 = "AGACCAAGTCTCTGCTACCGTACTCTGGACGAATCTCGTATGCCGTCTTCTGCTTGAAAA";  # fwd reads
my $il2 = "TGTAGAACCATGTCGTCAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAA";  # rev reads
while (length $il1 < $len) {
  # expand to min. length ($len) with "A"
  $il1 .= "A";
  $il2 .= "A";
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
my $total = 0;
my $count = 0;
my $vars = 0;
my $xlen = 0;
my $xover = 0;
while (my $line = <IP>) {
  chomp $line;
  my @spl = split(" ", $line);

  # skip if from one primer only
  next if ($spl[0] ne "ipcress:");
  next if ($spl[10] eq "single_A" || $spl[10] eq "single_B");

  # skip if no genomic segment loaded or no primers loaded
  my @div = split(':', $spl[1]);
  if (! exists $chr{$div[0]}) {
    print "Warning! No sequence loaded for reference $div[0]\n";
    next;
  }
  if ((! exists $fwd{$spl[2]}) || (! exists $rev{$spl[2]})) {
    print "Warning! No primers loaded for amplicon $spl[2]\n";
    next;
  }

  # get primers
  my $pr5 = ($spl[4] eq 'B' ? $rev{$spl[2]} : $fwd{$spl[2]});
  my $pr3 = ($spl[4] eq 'B' ? $fwd{$spl[2]} : $rev{$spl[2]});

  # score primer-genome alignment
  my $scoreF; my $scoreR;
  if (scalar @spl < 13) {
    # not scored before -- do it now
    my $fwdP = substr($chr{$div[0]}, $spl[5], length $pr5);
    my $revP = revComp(substr($chr{$div[0]}, $spl[8], length $pr3));
    $fwdP =~ tr/a-z/A-Z/;
    $revP =~ tr/a-z/A-Z/;
    $scoreF = scoreAlign($fwdP, $pr5);
    $scoreR = scoreAlign($revP, $pr3);
  } else {
    $scoreF = $spl[11];
    $scoreR = $spl[12];
  }

  # create simulated reads
  my $tot = 0;
  if ($scoreF >= $pct && $scoreR >= $pct) {

    # calculate proportion for these primers --
    #   $pct => 1%; 1.00 => 100%
    my $propF = (0.99 * $scoreF + 0.01 - $pct) / (1 - $pct);
    my $propR = (0.99 * $scoreR + 0.01 - $pct) / (1 - $pct);
    my $prop = $propF * $propR;
    $tot = int($prop * $num + 0.5);  # number of reads to create
  }
  next if ($tot < 1);

  # skip if length outside designated lengths
  if ($spl[3] < $min || $spl[3] > $max) {
    $xlen += $tot;
    next;
  }

  # skip if overlapping amplicon
  if (exists $exp{$spl[2]}) {
    my @cut = split("\t", $exp{$spl[2]});
    if (scalar @cut > 2 && $div[0] eq $cut[0] &&
        (($spl[5] > $cut[1] - 100 && $spl[5] < $cut[2] + 100) ||
        ($spl[8] > $cut[1] - 100) && ($spl[8] < $cut[2] + 100))) {
      # exact match is OK, of course
      if ($spl[5] != $cut[1] || $spl[8] != $cut[2]) {
        $xover += $tot;
        next;
      }
    }
  }

  # create base sequence
  my $head = "amp=$spl[2] ref=$div[0] pos=$spl[5] len=$spl[3]" .
    " str=" . ($spl[10] eq "forward" ? '+' : '-');
  my $st = $spl[5] + length $pr5;  # target starting position
  my $flseq = $pr5 . substr($chr{$div[0]}, $st, $spl[8] - $st)
      . revComp($pr3);  # full-length amplicon sequence
  $flseq =~ tr/a-z/A-Z/;
  my $qual = "\n+\n" . 'H' x $len . "\n";

  # load variants that apply to this amplicon
  my %vr;
  foreach my $pos (keys %{$var{$div[0]}}) {
    my @cut = split("\t", $var{$div[0]}{$pos});
    if (length $cut[0] == length $cut[1]) {
      # SNP only -- no MNP support (right now)
      $vr{$pos} = $var{$div[0]}{$pos} if ($pos > $st && $pos < $spl[8]+1);
    } elsif (length $cut[0] > length $cut[1]) {
      # deletion
      $vr{$pos} = $var{$div[0]}{$pos} if ($pos+1 > $st && $pos-1+length $cut[0] < $spl[8]+1);
    } else {
      # insertion
      $vr{$pos} = $var{$div[0]}{$pos} if ($pos+1 > $st && $pos < $spl[8]+1);
    }
  }

  # create variant sequence, from 3' end (for in/del-position issue)
  my $varseq = $flseq;
  my $mut = " var=";
  foreach my $pos (sort {$b <=> $a} keys %vr) {
    my @cut = split("\t", $vr{$pos});
    my $base = substr($varseq, $pos-$spl[5]-1, length $cut[0]);
    if ($base ne $cut[0]) {
      print "Warning! Sequence $base does not match listed ",
        "variant $cut[0] > $cut[1] at $div[0], $pos\n";
    } else {
      substr($varseq, $pos-$spl[5]-1, length $cut[0], $cut[1]);
      $mut .= "$pos,$cut[0]>$cut[1];";
    }
  }

  # reverse-complement sequences if necessary
  if ($spl[4] eq 'B') {
    $flseq = revComp($flseq);
    $varseq = revComp($varseq) if ($per);
  }

  # create PE reads
  my $seq1; my $seq2;
  if ($spl[3] < $len) {
    # short amplicon -- need to fill in 3' end with "adapter" seq.
    $seq1 = $flseq . substr($il1, 0, $len - $spl[3]);
    $seq2 = revComp($flseq) . substr($il2, 0, $len - $spl[3]);
  } else {
    $seq1 = substr($flseq, 0, $len);
    $seq2 = revComp(substr($flseq, -$len));
  }

  # create variant reads as well
  my $varseq1; my $varseq2;
  if ($varseq ne $flseq) {
    if (length $varseq < $len) {
      $varseq1 = $varseq . substr($il1, 0, $len - length $varseq);
      $varseq2 = revComp($varseq) . substr($il2, 0, $len - length $varseq);
    } else {
      $varseq1 = substr($varseq, 0, $len);
      $varseq2 = revComp(substr($varseq, -$len));
    }
  }

  # print reads
  my $vs = $per * $tot / 100;
  for (my $x = 0; $x < $tot; $x++) {
    my $beg = "\@read$count ";
    if ($x < $vs && $varseq ne $flseq) {
      print FQ1 $beg, $head, $mut, "\n", $varseq1, $qual;
      print FQ2 $beg, $head, $mut, "\n", $varseq2, $qual;
      $vars++;
    } else {
      print FQ1 $beg, $head, "\n", $seq1, $qual;
      print FQ2 $beg, $head, "\n", $seq2, $qual;
    }
    $count++;
  }
  $total++;

}
close IP;
close FQ1;
close FQ2;

print "Valid amplicons: $total",
  "\nReads simulated: $count",
  "\n";
print "  with variants: $vars\n" if (scalar @ARGV > 10);
print "Reads excluded due to length: $xlen\n" if (scalar @ARGV > 8);
print "Reads excluded due to overlapping: $xover\n" if (scalar @ARGV > 9);

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
