#!/usr/bin/perl

# John M. Gaspar (jmgaspar@gwu.edu)
# Jan. 2015

# Add "HP" tag to VCF that lists homopolymer length
#   adjacent to each variant.

# Analyzes only 10bp on either side of variant --
#   cf. lines 187-89, or <logfile> output.

use strict;
use warnings;

sub usage {
  print q(Usage: perl addHPtoVCF.pl  <infile>  <genome>  <outfile>  [<logfile>]
  Required:
    <infile>   Input VCF file -- should list one variant
                 per line, with INFO field "CIGAR"
    <genome>   Fasta file of reference genome
    <outfile>  Output VCF file
  Optional:
    <logfile>  Verbose log file (if selected, input VCF
                 should have FORMAT field "AF")
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

# open files
open(VCF, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(GEN, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
$/ = '>';
my $waste = <GEN>;
$/ = "\n";
open(OUT, ">$ARGV[2]");
if (scalar @ARGV > 3) {
  open(LOG, ">$ARGV[3]");
  print LOG "#CHROM\tPOS\tREF\tALT\tAF\tCIGAR\tGenomeSegment\tHP\n";
}

# analyze VCF file
my $chr = "";
my $seq = "";
my $pr = 1;  # flag for header printing
while (my $line = <VCF>) {
  if (substr($line, 0, 1) eq '#') {
    if ($pr && substr($line, 0, 8) eq '##FORMAT') {
      print OUT "##INFO=<ID=HP,Number=1,Type=Integer,",
        "Description=\"Length of adjacent homopolymer that matches variant\">\n";
      $pr = 0;
    }
    print OUT $line;
    next;
  }
  chomp $line;
  my @spl = split("\t", $line);
  die "Error! $ARGV[0] is improperly formatted\n" if (scalar @spl < 10);

  # skip if multiple alleles on one line
  my @b1 = split(',', $spl[3]);
  my @b2 = split(',', $spl[4]);
  if (scalar @b1 > 1 || scalar @b2 > 1) {
    print "Warning! Multiple variant alleles in VCF record:\n$line\n";
    print OUT "$line\n";
    next;
  }

  # determine allele frequency
  my $af = 0;
  if (scalar @ARGV > 3) {
    my @div = split(':', $spl[8]);
    my $idx = -1;
    for (my $x = 0; $x < scalar @div; $x++) {
      if ($div[$x] eq "AF") {
        $idx = $x;
        last;
      }
    }
    my @cut = split(':', $spl[9]);
    die "Error! Cannot find \"AF\" in VCF record:\n$line\n"
      if ($idx == -1 || $idx >= scalar @cut);
    $af = $cut[$idx];
  }

  # determine position of variant(s) using CIGAR
  if ($spl[7] !~ m/CIGAR\=([0-9DIMX]+)/) {
    die "Error! Cannot find \"CIGAR\" in VCF record:\n$line\n";
  }
  my $cig = $1;
  my @loc = ();
  my $pos = 0; my $sub = 0;
  my $alt = '*';
  my $prev = "";
  my $base = "";  # different nt
  my $flag = 0;   # 1 if complex variant
  my $del = "";   # deleted bases if variant is a deletion
  my $hit = 0;    # length of homopolymer
  while ($cig =~ m/(\d+)([IDMX])/g) {
    if ($2 ne 'M') {
      if ($prev) {
        $alt .= " $prev *";
        $prev = "";
      }
      for (my $x = 0; $x < $1; $x++) {
        push @loc, $spl[1] + $pos + $x;
        if ($2 eq 'D') {
          if ($base) {
            $flag = 1 if ($base ne substr($spl[3], $pos+$x, 1));
          } else {
            $del .= substr($spl[3], $pos+$x, 1);
          }
          $alt .= substr($spl[3], $pos+$x, 1);
        } else {
          if ($base) {
            $flag = 1 if ($base ne substr($spl[4], $sub+$x, 1));
          } else {
            for (my $y = 0; $y < length $del; $y++) {
              if (substr($spl[4], $sub+$x, 1) ne substr($del, $y, 1)) {
                $flag = 1;  # complex variant, base does not match del
                last;
              }
            }
            $base = substr($spl[4], $sub+$x, 1);
          }
          $alt .= substr($spl[4], $sub+$x, 1) if ($2 ne 'I');
        }
      }
      $alt .= '*';
    } else {
      if ($pos) {
        $prev = substr($spl[3], $pos, $1);
      }
      if ($base && length $' != 0) {
        # matching bases in complex variant must match variant bases
        for (my $x = 0; $x < $1; $x++) {
          $flag = 1 if ($base ne substr($spl[4], $sub+$x, 1));
          $hit++;  # add to homopolymer length
        }
      }
    }
    last if ($flag);
    $pos += $1 if ($2 ne 'I');
    $sub += $1 if ($2 ne 'D');
  }
  if (scalar @loc == 0) {
    die "Error! No variant position in $line\n";
  }

  if ($flag) {
    # complex variant
    print LOG "$spl[0]\t$spl[1]\t$spl[3]\t$spl[4]",
      "\t$af\t$cig\tnot checked\n" if (scalar @ARGV > 3);
    $hit = 0;
  } else {

    # if not already loaded, get chromosome from genome
    if ($spl[0] ne $chr) {
      $chr = $spl[0];
      $seq = "";
      for (my $x = 0; $x < 2; $x++) {
        local $/ = '>';
        while (my $chunk = <GEN>) {
          chomp $chunk;
          my @div = split("\n", $chunk);
          my $ch = shift @div;
          if ($spl[0] eq $ch) {
            $seq = join("", @div);
            last;
          }
        }
        if (! $seq) {
          # if not found, try again
          close GEN;
          open(GEN, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
          my $waste = <GEN>;
        } else {
          last;
        }
      }
    }
    if (! $seq) {
      die "Error! Cannot find chromosome $chr in $ARGV[1]\n";
    }

    # judge chrom segments (10bp on either side of variant)
    my $seg1 = substr($seq, $loc[0]-11, 10);
    my $seg2 = substr($seq, $loc[$#loc], 10);
    $seg1 =~ tr/a-z/A-Z/;
    $seg2 =~ tr/a-z/A-Z/;

    # add homopolymer lengths
    if ($del) {
      # for deletion, do not need to match $base
      $seg1 =~ m/((.)\2*)$/;
      my $hit1 = length $1;
      $seg2 =~ m/^($2*)/;
      $hit1 += length $1;

      # check other side
      $seg2 =~ m/^((.)\2*)/;
      my $hit2 = length $1;
      $seg1 =~ m/($2*)$/;
      $hit2 += length $1;

      $hit += ($hit2 > $hit1 ? $hit2 : $hit1);  # save maximum

    } else {
      $seg1 =~ m/($base*)$/;
      $hit += length $1;
      $seg2 =~ m/^($base*)/;
      $hit += length $1;
    }
    print LOG "$spl[0]\t$spl[1]\t$spl[3]\t$spl[4]\t$af",
      "\t$cig\t$seg1 $alt $seg2\t$hit\n" if (scalar @ARGV > 3);
  }

  # print new VCF record
  $spl[7] .= ";HP=$hit";
  print OUT join("\t", @spl), "\n";

}
close GEN;
close VCF;
close OUT;
close LOG if (scalar @ARGV > 3);
