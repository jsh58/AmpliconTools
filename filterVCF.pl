#!/usr/bin/perl

# John M. Gaspar
# July 2015

# Filter a VCF file, according to the following options
#   (in the order in which they are evaluated):
#   - limiting variants to specified target regions
#   - minimum read depth
#   - minimum number of variant allele observations
#   - minimum allele frequency
#   - minimum allele frequency for in/dels
#   - minimum allele frequency for C:G>T:A substitutions
#   - minimum allele frequency for variants in
#       homopolymer runs

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterVCF.pl  <infile>  <outfile>  [options]
  Required:
    <infile>     Input VCF file -- should list one variant per line,
                   with INFO fields CIGAR, TYPE, HP, and FORMAT fields
                   AF, AO, DP, RO (produced by makeVCF.pl, addHPtoVCF.pl)
    <outfile>    Output VCF file
  Filtering options (in order of precedence):
    -b  <file>   BED file listing locations of primers -- variants
                   outside of target regions will be eliminated
    -d  <int>    Minimum read depth
    -o  <int>    Minimum variant allele observations
    -a  <float>  Minimum allele frequency [0-1]
    -i  <float>  Minimum allele frequency for in/dels [0-1]
    -m  <float>  Minimum allele frequency for C:G>T:A substitutions
                   (possible deamination artifacts) [0-1]
    -p  <str>    For filtering variants in homopolymer runs. Argument
                   should contain HP length, followed by a minimum allele
                   frequency for variants (or two values, for subs and
                   in/dels separately), comma-separated.  For example:
                     -p 5,0.2      # removes all variants with
                                   #   HP >= 5 and AB < 0.2
                     -p 7,0.5,1.1  # removes all substitutions with
                                   #   HP >= 7 and AB < 0.5; removes
                                   #   all in/dels with HP >= 7
                   Multiple HP ranges can be chosen (colon-separated):
                     -p 5,0.1,0.2:8,0.4,0.5
                             # for HP = [5, 6, 7]: removes subs with
                             #   AB < 0.1 and in/dels with AB < 0.2;
                             # for HP >= 8: removes subs with AB < 0.4
                             #   and in/dels with AB < 0.5
  Other option:
    -ve          Option to print summary counts to STDOUT
);
  exit;
}

usage() if (scalar @ARGV < 2 || $ARGV[0] eq "-h");

# open files
open(IN, $ARGV[0]) || die "Cannot open $ARGV[0]\n";

# get other command-line parameters
my @opts;
my %reg;  # target regions
my %minHP; my $maxHP = 100;  # maximum HP value
my $minAF = 0; my $minAO = 0; my $minDP = 0;
my $minID = 0; my $minCT = 0;
my $verb = 0;
for (my $x = 2; $x < scalar @ARGV; $x++) {
  if ($ARGV[$x] eq "-h") {
    usage();
  } elsif ($ARGV[$x] eq "-ve") {
    $verb = 1;
  } elsif ($x < scalar @ARGV - 1) {
    if ($ARGV[$x] eq "-a") {
      $minAF = $ARGV[++$x];
      push @opts, "AF >= $minAF";
    } elsif ($ARGV[$x] eq "-o") {
      $minAO = $ARGV[++$x];
      push @opts, "AO >= $minAO";
    } elsif ($ARGV[$x] eq "-d") {
      $minDP = $ARGV[++$x];
      push @opts, "DP >= $minDP";
    } elsif ($ARGV[$x] eq "-i") {
      $minID = $ARGV[++$x];
      push @opts, "AF(in/del) >= $minID";
    } elsif ($ARGV[$x] eq "-m") {
      $minCT = $ARGV[++$x];
      push @opts, "AF(C:G>T:A) >= $minCT";
    } elsif ($ARGV[$x] eq "-p") {

      # parse homopolymer run parameter
      my @spl = split(':', $ARGV[++$x]);
      my $min = $maxHP;
      for (my $y = 0; $y < scalar @spl; $y++) {
        my @div = split(',', $spl[$y]);
        die "Error! HP value is not an integer: $spl[$y]\n"
          if ($div[0] != int $div[0]);
        my $hp = int(shift @div);
        die "Error! HP value improperly formatted: $spl[$y]\n"
          if (scalar @div != 1 && scalar @div != 2);
        $minHP{$hp} = join(',', @div);
        $min = $hp if ($hp < $min);
      }
      my $prev = "";
      for (my $y = $min; $y < $maxHP; $y++) {
        if (exists $minHP{$y}) {
          $prev = $minHP{$y};
        } else {
          $minHP{$y} = $prev;
        }
      }
      push @opts, "HP: $ARGV[$x]";

    } elsif ($ARGV[$x] eq "-b") {
      open(BED, $ARGV[++$x]) || die "Cannot open $ARGV[$x]\n";
      push @opts, "bedfile = $ARGV[$x]";

      # load target regions (convert to 1-based)
      my %pos;
      while (my $line = <BED>) {
        chomp $line;
        my @spl = split("\t", $line);
        if (scalar @spl < 4) {
          print "Warning! Improperly formatted line in $ARGV[$x]: $line\n  ",
            "Need chromName, chromStart, chromEnd, and ampliconName (tab-delimited)\n";
          next;
        }
        if (exists $pos{$spl[3]}) {
          my @div = split("\t", $pos{$spl[3]});
          if ($div[0] ne $spl[0]) {
            print "Warning! Skipping amplicon $spl[3] -- ",
              "located at chromosomes $spl[0] and $div[0]!?\n";
          } else {
            my $st; my $end;
            if ($spl[1] < $div[1]) {
              $st = $spl[2] + 1;
              $end = $div[1] + 1;
            } else {
              $st = $div[2] + 1;
              $end = $spl[1] + 1;
            }
            for (my $y = $st; $y < $end; $y++) {
              $reg{"$spl[0]\t$y"} = 1;
            }
          }
          delete $pos{$spl[3]};
        } else {
          $pos{$spl[3]} = "$spl[0]\t$spl[1]\t$spl[2]";
        }
      }
      close BED;
      die "Error! No target regions loaded from $ARGV[$x]\n"
        if (scalar keys %reg == 0);

    } else {
      die "Error! Unknown CL option: $ARGV[$x]\n";
    }
  } else {
    die "Error! Unknown CL option / parameter needed: $ARGV[$x]\n";
  }
}

# parse vcf file, produce output
my $var = 0; my $good = 0;  # counting variables
my $xAF = 0; my $xAO = 0; my $xDP = 0;
my $xID = 0; my $xCT = 0; my $xBed = 0;
my $xHP = 0;
my $pr = 1;  # flag for header printing
open(OUT, ">$ARGV[1]");
while (my $line = <IN>) {
  if (substr($line, 0, 1) eq '#') {
    if ($pr && substr($line, 0, 6) eq '##INFO'
        && scalar @opts > 0) {
      # print filter options to header
      @opts = sort @opts;
      print OUT "##filter=\"";
      for (my $x = 0; $x < scalar @opts - 1; $x++) {
        print OUT "$opts[$x]; "
      }
      print OUT "$opts[$#opts]\"\n";
      $pr = 0;
    }
    print OUT $line;
    next;
  }
  $var++;

  # process variant
  chomp $line;
  my @spl = split("\t", $line);
  die "Error! $ARGV[0] is improperly formatted\n" if (scalar @spl < 10);

  # skip if multiple alleles on one line
  my @b1 = split(',', $spl[3]);
  my @b2 = split(',', $spl[4]);
  if (scalar @b1 > 1 || scalar @b2 > 1) {
    print "Warning! Multiple variant alleles in VCF record:\n$line\n";
    print OUT "$line\n";
    $good++;
    next;
  }

  # check if variants are outside of target regions
  if (scalar keys %reg > 0) {

    # determine position of variant(s) using CIGAR
    if ($spl[7] !~ m/CIGAR\=([0-9DIMX]+)/) {
      die "Error! Cannot find \"CIGAR\" in VCF record:\n$line\n";
    }
    my $cig = $1;
    my $targ = 0;  # flag for location within a target
    my $pos = 0;
    while ($cig =~ m/(\d+)([IDMX])/g) {
      if ($2 ne 'M') {
        for (my $y = 0; $y < $1; $y++) {
          my $loc = $spl[1] + $pos + $y;
          if (exists $reg{"$spl[0]\t$loc"}) {
            $targ = 1;
            last;
          }
          last if ($2 eq 'I');
        }
      } 
      last if ($targ);
      $pos += $1 if ($2 ne 'I');
    }

    if (! $targ) {
      $xBed++;
      next;
    }
  }

  # load other info for variant: HP, TYPE
  my $type = 0;  # 0 for sub, 1 for ins, 2 for del
  if ($spl[7] =~ m/TYPE\=ins/) {
    $type = 1;
  } elsif ($spl[7] =~ m/TYPE\=del/) {
    $type = 2;
  } elsif ($spl[7] !~ m/TYPE\=sub/) {
    die "Error! Unknown TYPE in VCF record:\n$line\n";
  }
  my @hp = (0);
  if (scalar keys %minHP > 0) {
    if ($spl[7] !~ m/HP\=(\d+)/) {
      die "Error! Cannot find \"HP\" in VCF record:\n$line\n";
    }
    if (exists $minHP{$1}) {
      @hp = split(',', $minHP{$1});
    }
  }

  # load FORMAT fields: AF, AO, DP
  my @div = split(':', $spl[8]);
  my $idxAF = -1; my $idxAO = -1; my $idxDP = -1;
  for (my $x = 0; $x < scalar @div; $x++) {
    if ($div[$x] eq "AF") {
      $idxAF = $x;
    } elsif ($div[$x] eq "AO") {
      $idxAO = $x;
    } elsif ($div[$x] eq "DP") {
      $idxDP = $x;
    }
  }
  my @cut = split(':', $spl[9]);
  die "Error! Cannot find \"AF\" in VCF record:\n$line\n"
    if ($idxAF == -1 || $idxAF >= scalar @cut);
  die "Error! Cannot find \"AO\" in VCF record:\n$line\n"
    if ($idxAO == -1 || $idxAO >= scalar @cut);
  die "Error! Cannot find \"DP\" in VCF record:\n$line\n"
    if ($idxDP == -1 || $idxDP >= scalar @cut);
  my $af = $cut[$idxAF];
  my $ao = $cut[$idxAO];
  my $dp = $cut[$idxDP];

  # check other filtering criteria
  if ($dp < $minDP) {
    $xDP++;
  } elsif ($ao < $minAO) {
    $xAO++;
  } elsif ($af < $minAF) {
    $xAF++;
  } elsif ($type && $af < $minID) {
    $xID++;
  } elsif ( (($spl[3] eq 'C' && $spl[4] eq 'T') ||
             ($spl[3] eq 'G' && $spl[4] eq 'A'))
            && $af < $minCT) {
    $xCT++;
  } elsif ( ((scalar @hp == 1 || (!$type)) && $af < $hp[0])
      || (scalar @hp > 1 && $type && $af < $hp[1]) ) {
    $xHP++;
  } else {
    print OUT "$line\n";
    $good++;
  }
}
close IN;
close OUT;

if ($verb) {
  print "Variants in $ARGV[0]: $var",
    "\nVariants printed to $ARGV[1]: $good",
    "\nVariants removed: ", $var - $good;
  print "\n  Outside of target regions: $xBed" if (scalar keys %reg > 0);
  print "\n  Read depth less than $minDP: $xDP" if ($minDP);
  print "\n  Variant allele observations less than $minAO: $xAO" if ($minAO);
  print "\n  Allele frequency less than $minAF: $xAF" if ($minAF);
  print "\n  In/del allele frequency less than $minID: $xID" if ($minID);
  print "\n  C:G>T:A allele frequency less than $minCT: $xCT" if ($minCT);
  print "\n  Variant in/near homopolymer run: $xHP" if (scalar keys %minHP > 0);
  print "\n";
}
