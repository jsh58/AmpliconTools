#!/usr/bin/perl

# John M. Gaspar
# June 2015

# Filter singleton reads after separate primer removal.
# Options to remove duplicates and chimeras.

use strict;
use warnings;

sub usage {
  print q(Usage: perl filterSingle.pl  <in1>  <in2>  <out>  [options]
  Required:
    <in1>   Input FASTQ file #1, with primers removed and amplicon identification
              in header (produced by removePrimer)
    <in2>   Input FASTQ file #2, with primers removed and amplicon identification
              in header (produced by removePrimer)
    <out>   Output FASTQ file with filtered singletons
  Options for duplicates (in order of precedence):
     -c     Option to remove both singletons that had different primers
              removed (like they were derived from a PCR chimera)
     -b     Option to favor singletons that have both primers removed
              over those that have not
     -q     Option to keep only higher quality singleton
  Other options:
     -ve    Option to print summary counts to STDOUT
);
  exit;
}

usage() if (scalar @ARGV < 3 || $ARGV[0] eq "-h");

# open files
open(FQ1, $ARGV[0]) || die "Cannot open $ARGV[0]\n";
open(FQ2, $ARGV[1]) || die "Cannot open $ARGV[1]\n";
open(OUT, ">$ARGV[2]");

# save options
my $chim = 0;
my $bth = 0;
my $qul = 0;
my $verb = 0;
for (my $x = 3; $x < scalar @ARGV; $x++) {
  if ($ARGV[$x] eq "-c") {
    $chim = 1;
  } elsif ($ARGV[$x] eq "-b") {
    $bth = 1;
  } elsif ($ARGV[$x] eq "-q") {
    $qul = 1;
  } elsif ($ARGV[$x] eq "-ve") {
    $verb = 1;
  }
}

# load info from first file
my %pr;   # which amplicon the read came from
my %qual; # avg. quality score
my %seq;  # header, sequence, and qual scores
my %both; # boolean if read had both primers removed
my $count = 0;
while (my $q = <FQ1>) {
  next if (substr($q, 0, 1) ne "@");
  my $w = <FQ1>;
  my $e = <FQ1>;
  my $r = <FQ1>;
  my @spl = split(" ", $q);
  die "Error! Read $spl[0] is duplicated in $ARGV[0]\n"
    if (exists $seq{$spl[0]});

  # determine amplicon ID
  my $id = "";
  my $bot = 0;
  for (my $x = 1; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $id = $spl[$x-1];
      $bot = 1 if ($x + 1 <= $#spl && $spl[$x+1] eq "both");
      last;
    }
  }
  die "Error! $ARGV[0] is improperly formatted:\n",
    "  no amplicon ID in $q\n" if (!$id);

  $pr{$spl[0]} = $id;
  $seq{$spl[0]} = $q . $w . $e . $r;
  $both{$spl[0]} = $bot;

  # save avg qual score
  if ($qul) {
    chomp $r;
    my $tot = 0;
    for (my $x = 0; $x < length $r; $x++) {
      $tot += ord(substr($r, $x, 1)) - 33; 
    }
    $qual{$spl[0]} = $tot / length $r;
  }

  $count++;
}
close FQ1;
print "Reads in $ARGV[0]: $count\n" if ($verb);

# parse second file
my $print = 0; my $crem = 0;
my $brem = 0; my $qrem = 0;  # counting variables
$count = 0;
my %dup;  # to check for duplicates
while (my $q = <FQ2>) {
  next if (substr($q, 0, 1) ne "@");
  my $w = <FQ2>;
  my $e = <FQ2>;
  my $r = <FQ2>;
  my @spl = split(" ", $q);
  die "Error! Read $spl[0] is duplicated in $ARGV[1]\n"
    if (exists $dup{$spl[0]});
  $dup{$spl[0]} = 1;

  # determine amplicon ID
  my $id = "";
  my $bot = 0;
  for (my $x = 1; $x < scalar @spl; $x++) {
    if ($spl[$x] eq "fwd" || $spl[$x] eq "rev") {
      $id = $spl[$x-1];
      $bot = 1 if ($x + 1 <= $#spl && $spl[$x+1] eq "both");
      last;
    }
  }
  die "Error! $ARGV[1] is improperly formatted:\n",
    "  no amplicon ID in $q\n" if (!$id);
  $count++;

  # check for duplicate
  if (exists $seq{$spl[0]}) {

    # skip chimeras
    if ($chim && $pr{$spl[0]} ne $id) {
      delete $seq{$spl[0]};
      $crem += 2;
      next;
    }

    # check if one has both primers removed
    if ($bth) {
      if ($both{$spl[0]} && !$bot) {
        print OUT $seq{$spl[0]};      
        delete $seq{$spl[0]};
        $print++;
        $brem++;
        next;
      } elsif ($bot && !$both{$spl[0]}) {
        print OUT "$q$w$e$r";
        delete $seq{$spl[0]};
        $print++;
        $brem++;
        next;
      }
    }

    # compare quality scores
    if ($qul) {
      chomp $r;
      my $tot = 0;
      for (my $x = 0; $x < length $r; $x++) {
        $tot += ord(substr($r, $x, 1)) - 33; 
      }
      # print higher avg -- if tie, choose 1st
      if ($tot / length $r > $qual{$spl[0]}) {
        print OUT "$q$w$e$r\n";
      } else {
        print OUT $seq{$spl[0]};
      }
      $print++;
      $qrem++;
    } else {
      # print both
      print OUT "$seq{$spl[0]}$q$w$e$r";
      $print += 2;
    }
    delete $seq{$spl[0]};

  } else {
    # not a duplicate
    print OUT "$q$w$e$r";
    $print++;
  }

}
close FQ2;

# print remaining singletons from 1st file
foreach my $re (sort keys %seq) {
  print OUT $seq{$re};
  $print++;
}
close OUT;

if ($verb) {
  print "Reads in $ARGV[1]: $count",
    "\nReads printed to $ARGV[2]: $print";
  print "\nReads removed for being chimeras: $crem" if ($chim);
  print "\nReads removed for not having both primers: $brem" if ($bth);
  print "\nReads removed for being lower quality: $qrem" if ($qul);
  print "\n";
}
