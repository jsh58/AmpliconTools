#!/bin/bash -e

# John M. Gaspar
# Dec. 2014

# This script is designed to detect variants in a sample that has
#   been analyzed by amplicon-based targeted resequencing.

# External software requirements:
#   - bowtie2 (2.2.3)   -- assumed to be in $PATH
#   - samtools (0.1.19) -- assumed to be in $PATH
#   - VarScan (2.3.7)   -- set location here:
VARSCAN="./VarScan.v2.3.7.jar"
if [ ! -f $VARSCAN ]; then
  echo "VarScan jar file not found"
  exit -1
fi

# check command-line arguments
if [ $# -lt 6 ]; then
  echo "Usage: ./`basename $0`  <FASTQ1>  <FASTQ2>  <BED>  <GEN>  <IDX>  <DIR>" 1>&2
  exit -1
fi

# input files
file1=$1          # Input FASTQ file, paired-end reads #1
file2=$2          # Input FASTQ file, paired-end reads #2
bed=$3            # BED file listing locations of primers
gen=$4            # reference genome (FASTA)
idx=$5            # bowtie2 index prefix (indexes will be generated if necessary)
dir=$6            # output directory

# check input files
if [[ ! -f $file1 || ! -f $file2 ]]; then
  echo "Input FASTQ files not found"
  exit -1
elif [ ! -f $bed ]; then
  echo "Input primer BED file not found"
  exit -1
elif [ ! -f $gen ]; then
  echo "Input reference genome not found"
  exit -1
fi

# determine if input FASTQ files are gzip compressed
gz=""
if [[ ${file1:(-3)} == ".gz" &&
    ${file2:(-3)} == ".gz" ]]; then
  gz=".gz"
fi

# set home directory
HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# retrieve primer-target sequences
prim=primers.txt
if [ ! -f $prim ]; then
  perl ${HOME_DIR}/getPrimers.pl $bed $gen $prim
fi

# stitch together reads
echo "Stitching reads"
tr1=join.fastq$gz
tr2=un1.fastq$gz
tr3=un2.fastq$gz
stParam="-m 20 -p 0.1 -d"  # min overlap 20, 10% allowed mismatches, dovetailing
${HOME_DIR}/stitch -1 $file1 -2 $file2 -o $tr1 -u1 $tr2 -u2 $tr3 $stParam

# remove primers, with -rq
echo "Removing primers"
log1=joinlog.txt
tr0=join-pr.fastq$gz
tr4=join-nopr.fastq$gz
rpParam="-fp -1,1 -rp -1,1 -ef 2 -er 2"  # allowing 2 subs, can start at +/- 1
${HOME_DIR}/removePrimer -i $tr1 -p $prim -o $tr0 $rpParam -rq -l $log1 -w $tr4  # require both primers

# retrieve reads whose primers weren't found
echo "Getting failure reads"
tr5=nopr1.fastq$gz
tr6=nopr2.fastq$gz
python ${HOME_DIR}/getReads.py $tr4 $file1 $tr5
python ${HOME_DIR}/getReads.py $tr4 $file2 $tr6
# cat with unjoined reads
cat $tr2 >> $tr5
cat $tr3 >> $tr6

# removePrimer individually
echo "Removing primers individually"
tr7=nopr1-pr.fastq$gz
tr8=nopr2-pr.fastq$gz
log2=nopr1log.txt
log3=nopr2log.txt
rpParam2="-rl 16 -el 1 -b $bed -bp -1,1"  # more options to find second primer
${HOME_DIR}/removePrimer -i $tr5 -p $prim -o $tr7 $rpParam $rpParam2 -l $log2
${HOME_DIR}/removePrimer -i $tr6 -p $prim -o $tr8 $rpParam $rpParam2 -l $log3

# filter singletons
tr9=noprcomb.fastq$gz
fsParam="-b -q -c"  # prefer both primers removed, higher quality read, no chimeras
perl ${HOME_DIR}/filterSingle.pl $tr7 $tr8 $tr9 $fsParam

# quality trim
echo "Quality filtering"
out1=joined.fastq$gz
qtParam="-t 30 -n 20"  # min avg qual 30; min len 20; no window filtering
${HOME_DIR}/qualTrim -i $tr0 -o $out1 $qtParam
tr10=noprcomb-qt.fastq$gz
${HOME_DIR}/qualTrim -i $tr9 -o $tr10 $qtParam

# cat joined and singletons
out2=combined.fastq$gz
cat $out1 $tr10 > $out2

# check for bowtie2 indexes
gen2=tempGenome.fa
if [[ ! -f $idx.1.bt2 || ! -f $idx.2.bt2 ||
    ! -f $idx.3.bt2 || ! -f $idx.4.bt2 ||
    ! -f $idx.rev.1.bt2 || ! -f $idx.rev.2.bt2 ]]; then
  echo "Building bowtie2 indexes"
  if [ ${gen:(-3)} == ".gz" ]; then
    echo "  Gunzipping genome"
    gunzip -c $gen > $gen2
    bowtie2-build $gen2 $idx
  else
    bowtie2-build $gen $idx
  fi
fi

# map with bowtie2
echo "Mapping with bowtie2"
out3=combined.sam
bwtParam="-D 200 -N 1 -L 18 -i S,1,0.50 -k 20"
proc=1   # number of processors
bowtie2 -x $idx -U $out2 -S $out3 $bwtParam -p $proc

# find length variants
echo "Finding length variants"
tr11=len1.txt
perl ${HOME_DIR}/findLengthVars.pl $out1 $bed $tr11
# the output file (len1.txt) can be edited to exclude certain length variants
tr12=len2.fastq
perl ${HOME_DIR}/getLengthVars.pl $out1 $tr11 $tr12
len3=realign.txt
tr13=len3v.txt
perl ${HOME_DIR}/alignLengthVars.pl $tr12 $prim $bed $len3 $tr13 $gen

# check alternative mapping sites
out4=altMapping.txt
if [ ! -f $out4 ]; then
  echo "Checking alternative mapping sites"
  perl ${HOME_DIR}/checkAltMapping.pl $out2 $out3 $prim $bed $gen $out4
fi

# at this point, one can combine the altMapping results for
#   multiple samples using combAltMapping.pl, then proceed
#   with the new altMapping.txt as $out4

# filter SAM -- multi-mapping and realignment of length variants
echo "Filtering SAM"
out5=combinedFiltered.sam
log4=realign.log
perl ${HOME_DIR}/filterSAM.pl $out2 $bed $out4 $out3 $out5 $len3 $log4

# convert SAM to sorted BAM
echo "Converting SAM to sorted BAM"
out6=combinedFiltered.bam
samtools view -bS $out5 | samtools sort -f - $out6

# call variants
echo "Calling variants"
out7=combinedFiltered.pileup
if [ ${gen:(-3)} == ".gz" ]; then
  if [ ! -f $gen2 ]; then
    echo "  Gunzipping genome for samtools mpileup"
    gunzip -c $gen > $gen2
  fi
  samtools mpileup -B -Q 0 -d 100000 -f $gen2 $out6 > $out7
else
  samtools mpileup -B -Q 0 -d 100000 -f $gen $out6 > $out7
fi
qual=30  # min quality score to count a base
tr15=combFil.snp
java -jar ${VARSCAN} pileup2snp --min-avg-qual $qual --min-coverage 0 --min-var-freq 0.01 --variants < $out7 > $tr15
tr16=combFil.indel
java -jar ${VARSCAN} pileup2indel --min-avg-qual $qual --min-coverage 0 --min-var-freq 0.01 --variants < $out7 > $tr16

# make VCF, filter
echo "Producing VCF"
tr17=temp.vcf
perl ${HOME_DIR}/makeVCF.pl $tr15 $tr16 $out7 $tr17 $qual $dir
tr18=temp2.vcf
perl ${HOME_DIR}/addHPtoVCF.pl $tr17 $gen $tr18

# filter variants: remove variants outside target regions
out8=combinedFiltered.vcf
perl ${HOME_DIR}/filterVCF.pl $tr18 $out8 -b $bed
  # other filtering options for filterVCF.pl:
  #   min. depth 20                          -d 20
  #   min. variant allele observations 10    -o 10
  #   min. allele freq. (AF) 0.02            -a 0.02
  #   min. AF for in/dels 0.05               -i 0.05
  #   min. AF for C:G>T:A 0.1                -m 0.1
  #   min. AF for variants in/near homopolymer runs (sliding scale):
  #       -p 4,0.05,0.1:5,0.1,0.2:6,0.2,0.3:7,0.3,0.4:8,0.4,0.5

# remove extra files
rm $tr0 $tr1 $tr2 $tr3 $tr4 $tr5 $tr6 $tr7 $tr8 $tr9 $tr10 \
  $tr11 $tr12 $tr13 $tr15 $tr16 $tr17 $tr18
if [ -f $gen2 ]; then
  rm $gen2
fi

# move files to output directory
if [ ! -d $dir ]; then
  mkdir $dir
fi
if [ $dir != "." ]; then
  mv $out1 $out2 $out3 $out4 $out5 $out6 $out7 $out8 \
    $log1 $log2 $log3 $log4 $len3 $dir
fi
