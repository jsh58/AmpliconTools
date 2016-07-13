#!/usr/bin/python

# John M. Gaspar
# July 2016

# Retrieve a set of reads from a fastq file.

import sys
import gzip

def usage():
  print 'Usage: python getReads.py  <infile1>  <infile2>  <outfile>'
  print '  Required:'
  print '    <infile1>  FASTQ file listing reads to be retrieved'
  print '                 (used only to get read headers)'
  print '    <infile2>  FASTQ file listing original reads'
  print '    <outfile>  Output file for reads'
  print '  Note: Input FASTQ files may be gzip compressed (with ".gz"'
  print '    extension). <outfile> will be compressed if <infile2> is.'
  sys.exit(-1)

def openRead(filename):
  '''
  Open the given file for reading. Determine if file
    is gzip compressed based on extension ".gz".
  '''
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
      gz = 1
    else:
      f = open(filename, 'r')
      gz = 0
  except IOError:
    sys.stderr.write('Error! Cannot open %s\n' % filename)
    sys.exit(-1)
  return f, gz

def openWrite(filename, gz):
  '''
  Open given file for writing, gzip compressed or not.
  '''
  try:
    if gz:
      if filename[-3:] != '.gz':
        filename += '.gz'
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

def loadHeaders(f):
  '''
  Load fastq headers from given file.
  '''
  d = {}
  while True:
    line = f.readline()
    if not line: break
    if line[0] != '@': continue
    spl = line.rstrip().split()
    d[spl[0]] = 1
    for i in range(3):
      line = f.readline()
  return d

def processReads(fIn, fOut, headers):
  '''
  Write reads matching headers from fIn to fOut.
  '''
  while True:
    line = fIn.readline()
    if not line: break
    if line[0] != '@': continue
    spl = line.rstrip().split()
    if spl[0] in headers:
      fOut.write(line)
      for i in range(3):
        line = fIn.readline()
        fOut.write(line)
    else:
      for i in range(3):
        line = fIn.readline()

def main():
  '''
  Main.
  '''
  args = sys.argv[1:]
  if len(args) < 3 or args[0] == '-h':
    usage()

  # open files
  f1, dummy = openRead(args[0])
  f2, gz = openRead(args[1])
  fOut = openWrite(args[2], gz)

  # save read headers from f1
  headers = loadHeaders(f1)
  f1.close()

  # copy matching reads from f2
  processReads(f2, fOut, headers)
  f2.close()
  fOut.close()

if __name__ == '__main__':
  main()
