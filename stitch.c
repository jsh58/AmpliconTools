/*
  John Gaspar
  April 2015

  Stitching paired-end reads together.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "stitch.h"

/* void usage()
 * Prints usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./stitch {%s <file> %s <file>", FIRST, SECOND);
  fprintf(stderr, " %s <file>} [optional parameters]\n", OUTFILE);
  fprintf(stderr, "Required parameters:\n");
  fprintf(stderr, "  %s  <file>       Input FASTQ file with reads from forward direction\n", FIRST);
  fprintf(stderr, "  %s  <file>       Input FASTQ file with reads from reverse direction\n", SECOND);
  fprintf(stderr, "  %s  <file>       Output FASTQ file for stitched reads\n", OUTFILE);
  fprintf(stderr, "  Note: Both input files can be gzip compressed (with \"%s\"\n", GZEXT);
  fprintf(stderr, "    extensions), in which case the output FASTQ file(s) will also\n");
  fprintf(stderr, "    be gzip compressed.\n");
  fprintf(stderr, "    Also, reads in both input files can be trimmed of poor quality\n");
  fprintf(stderr, "    bases prior to using this program, but, since the stitched read\n");
  fprintf(stderr, "    is defined by the 5' ends of the PE reads, one should be wary\n");
  fprintf(stderr, "    of trimming them at that end.\n");
  fprintf(stderr, "Optional parameters:\n");
  fprintf(stderr, "  %s  <file>       Log file for stitching results\n", LOGFILE);
  fprintf(stderr, "  %s <file>       FASTQ file containing non-stitched forward reads\n", UNFILE1);
  fprintf(stderr, "  %s <file>       FASTQ file containing non-stitched reverse reads\n", UNFILE2);
  fprintf(stderr, "  %s  <int>        Minimum overlap of the paired-end reads (def. 20)\n", OVERLAP);
  fprintf(stderr, "  %s  <float>      Mismatches to allow in the overlapped region\n", MISMATCH);
  fprintf(stderr, "                     (in [0-1), a fraction of the overlap length; def. 0)\n");
  fprintf(stderr, "  %s               Option to check for dovetailing of the reads\n", DOVEOPT);
  fprintf(stderr, "                     (min. overlap from %s still applies)\n", OVERLAP);
  fprintf(stderr, "  %s <file>       Log file for dovetailed reads only\n", DOVEFILE);
  fprintf(stderr, "  %s               Option to produce shortest stitched read, given\n", MAXOPT);
  fprintf(stderr, "                     multiple overlapping possibilities (by default,\n");
  fprintf(stderr, "                     the longest stitched read is produced)\n");
  fprintf(stderr, "  %s              Option to print counts of stitching results to stdout\n", VERBOSE);
  exit(-1);
}

/* int error()
 * Prints an error message.
 */
int error(char* msg, int err) {
  char* msg2;
  if (err == ERROPEN) msg2 = MERROPEN;
  else if (err == ERRCLOSE) msg2 = MERRCLOSE;
  else if (err == ERROPENW) msg2 = MERROPENW;
  else if (err == ERRUNK) msg2 = MERRUNK;
  else if (err == ERRMEM) msg2 = MERRMEM;
  else if (err == ERRSEQ) msg2 = MERRSEQ;
  else if (err == ERRQUAL) msg2 = MERRQUAL;
  else if (err == ERRHEAD) msg2 = MERRHEAD;
  else if (err == ERRINT) msg2 = MERRINT;
  else if (err == ERRFLOAT) msg2 = MERRFLOAT;
  else if (err == ERRPARAM) msg2 = MERRPARAM;
  else if (err == ERROVER) msg2 = MERROVER;
  else if (err == ERRMISM) msg2 = MERRMISM;
  else msg2 = DEFERR;

  fprintf(stderr, "Error! %s%s\n", msg, msg2);
  return -1;
}

/* void* memalloc()
 * Allocates a heap block.
 */
void* memalloc(int size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* float getFloat(char*)
 * Converts the given char* to a float.
 */
float getFloat(char* in) {
  char** endptr = NULL;
  float ans = strtof(in, endptr);
  if (endptr != '\0')
    exit(error(in, ERRFLOAT));
  return ans;
}

/* int getInt(char*)
 * Converts the given char* to an int.
 */
int getInt(char* in) {
  char** endptr = NULL;
  int ans = (int) strtol(in, endptr, 10);
  if (endptr != '\0')
    exit(error(in, ERRINT));
  return ans;
}

/* char rc(char)
 * Returns the complement of the given base.
 */
char rc(char in) {
  char out;
  if (in == 'A') out = 'T';
  else if (in == 'T') out = 'A';
  else if (in == 'C') out = 'G';
  else if (in == 'G') out = 'C';
  else if (in == 'N') out = 'N';
  else exit(error("", ERRUNK));
  return out;
}

/* char* getLine()
 * Reads the next line from a file.
 */
char* getLine(char* line, int size, File in, int gz) {
  if (gz)
    return gzgets(in.gzf, line, size);
  else
    return fgets(line, size, in.f);
}

/* void copyStr()
 * Copy a sequence/quality score.
 * Reverse (REV) or rev-comp (RC) if needed (3rd param).
 */
void copyStr(char* out, char* in, int rev) {
  int i = 0;
  if (rev == FWD)
    for (; in[i] != '\n' && in[i] != '\0'; i++)
      out[i] = in[i];
  else {
    int j = 0;
    for (; in[j] != '\n' && in[j] != '\0'; j++) ;
    for (j--; j > -1; j--)
      out[i++] = (rev == REV ? in[j] : rc(in[j]));
  }
  out[i] = '\0';
}

/* int getSeq()
 * Read sequence and quality scores from a fastq file.
 */
int getSeq(File in, char* line, char* seq, char* qual,
    int nSeq, int nQual, int gz) {
  for (int i = 0; i < 3; i++) {
    if (getLine(line, MAX_SIZE, in, gz) == NULL)
      exit(error("", ERRSEQ));
    if (i == 0)
      copyStr(seq, line, nSeq);
    else if (i == 2)
      copyStr(qual, line, nQual);
  }
  int len = strlen(seq);
  if (len != strlen(qual))
    exit(error("", ERRQUAL));
  return len;
}

/* float compare()
 * Compare two sequences. Return the percent mismatch.
 */
float compare(char* seq1, char* seq2, int length,
    float mismatch, int overlap) {
  int mis = 0;       // number of mismatches
  int len = length;  // length of overlap, not counting Ns
  float allow = len * mismatch;
  for (int i = 0; i < length; i++) {
    // do not count Ns
    if (seq1[i] == 'N' || seq2[i] == 'N') {
      if (--len < overlap || mis > len * mismatch)
        return NOTMATCH;
      allow = len * mismatch;
    } else if (seq1[i] != seq2[i] && ++mis > allow)
      return NOTMATCH;
  }
  return (float) mis / len;
}

/* int findPos()
 * Find optimal overlapping position.
 */
int findPos (char* seq1, char* seq2, char* qual1,
    char* qual2, int len1, int len2, int overlap,
    int dovetail, float mismatch, int maxLen,
    float* best) {
  int pos = len1 - overlap + 1;  // position of match
  for (int i = len1 - overlap; i > -1; i--) {
    if (len1 - i > len2 && !dovetail)
      break;
    float res = compare(seq1 + i, seq2,
      len1-i < len2 ? len1-i : len2, mismatch, overlap);
    if (res < *best || (res == *best && !maxLen)) {
      *best = res;
      pos = i;
    }
    if (res == 0.0f && maxLen)
      return pos;  // shortcut for exact match
  }

  // check for dovetailing
  if (dovetail) {
    for (int i = 1; i < len2 - overlap + 1; i++) {
      float res = compare(seq1, seq2 + i,
        len2-i < len1 ? len2-i : len1, mismatch, overlap);
      if (res < *best || (res == *best && !maxLen)) {
        *best = res;
        pos = -i;
      }
      if (res == 0.0f && maxLen)
        return pos;  // shortcut for exact match
    }
  }

  return pos;
}

/* void createSeq()
 * Create stitched sequence (into seq1, qual1).
 */
void createSeq(char* seq1, char* seq2, char* qual1, char* qual2,
    int len1, int len2, int pos) {
  int len = len2 + pos;  // length of stitched sequence
  for (int i = 0; i < len; i++) {
    if (i - pos < 0)
      continue;
    // disagreements favor higher quality score or
    //   equal quality score that is closer to 5' end
    else if (i >= len1 ||
        (seq1[i] != seq2[i-pos] && (qual1[i] < qual2[i-pos] ||
        (qual1[i] == qual2[i-pos] && i >= len2 - i + pos)))) {
      seq1[i] = seq2[i-pos];
      qual1[i] = qual2[i-pos];
    } else if (qual1[i] < qual2[i-pos])
      qual1[i] = qual2[i-pos];
  }
  seq1[len] = '\0';
  qual1[len] = '\0';
}

/* void printRes()
 * Print stitched read.
 */
void printRes(File out, File log, int logOpt, File dove,
    int doveOpt, char* header, char* seq1, char* seq2,
    char* qual1, char* qual2, int len1, int len2,
    int pos, float best, int gz) {
  // log result
  if (logOpt) {
    fprintf(log.f, "%s\t%d\t%d\t", header,
      pos < 0 ? (len2+pos < len1 ? len2+pos : len1) :
      (len1-pos < len2 ? len1-pos : len2), len2 + pos);
    best ? fprintf(log.f, "%.3f", best) : fprintf(log.f, "0");
    fprintf(log.f, "\n");
  }

  // log 3' overhangs of dovetailed sequence(s)
  if (doveOpt && (len1 > len2 + pos || pos < 0)) {
    fprintf(dove.f, "%s\t%s\t", header, len1 > len2 + pos ?
      seq1 + len2 + pos : "-");
    if (pos < 0)
      for (int i = -1; i - pos > -1; i--)
        fprintf(dove.f, "%c", rc(seq2[i - pos]));
    else
      fprintf(dove.f, "-");
    fprintf(dove.f, "\n");
  }

  // print stitched sequence
  createSeq(seq1, seq2, qual1, qual2, len1, len2, pos);
  gz ? gzprintf(out.gzf, "@%s\n%s\n+\n%s\n", header, seq1, qual1)
    : fprintf(out.f, "@%s\n%s\n+\n%s\n", header, seq1, qual1);
}

/* void printFail()
 * Print stitch failure reads.
 */
void printFail(File un1, File un2, int unOpt,
    File log, int logOpt, char* header, char* head1,
    char* head2, char* seq1, char* seq2, char* qual1,
    char* qual2, int len, int gz) {
  if (logOpt)
    fprintf(log.f, "%s\tn/a\n", header);
  if (unOpt) {
    gz ? gzprintf(un1.gzf, "@%s\n%s\n+\n%s\n", head1, seq1, qual1)
      : fprintf(un1.f, "@%s\n%s\n+\n%s\n", head1, seq1, qual1);
    // put rev sequence back
    gz ? gzprintf(un2.gzf, "@%s\n", head2)
      : fprintf(un2.f, "@%s\n", head2);
    for (int i = len - 1; i > -1; i--)
      gz ? gzputc(un2.gzf, rc(seq2[i])) : putc(rc(seq2[i]), un2.f);
    gz ? gzprintf(un2.gzf, "\n+\n") : fprintf(un2.f, "\n+\n");
    for (int i = len - 1; i > -1; i--)
      gz ? gzputc(un2.gzf, qual2[i]) : putc(qual2[i], un2.f);
    gz ? gzputc(un2.gzf, '\n') : putc('\n', un2.f);
  }
}

/* int readFile()
 * Parses the input file. Produces the output file(s).
 */
int readFile(File in1, File in2, File out,
    File un1, File un2, int unOpt, File log,
    int logOpt, int overlap, int dovetail, File dove,
    int doveOpt, float mismatch, int maxLen, int* stitch,
    int* fail, int gz) {

  char* line = (char*) memalloc(MAX_SIZE);
  char* head1 = (char*) memalloc(MAX_SIZE);
  char* seq1 = (char*) memalloc(MAX_SIZE);
  char* qual1 = (char*) memalloc(MAX_SIZE);
  char* head2 = (char*) memalloc(MAX_SIZE);
  char* seq2 = (char*) memalloc(MAX_SIZE);
  char* qual2 = (char*) memalloc(MAX_SIZE);
  char* header = (char*) memalloc(MAX_SIZE); // consensus header

  int count = 0;
  while (getLine(line, MAX_SIZE, in1, gz) != NULL) {
    count++;

    // save headers
    int i;
    for (i = 0; line[i + 1] != '\n'; i++)
      head1[i] = line[i + 1];
    head1[i] = '\0';
    if (getLine(line, MAX_SIZE, in2, gz) == NULL)
      exit(error("", ERRSEQ));
    for (i = 0; line[i + 1] != '\n'; i++)
      head2[i] = line[i + 1];
    head2[i] = '\0';

    // make sure headers match (up to first space character),
    // save consensus header too
    int ok = 0;
    int j;
    for (j = 0; j < i; j++) {
      if (head1[j] != head2[j]) {
        if (ok)
          break;
        exit(error(head1, ERRHEAD));
      } else if (head1[j] == ' ')
        ok = 1;  // headers match
      header[j] = head1[j];
    }
    if (header[j - 1] == ' ')
      header[j - 1] = '\0'; // removing trailing space
    else
      header[j] = '\0';

    // save sequences and quality scores for the reads
    int len1 = getSeq(in1, line, seq1, qual1, FWD, FWD, gz);
    int len2 = getSeq(in2, line, seq2, qual2, RC, REV, gz);

    // stitch reads, print result
    float best = 1.0f;
    int pos = findPos(seq1, seq2, qual1, qual2, len1, len2,
      overlap, dovetail, mismatch, maxLen, &best);
    if (pos == len1 - overlap + 1) {
      printFail(un1, un2, unOpt, log, logOpt, header, head1,
        head2, seq1, seq2, qual1, qual2, len2, gz);
      (*fail)++;
    } else {
      printRes(out, log, logOpt, dove, doveOpt, header, seq1,
        seq2, qual1, qual2, len1, len2, pos, best, gz);
      (*stitch)++;
    }
  }

  // free memory
  free(line);
  free(head1);
  free(seq1);
  free(qual1);
  free(head2);
  free(seq2);
  free(qual2);
  free(header);
  return count;
}

/* void openWrite()
 * Open a file for writing.
 */
void openWrite(char* outFile, File* out, int gz) {
  if (gz) {
    if (!strcmp(outFile + strlen(outFile) - strlen(GZEXT), GZEXT))
      out->gzf = gzopen(outFile, "w");
    else {
      // add ".gz" to outFile
      char* outFile2 = memalloc(strlen(outFile) +
        strlen(GZEXT) + 1);
      strcpy(outFile2, outFile);
      strcat(outFile2, GZEXT);
      out->gzf = gzopen(outFile2, "w");
      free(outFile2);
    }
    if (out->gzf == NULL)
      exit(error(outFile, ERROPENW));
  } else {
    out->f = fopen(outFile, "w");
    if (out->f == NULL)
      exit(error(outFile, ERROPENW));
  }
}

/* void openRead()
 * Open a file for reading.
 */
void openRead(char* inFile, File* in, int gz) {
  if (gz) {
    in->gzf = gzopen(inFile, "r");
    if (in->gzf == NULL)
      exit(error(inFile, ERROPEN));
  } else {
    in->f = fopen(inFile, "r");
    if (in->f == NULL)
      exit(error(inFile, ERROPEN));
  }
}

/* void openFiles()
 * Opens the files to run the program.
 */
void openFiles(char* outFile, File* out,
    char* inFile1, File* in1, char* inFile2,
    File* in2, char* unFile1, File* un1,
    char* unFile2, File* un2, char* logFile,
    File* log, char* doveFile, File* dove,
    int dovetail, int gz) {
  // open required files
  openRead(inFile1, in1, gz);
  openRead(inFile2, in2, gz);
  openWrite(outFile, out, gz);

  // open optional files
  if (unFile1 != NULL && unFile2 != NULL) {
    openWrite(unFile1, un1, gz);
    openWrite(unFile2, un2, gz);
  }
  if (logFile != NULL) {
    openWrite(logFile, log, 0);
    fprintf(log->f, "Read\tOverlapLen\tStitchedLen\tMismatch\n");
  }
  if (dovetail && doveFile != NULL) {
    openWrite(doveFile, dove, 0);
    fprintf(dove->f, "Read\tDovetailFwd\tDovetailRev\n");
  }
}

/* void getParams()
 * Parses the command line.
 */
void getParams(int argc, char** argv) {

  char* outFile = NULL, *inFile1 = NULL, *inFile2 = NULL,
    *unFile1 = NULL, *unFile2 = NULL, *logFile = NULL,
    *doveFile = NULL;
  int overlap = DEFOVER, dovetail = 0, maxLen = 1;
  int verbose = 0;
  float mismatch = DEFMISM;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (!strcmp(argv[i], MAXOPT))
      maxLen = 0;
    else if (!strcmp(argv[i], DOVEOPT))
      dovetail = 1;
    else if (!strcmp(argv[i], VERBOSE))
      verbose = 1;
    else if (i < argc - 1) {
      if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], FIRST))
        inFile1 = argv[++i];
      else if (!strcmp(argv[i], SECOND))
        inFile2 = argv[++i];
      else if (!strcmp(argv[i], UNFILE1))
        unFile1 = argv[++i];
      else if (!strcmp(argv[i], UNFILE2))
        unFile2 = argv[++i];
      else if (!strcmp(argv[i], LOGFILE))
        logFile = argv[++i];
      else if (!strcmp(argv[i], DOVEFILE))
        doveFile = argv[++i];
      else if (!strcmp(argv[i], OVERLAP))
        overlap = getInt(argv[++i]);
      else if (!strcmp(argv[i], MISMATCH))
        mismatch = getFloat(argv[++i]);
      else
        exit(error(argv[i], ERRPARAM));
    } else
      usage();
  }

  // check for parameter errors
  if (outFile == NULL || inFile1 == NULL || inFile2 == NULL)
    usage();
  if (overlap <= 0)
    exit(error("", ERROVER));
  if (mismatch < 0.0f || mismatch >= 1.0f)
    exit(error("", ERRMISM));

  // determine if inputs are gzip compressed
  int gz = 0;
  if (!strcmp(inFile1 + strlen(inFile1) - strlen(GZEXT), GZEXT) &&
      !strcmp(inFile2 + strlen(inFile2) - strlen(GZEXT), GZEXT) )
    gz = 1;

  // open files
  File out, in1, in2, un1, un2, log, dove;
  openFiles(outFile, &out, inFile1, &in1, inFile2, &in2,
    unFile1, &un1, unFile2, &un2, logFile, &log,
    doveFile, &dove, dovetail, gz);

  // read file
  int stitch = 0, fail = 0;  // counting variables
  int count = readFile(in1, in2, out, un1, un2,
    unFile1 != NULL && unFile2 != NULL, log, logFile != NULL,
    overlap, dovetail, dove, dovetail && doveFile != NULL,
    mismatch, maxLen, &stitch, &fail, gz);

  if (verbose) {
    printf("Reads analyzed: %d\n", count);
    printf("  Successfully stitched: %d\n", stitch);
    printf("  Stitch failures: %d\n", fail);
  }

  // close files
  if ( ( gz && ( gzclose(in1.gzf) != Z_OK ||
      gzclose(in2.gzf) != Z_OK || gzclose(out.gzf) != Z_OK ||
      (unFile1 != NULL && unFile2 != NULL &&
      (gzclose(un1.gzf) != Z_OK || gzclose(un2.gzf) != Z_OK) ) ) ) ||
      ( ! gz && ( fclose(in1.f) || fclose(in2.f) || fclose(out.f) ||
      (unFile1 != NULL && unFile2 != NULL &&
      (fclose(un1.f) || fclose(un2.f)) ) ) ) ||
      (logFile != NULL && fclose(log.f)) ||
      (dovetail && doveFile != NULL && fclose(dove.f)) )
    exit(error("", ERRCLOSE));
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getParams(argc, argv);
  return 0;
}
