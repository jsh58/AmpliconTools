/*
  John Gaspar
  April 2015

  Stitching paired-end reads together.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stitch.h"

/* void usage()
 * Prints usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./stitch {%s <file> %s <file>", FIRST, SECOND);
  fprintf(stderr, " %s <file>} [optional parameters]\n", OUTFILE);
  fprintf(stderr, "Required parameters:\n");
  fprintf(stderr, "  %s  <file>       Input file containing reads from forward direction\n", FIRST);
  fprintf(stderr, "  %s  <file>       Input file containing reads from reverse direction\n", SECOND);
  fprintf(stderr, "                   NOTE: Reads in both input files can be trimmed of poor\n");
  fprintf(stderr, "                     quality bases at the 3' end, but, since the stitched\n");
  fprintf(stderr, "                     read is defined by the 5' ends of the PE reads, one\n");
  fprintf(stderr, "                     should be wary of trimming them at that end.\n");
  fprintf(stderr, "  %s  <file>       Output file for stitched reads\n", OUTFILE);
  fprintf(stderr, "Optional parameters:\n");
  fprintf(stderr, "  %s  <file>       Log file for stitching results\n", LOGFILE);
  fprintf(stderr, "  %s <file>       File containing non-stitched forward reads\n", UNFILE1);
  fprintf(stderr, "  %s <file>       File containing non-stitched reverse reads\n", UNFILE2);
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
int getSeq(FILE* in, char* line, char* seq,
    char* qual, int nSeq, int nQual) {
  for (int i = 0; i < 3; i++) {
    if (fgets(line, MAX_SIZE, in) == NULL)
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
void printRes(FILE* out, FILE* log, FILE* dove,
    char* header, char* seq1, char* seq2, char* qual1,
    char* qual2, int len1, int len2, int pos, float best) {
  // log result
  if (log != NULL) {
    fprintf(log, "%s\t%d\t%d\t", header,
      pos < 0 ? (len2+pos < len1 ? len2+pos : len1) :
      (len1-pos < len2 ? len1-pos : len2),
      len2 + pos);
    best ? fprintf(log, "%.3f", best) : fprintf(log, "0");
    fprintf(log, "\n");
  }

  // log 3' overhangs of dovetailed sequence(s)
  if (dove != NULL && (len1 > len2 + pos || pos < 0)) {
    fprintf(dove, "%s\t%s\t", header, len1 > len2 + pos ?
      seq1 + len2 + pos : "-");
    if (pos < 0)
      for (int i = -1; i - pos > -1; i--)
        fprintf(dove, "%c", rc(seq2[i - pos]));
    else
      fprintf(dove, "-");
    fprintf(dove, "\n");
  }

  // print stitched sequence
  createSeq(seq1, seq2, qual1, qual2, len1, len2, pos);
  fprintf(out, "@%s\n%s\n+\n%s\n", header, seq1, qual1);
}

/* void printFail()
 * Print stitch failure reads.
 */
void printFail(FILE* un1, FILE* un2, FILE* log,
    char* header, char* head1, char* head2, char* seq1,
    char* seq2, char* qual1, char* qual2, int len) {
  if (log != NULL)
    fprintf(log, "%s\tn/a\n", header);
  if (un1 != NULL && un2 != NULL) {
    fprintf(un1, "@%s\n%s\n+\n%s\n", head1, seq1, qual1);
    // put rev sequence back
    fprintf(un2, "@%s\n", head2);
    for (int i = len - 1; i > -1; i--)
      fprintf(un2, "%c", rc(seq2[i]));
    fprintf(un2, "\n+\n");
    for (int i = len - 1; i > -1; i--)
      fprintf(un2, "%c", qual2[i]);
    fprintf(un2, "\n");
  }
}

/* int readFile()
 * Parses the input file. Produces the output file(s).
 */
int readFile(FILE* in1, FILE* in2, FILE* out,
    FILE* un1, FILE* un2, FILE* log, int overlap,
    int dovetail, FILE* dove, float mismatch,
    int maxLen, int* stitch, int* fail) {

  char* line = (char*) memalloc(MAX_SIZE);
  char* head1 = (char*) memalloc(MAX_SIZE);
  char* seq1 = (char*) memalloc(MAX_SIZE);
  char* qual1 = (char*) memalloc(MAX_SIZE);
  char* head2 = (char*) memalloc(MAX_SIZE);
  char* seq2 = (char*) memalloc(MAX_SIZE);
  char* qual2 = (char*) memalloc(MAX_SIZE);
  char* header = (char*) memalloc(MAX_SIZE); // consensus header

  int count = 0;
  while (fgets(line, MAX_SIZE, in1) != NULL) {
    count++;

    // save headers
    int i;
    for (i = 0; line[i + 1] != '\n'; i++)
      head1[i] = line[i + 1];
    head1[i] = '\0';
    if (fgets(line, MAX_SIZE, in2) == NULL)
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
        else
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
    int len1 = getSeq(in1, line, seq1, qual1, FWD, FWD);
    int len2 = getSeq(in2, line, seq2, qual2, RC, REV);

    // stitch reads, print result
    float best = 1.0f;
    int pos = findPos(seq1, seq2, qual1, qual2, len1, len2,
      overlap, dovetail, mismatch, maxLen, &best);
    if (pos == len1 - overlap + 1) {
      printFail(un1, un2, log, header, head1, head2, seq1,
        seq2, qual1, qual2, len2);
      (*fail)++;
    } else {
      printRes(out, log, dove, header, seq1, seq2, qual1,
        qual2, len1, len2, pos, best);
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

/* FILE* openWrite()
 * Opens a file for writing.
 */
FILE* openWrite(char* outFile) {
  FILE* out = fopen(outFile, "w");
  if (out == NULL)
    exit(error(outFile, ERROPENW));
  return out;
}

/* void openFiles()
 * Opens the files to run the program.
 */
void openFiles(char* outFile, FILE** out,
    char* inFile1, FILE** in1, char* inFile2, FILE** in2,
    char* unFile1, FILE** un1, char* unFile2, FILE** un2,
    char* logFile, FILE** log, char* doveFile, FILE** dove,
    int dovetail) {
  // open required files
  *out = openWrite(outFile);
  *in1 = fopen(inFile1, "r");
  if (*in1 == NULL)
    exit(error(inFile1, ERROPEN));
  *in2 = fopen(inFile2, "r");
  if (*in2 == NULL)
    exit(error(inFile2, ERROPEN));

  // open optional files
  if (unFile1 != NULL && unFile2 != NULL) {
    *un1 = openWrite(unFile1);
    *un2 = openWrite(unFile2);
  }
  if (logFile != NULL) {
    *log = openWrite(logFile);
    fprintf(*log, "Read\tOverlapLen\tStitchedLen\tMismatch\n");
  }
  if (dovetail && doveFile != NULL) {
    *dove = openWrite(doveFile);
    fprintf(*dove, "Read\tDovetailFwd\tDovetailRev\n");
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

  // open files
  FILE* out = NULL, *in1 = NULL, *in2 = NULL, *un1 = NULL,
    *un2 = NULL, *log = NULL, *dove = NULL;
  openFiles(outFile, &out, inFile1, &in1, inFile2, &in2,
    unFile1, &un1, unFile2, &un2, logFile, &log,
    doveFile, &dove, dovetail);

  // read file
  int stitch = 0, fail = 0;  // counting variables
  int count = readFile(in1, in2, out, un1, un2, log,
    overlap, dovetail, dove, mismatch, maxLen,
    &stitch, &fail);

  if (verbose) {
    printf("Reads analyzed: %d\n", count);
    printf("  Successfully stitched: %d\n", stitch);
    printf("  Stitch failures: %d\n", fail);
  }

  // close files
  if (fclose(out) || fclose(in1) || fclose(in2) ||
      (un1 != NULL && fclose(un1)) || (un2 != NULL && fclose(un2)) ||
      (log != NULL && fclose(log)) ||
      (dovetail && doveFile != NULL && fclose(dove)))
    exit(error("", ERRCLOSE));
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getParams(argc, argv);
  return 0;
}
