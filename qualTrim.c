/*
  John Gaspar
  July 2014

  Quality trimming a fastq file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "qualTrim.h"

/* void usage()
 * Print usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./qualTrim {%s <file> ", INFILE);
  fprintf(stderr, "%s <file>} [optional parameters]\n", OUTFILE);
  fprintf(stderr, "Required parameters:\n");
  fprintf(stderr, "  %s <file>   Input FASTQ file with Sanger-scale quality\n", INFILE);
  fprintf(stderr, "                scores (phred + 33). Can be gzip compressed,\n");
  fprintf(stderr, "                with \"%s\" extension\n", GZEXT);
  fprintf(stderr, "  %s <file>   Output FASTQ file (will be gzip compressed if\n", OUTFILE);
  fprintf(stderr, "                input is)\n");
  fprintf(stderr, "Optional parameters:\n");
  fprintf(stderr, "  %s <int>    Window length\n", WINDOWLEN);
  fprintf(stderr, "  %s <float>  Minimum avg. quality in the window\n", WINDOWAVG);
  fprintf(stderr, "  %s <float>  Minimum avg. quality for the full read\n", QUALAVG);
  fprintf(stderr, "                (after any window truncations)\n");
  fprintf(stderr, "  %s <int>    Minimum length of a read\n", MINLEN);
  fprintf(stderr, "                (after any window truncations)\n");
  fprintf(stderr, "  %s          Option to trim reads only at 5' end\n", FIVEOPT);
  fprintf(stderr, "  %s          Option to trim reads only at 3' end\n", THREEOPT);
  fprintf(stderr, "  %s         Option to print counts of results to stdout\n", VERBOSE);
  exit(-1);
}

/* int error()
 * Print an error message.
 */
int error(char* msg, int err) {
  char* msg2;
  if (err == ERROPEN) msg2 = MERROPEN;
  else if (err == ERRCLOSE) msg2 = MERRCLOSE;
  else if (err == ERROPENW) msg2 = MERROPENW;
  else if (err == ERRMEM) msg2 = MERRMEM;
  else if (err == ERRPARAM) msg2 = MERRPARAM;
  else if (err == ERRSEQ) msg2 = MERRSEQ;
  else if (err == ERRINT) msg2 = MERRINT;
  else if (err == ERRFLOAT) msg2 = MERRFLOAT;
  else msg2 = DEFERR;

  fprintf(stderr, "Error! %s%s\n", msg, msg2);
  return -1;
}

/* void memalloc()
 * Allocate memory from the heap.
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

/* int getEnd()
 * Determine the 3' end.
 */
int getEnd(char* line, int len, float qual, int end) {
  int last = 0;
  int i;
  for (i = 1; i < len; i++) {
    float sum = 0.0f;
    int j = end - 1;
    for (int k = 0; k < i; k++)
      sum += line[j - k] - OFFSET;
    if (sum / i < qual)
      last = j - i + 1;
    else if (last)
      return last;
  }
  for (i = end - 1; i > len - 2; i--) {
    float sum = 0.0f;
    for (int j = 0; j < len; j++)
      sum += line[i - j] - OFFSET;
    if (sum / len < qual)
      last = i - len + 1;
    else if (last)
      return last;
    else
      break;
  }
  return i == len - 2 ? last : end;
}

/* int getStart()
 * Determine the 5' end.
 */
int getStart(char* line, int len, float qual, int end) {
  int st = 0;
  int i;
  for (i = 1; i < len; i++) {
    float sum = 0.0f;
    for (int k = 0; k < i; k++)
      sum += line[k] - OFFSET;
    if (sum / i < qual)
      st = i;
    else if (st)
      return st;
  }
  for (i = 0; i < end - len + 1; i++) {
    float sum = 0.0f;
    for (int j = 0; j < len; j++)
      sum += line[i + j] - OFFSET;
    if (sum / len < qual)
      st = i + len;
    else if (st)
      return st;
    else
      break;
  }
  return st;
}

/* int trimQual()
 * Determine ends of the read.
 */
int trimQual(char* line, int len, float qual, int* end,
    int opt5, int opt3) {
  if (opt3)
    *end = getEnd(line, len, qual, *end);
  int st = 0;
  if (opt5)
    st = getStart(line, len, qual, *end);
  return st;
}

/* int checkQual()
 * Check average of all qual scores (from st to end).
 * Return 1 if OK, else 0.
 */
int checkQual(char* line, int st, int end, float avg) {
  float sum = 0.0f;
  for (int i = st; i < end; i++)
    sum += line[i] - OFFSET;
  return sum / (end - st) < avg ? 1 : 0;
}

/* void printOut()
 * Print the given string to a file.
 */
void printOut(union File out, const char* format,
    char* str, int gz) {
  if (gz)
    gzprintf(out.gzf, format, str);
  else
    fprintf(out.f, format, str);
}

/* char* getLine()
 * Reads the next line from a file.
 */
char* getLine(char* line, int size, union File in, int gz) {
  if (gz)
    return gzgets(in.gzf, line, size);
  else
    return fgets(line, size, in.f);
}

/* void readFile()
 * Control the I/O.
 */
void readFile(union File in, union File out, int len,
    float qual, float avg, int minLen, int opt5, int opt3,
    int gz, int verbose) {
  char* head = (char*) memalloc(MAX_SIZE);
  char* seq = (char*) memalloc(MAX_SIZE);
  char* line = (char*) memalloc(MAX_SIZE);

  int count = 0, elim = 0;
  while (getLine(head, MAX_SIZE, in, gz) != NULL) {
    if (head[0] != '@')
      continue;

    // load sequence and quality scores
    if (getLine(seq, MAX_SIZE, in, gz) == NULL)
      exit(error("", ERRSEQ));
    for (int i = 0; i < 2; i++)
      if (getLine(line, MAX_SIZE, in, gz) == NULL)
        exit(error("", ERRSEQ));
    int end = strlen(line) - 1;
    if (line[end] == '\n')
      line[end] = '\0';
    if (end < len) {
      elim++;
      continue;
    }

    // trim read
    int st = 0;
    if (len)
      st = trimQual(line, len, qual, &end, opt5, opt3);
    if (avg && checkQual(line, st, end, avg)) {
      elim++;
      continue;
    }

    // print output
    if (st < end && end - st >= minLen) {
      printOut(out, "%s", head, gz);
      for (int i = st; i < end; i++)
        gz ? gzputc(out.gzf, seq[i]) : putc(seq[i], out.f);
      printOut(out, "\n+\n", NULL, gz);
      for (int i = st; i < end; i++)
        gz ? gzputc(out.gzf, line[i]) : putc(line[i], out.f);
      printOut(out, "\n", NULL, gz);
      count++;
    } else
      elim++;
  }

  if (verbose)
    printf("Reads printed: %d\nReads eliminated: %d\n",
      count, elim);

  free(seq);
  free(line);
  free(head);
}

/* void openWrite()
 * Open a file for writing.
 */
void openWrite(char* outFile, union File* out, int gz) {
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

/* void openFiles()
 * Open input and output files.
 */
void openFiles(char* outFile, union File* out,
    char* inFile, union File* in, int gz) {
  if (gz) {
    in->gzf = gzopen(inFile, "r");
    if (in->gzf == NULL)
      exit(error(inFile, ERROPEN));
  } else {
    in->f = fopen(inFile, "r");
    if (in->f == NULL)
      exit(error(inFile, ERROPEN));
  }
  openWrite(outFile, out, gz);
}

/* void getParams()
 * Get command-line parameters.
 */
void getParams(int argc, char** argv) {

  char* outFile = NULL, *inFile = NULL;
  int windowLen = 0, minLen = 0, opt5 = 1, opt3 = 1;
  int verbose = 0;
  float windowAvg = 0.0f, qualAvg = 0.0f;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (!strcmp(argv[i], FIVEOPT))
      opt3 = 0;
    else if (!strcmp(argv[i], THREEOPT))
      opt5 = 0;
    else if (!strcmp(argv[i], VERBOSE))
      verbose = 1;
    else if (i < argc - 1) {
      if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], INFILE))
        inFile = argv[++i];
      else if (!strcmp(argv[i], WINDOWLEN))
        windowLen = getInt(argv[++i]);
      else if (!strcmp(argv[i], WINDOWAVG))
        windowAvg = getFloat(argv[++i]);
      else if (!strcmp(argv[i], QUALAVG))
        qualAvg = getFloat(argv[++i]);
      else if (!strcmp(argv[i], MINLEN))
        minLen = getInt(argv[++i]);
      else
        exit(error(argv[i], ERRPARAM));
    } else
      exit(error(argv[i], ERRPARAM));
  }

  if (outFile == NULL || inFile == NULL)
    usage();

  // process file
  union File out, in;
  int gz = 0;
  if (!strcmp(inFile + strlen(inFile) - strlen(GZEXT), GZEXT))
    gz = 1;
  openFiles(outFile, &out, inFile, &in, gz);
  readFile(in, out, windowLen, windowAvg, qualAvg,
    minLen, opt5, opt3, gz, verbose);

  if ( (gz && (gzclose(in.gzf) != Z_OK || gzclose(out.gzf) != Z_OK))
      || ( ! gz && (fclose(in.f) || fclose(out.f))) )
    exit(error("", ERRCLOSE));
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getParams(argc, argv);
  return 0;
}
