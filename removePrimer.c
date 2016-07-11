/*
  John Gaspar
  July 2014

  Removing primers from a fasta/fastq file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "removePrimer.h"

// global variables
static char* line;
static char* hline;
static Primer* primo;

/* void usage()
 * Prints usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./removePrimer {%s <file> %s <file>", INFILE, PRIMFILE);
  fprintf(stderr, " %s <file>} [optional parameters]\n", OUTFILE);
  fprintf(stderr, "Required parameters:\n");
  fprintf(stderr, "  %s  <file>       Input file containing reads from which to remove primers\n", INFILE);
  fprintf(stderr, "                     (in fasta or fastq format; can be gzip compressed,\n");
  fprintf(stderr, "                     with \"%s\" extension)\n", GZEXT);
  fprintf(stderr, "  %s  <file>       Input file listing primer sequences, one set (forward and\n", PRIMFILE);
  fprintf(stderr, "                     reverse) per line, comma- or tab-delimited. For example:\n");
  fprintf(stderr, "                       341F-926R,CCTACGGGAGGCAGCAG,AAACTCAAAKGAATTGACGG\n");
  fprintf(stderr, "                       968F-1401R,AACGCGAAGAACCTTAC,CGTTCCCGGGCCTTGTACACACCG\n");
  fprintf(stderr, "                     The sequences can contain IUPAC ambiguous DNA codes.\n");
  fprintf(stderr, "                     NOTE: Both primers should be given with respect to the plus\n");
  fprintf(stderr, "                       strand, i.e. the sequence given for the reverse primer is\n");
  fprintf(stderr, "                       the reverse-complement of actual reverse primer\n");
  fprintf(stderr, "  %s  <file>       Output file for trimmed reads (same format and compression\n", OUTFILE);
  fprintf(stderr, "                     as input file containing reads [%s])\n", INFILE);
  fprintf(stderr, "Optional parameters:\n");
  fprintf(stderr, "  %s <int[,int]>  Position (or range of positions) at which to begin searching\n", FWDPOS);
  fprintf(stderr, "                     for the first primer (def. 0 [i.e. search will begin\n");
  fprintf(stderr, "                     at the first base of each read]). For example, specifying\n");
  fprintf(stderr, "                     %s -1,1 will allow the first primer to match starting\n", FWDPOS);
  fprintf(stderr, "                     at positions -1, 0, or 1\n");
  fprintf(stderr, "  %s <int[,int]>  Position (or range of positions) at which to begin searching\n", REVPOS);
  fprintf(stderr, "                     for the second primer (def. 0 [i.e. search will begin\n");
  fprintf(stderr, "                     at the last base of each read])\n");
  fprintf(stderr, "  %s <int>        Mismatches to the first primer to allow (def. 0)\n", MISALLOW);
  fprintf(stderr, "  %s <int>        Mismatches to the second primer to allow for full-length,\n", REVMIS);
  fprintf(stderr, "                     3' end checking (def. 0)\n");
  fprintf(stderr, "  %s <int>        Check also for internal matches of the second primer to the\n", REVLENGTH);
  fprintf(stderr, "                     read, using the specified length of the primer\n");
  fprintf(stderr, "  %s <int>        Mismatches to the second primer to allow for internal\n", REVLMIS);
  fprintf(stderr, "                     matching (def. 0)\n");
  fprintf(stderr, "  %s  <file>       Check also for minimal matches of the second primer to the\n", BEDFILE);
  fprintf(stderr, "                     3' end of the read, using the expected amplicon lengths\n");
  fprintf(stderr, "                     derived from the given BED file (no mismatches allowed)\n");
  fprintf(stderr, "  %s <int[,int]>  Position (or range of positions) at which to begin searching\n", BEDPOS);
  fprintf(stderr, "                     for the second primer based on expected length (def. 0)\n");
  fprintf(stderr, "  %s              Option to require second primer be found\n", REVOPT);
  fprintf(stderr, "  %s  <file>       Log file for counts of matches\n", LOGFILE);
  fprintf(stderr, "  %s  <file>       Output file for non-trimmed reads\n", WASTEFILE);
  fprintf(stderr, "  %s  <file>       Output file for trimmed reads with correct primers reattached\n", CORRFILE);
  fprintf(stderr, "                     (should only be used if specifying %s)\n", REVOPT);
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
  else if (err == ERRPRIM) msg2 = MERRPRIM;
  else if (err == ERRPREP) msg2 = MERRPREP;
  else if (err == ERRINT) msg2 = MERRINT;
  else if (err == ERRBED) msg2 = MERRBED;
  else if (err == ERRBEDA) msg2 = MERRBEDA;
  else if (err == ERRINVAL) msg2 = MERRINVAL;
  else msg2 = DEFERR;

  fprintf(stderr, "Error! %s%s\n", msg, msg2);
  return -1;
}

/* void freeMemory()
 * Frees allocated memory.
 */
void freeMemory(void) {
  Primer* temp;
  for (Primer* p = primo; p != NULL; ) {
    free(p->name);
    free(p->fwd);
    free(p->rev);
    free(p->frc);
    free(p->rrc);
    temp = p;
    p = p->next;
    free(temp);
  }
  free(line);
  free(hline);
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

/* char* getLine()
 * Reads the next line from a file.
 */
char* getLine(char* line, int size, File in, int gz) {
  if (gz)
    return gzgets(in.gzf, line, size);
  else
    return fgets(line, size, in.f);
}

/* int fastaOrQ ()
 * Determines, based on the first character, if
 *   a file is likely fasta or fastq.
 */
int fastaOrQ(File in, int gz) {
  while (getLine(line, MAX_SIZE, in, gz) != NULL)
    if (line[0] == '>')
      return 0;
    else if (line[0] == '@')
      return 1;
    else if (line[0] != '#')
      break;
  exit(error("", ERRUNK));
}

/* int ambig(char, char)
 * Checks ambiguous DNA bases.
 */
int ambig(char x, char y) {
  if (x == 'N' ||
      (x == 'W' && (y == 'A' || y == 'T')) ||
      (x == 'S' && (y == 'C' || y == 'G')) ||
      (x == 'M' && (y == 'A' || y == 'C')) ||
      (x == 'K' && (y == 'G' || y == 'T')) ||
      (x == 'R' && (y == 'A' || y == 'G')) ||
      (x == 'Y' && (y == 'C' || y == 'T')) ||
      (x == 'B' && (y == 'C' || y == 'G' || y == 'T')) ||
      (x == 'D' && (y == 'A' || y == 'G' || y == 'T')) ||
      (x == 'H' && (y == 'A' || y == 'C' || y == 'T')) ||
      (x == 'V' && (y == 'A' || y == 'C' || y == 'G')))
    return 0;
  return 1;
}

/* Primer* findPrim(char*)
 * Finds a primer match to the given sequence.
 */
Primer* findPrim(char* seq, int misAllow, int fwdSt, int fwdEnd,
    int* st, int* f) {
  for (Primer* p = primo; p != NULL; p = p->next) {
    for (int i = 0; i < 2; i++) {
      char* prim = (i ? p->rrc : p->fwd);
      // allow primer to match starting at diff. positions
      for (int off = fwdSt; off < fwdEnd; off++) {
        int mis = misAllow;
        int j;
        for (j = 0; prim[j] != '\0'; j++)
          if (j + off < 0)
            continue;
          else if (seq[j+off] == '\0' ||
              (prim[j] != seq[j+off] &&
              (prim[j] == 'A' || prim[j] == 'C' ||
              prim[j] == 'G' || prim[j] == 'T' ||
              ambig(prim[j], seq[j+off])) && --mis < 0))
            break;

        if (prim[j] == '\0') {
          if (seq[j+off] != '\0') {
            *st = j+off;
            *f = i;
            return p;
          } else
            return NULL;
        }
      }
    }
  }
  return NULL;
}

/* int checkRevLen()
 * Checks the seq for a match of the reverse primer
 *   based on the expected amplicon length.
 */
int checkRevLen(char* seq, char* rev, int st,
    int bedSt, int bedEnd) {
  // check only the 3' fragment, do not allow mismatches
  // allow primer to match starting at diff. positions
  int len = strlen(seq);
  for (int off = bedSt; off < bedEnd; off++) {
    if (st+off >= len)
      break;
    int j;
    for (j = 0; rev[j] != '\0' && seq[st+j+off] != '\0'; j++)
      if (rev[j] != seq[st+j+off] &&
          (rev[j] == 'A' || rev[j] == 'C' ||
          rev[j] == 'G' || rev[j] == 'T' ||
          ambig(rev[j], seq[st+j+off])))
        break;
    if (rev[j] == '\0' || seq[st+j+off] == '\0')
      return st+off;
  }
  return 0;
}

/* int checkRevInt()
 * Checks the seq for a match of the reverse primer
 *   internally.
 */
int checkRevInt(char* seq, char* rev, int st,
    int misAllow, int len) {
  int last = strlen(seq) - len + 1;
  for (int i = st; i < last; i++) {
    int mis = misAllow;
    int j;
    for (j = 0; j < len; j++)
      if (rev[j] != seq[i + j] &&
          (rev[j] == 'A' || rev[j] == 'C' ||
          rev[j] == 'G' || rev[j] == 'T' ||
          ambig(rev[j], seq[i + j])) && --mis < 0)
        break;
    if (j == len)
      return i;
  }
  return 0;
}

/* int checkRevEnd()
 * Checks the seq for a match of the reverse primer
 *   at the 3' end.
 */
int checkRevEnd(char* seq, char* rev, int misAllow,
    int revSt, int revEnd) {
  int primEnd = strlen(rev) - 1;
  int seqEnd = strlen(seq) - 1;
  // allow primer to match starting at diff. positions
  for (int off = revSt; off < revEnd; off++) {
    int mis = misAllow;
    int j;
    for (j = 0; j < primEnd + 1; j++) {
      int seqPos = seqEnd - j - off;
      int pos = primEnd - j;
      if (j + off < 0)
        continue;
      else if (seqPos < 0 ||
          (rev[pos] != seq[seqPos] &&
          (rev[pos] == 'A' || rev[pos] == 'C' ||
          rev[pos] == 'G' || rev[pos] == 'T' ||
          ambig(rev[pos], seq[seqPos])) && --mis < 0))
        break;
    }
    if (j == primEnd + 1) {
      return seqEnd - j - off + 1;  // last base of primer
    }
  }
  return 0;
}

/* int readFile()
 * Parses the input file. Produces the output file(s).
 */
int readFile(File in, File out, int misAllow, int* match,
    int* rcmatch, int fwdSt, int fwdEnd, int revSt, int revEnd,
    int bedSt, int bedEnd, File waste, int wasteOpt, int revMis,
    int revLen, int revLMis, int revOpt, File corr, int corrOpt,
    int gz) {
  // determine if input is fasta or fastq
  int aorq = fastaOrQ(in, gz);
  gz ? gzrewind(in.gzf) : rewind(in.f);

  int count = 0;
  while (getLine(hline, MAX_SIZE, in, gz) != NULL) {
    if (hline[0] == '#')
      continue;
    count++;
    if (getLine(line, MAX_SIZE, in, gz) == NULL)
      exit(error("", ERRSEQ));
    int len = strlen(line) - 1;
    if (line[len] == '\n')
      line[len] = '\0';

    int st = 0, end = 0, f = 0;
    Primer* p = findPrim(line, misAllow, fwdSt, fwdEnd, &st, &f);
    if (p != NULL) {
      (*match)++;
      f ? p->rcount++ : p->fcount++;

      // search for reverse primer
      // first, check 3' end
      char* rev = (f ? p->frc : p->rev);
      end = checkRevEnd(line, rev, revMis, revSt, revEnd);

      // check internal sequence
      if (!end && revLen) {
        int setLen = strlen(rev);
        if (setLen > revLen)
          setLen = revLen;
        end = checkRevInt(line, rev, st, revLMis, setLen);
      }

      // check based on amplicon length
      if (!end && p->len && st + p->len < strlen(line))
        end = checkRevLen(line, rev, st + p->len, bedSt, bedEnd);

      // evaluate outcome, produce output
      if (end <= st)
        end = 0;
      if (end)
        f ? p->rcountr++ : p->fcountr++;
      if (revOpt && !end) {
        // rev primer not found (and was required [revOpt])
        if (wasteOpt)
          gz ? gzprintf(waste.gzf, "%s%s\n", hline, line)
            : fprintf(waste.f, "%s%s\n", hline, line);
      } else {
        // print header
        for (int i = 0; hline[i] != '\0' && hline[i] != '\n'; i++)
          gz ? gzputc(out.gzf, hline[i]) : putc(hline[i], out.f);
        gz ? gzprintf(out.gzf, " %s%s%s\n", p->name,
          f ? REV : FWD, end ? BOTH : "")
          : fprintf(out.f, " %s%s%s\n", p->name,
          f ? REV : FWD, end ? BOTH : "");
        if (corrOpt) {
          for (int i = 0; hline[i] != '\0' && hline[i] != '\n'; i++)
            gz ? gzputc(corr.gzf, hline[i]) : putc(hline[i], corr.f);
          gz ? gzprintf(corr.gzf, " %s%s%s\n", p->name,
            f ? REV : FWD, end ? BOTH : "")
            : fprintf(corr.f, " %s%s%s\n", p->name,
            f ? REV : FWD, end ? BOTH : "");
        }
        // print sequence
        if (!end)
          end = len;
        else
          (*rcmatch)++;
        for (int i = st; i < end; i++)
          gz ? gzputc(out.gzf, line[i]) : putc(line[i], out.f);
        gz ? gzputc(out.gzf, '\n') : putc('\n', out.f);
        // reattach primers
        if (corrOpt) {
          gz ? gzprintf(corr.gzf, "%s", f ? p->rrc : p->fwd)
            : fprintf(corr.f, "%s", f ? p->rrc : p->fwd);
          for (int i = st; i < end; i++)
            gz ? gzputc(corr.gzf, line[i]) : putc(line[i], corr.f);
          gz ? gzprintf(corr.gzf, "%s\n", f ? p->frc : p->rev)
            : fprintf(corr.f, "%s\n", f ? p->frc : p->rev);
        }
      }
    } else if (wasteOpt)
      gz ? gzprintf(waste.gzf, "%s%s\n", hline, line)
        : fprintf(waste.f, "%s%s\n", hline, line);

    // read next 2 lines if fastq
    if (aorq) {
      for (int i = 0; i < 2; i++)
        if (getLine(line, MAX_SIZE, in, gz) == NULL)
          exit(error("", ERRSEQ));
        else if (p != NULL) {
          if (revOpt && !end) {
            if (wasteOpt)
              gz ? gzprintf(waste.gzf, "%s", line)
                : fprintf(waste.f, "%s", line);
          } else if (i) {
            for (int j = st; j < end; j++)
              gz ? gzputc(out.gzf, line[j]) : putc(line[j], out.f);
            gz ? gzputc(out.gzf, '\n') : putc('\n', out.f);
            if (corrOpt) {
              for (int j = 0; j < strlen(f ? p->rrc : p->fwd); j++)
                gz ? gzputc(corr.gzf, 'I') : putc('I', corr.f);
              for (int j = st; j < end; j++)
                gz ? gzputc(corr.gzf, line[j]) : putc(line[j], corr.f);
              for (int j = 0; j < strlen(f ? p->frc : p->rev); j++)
                gz ? gzputc(corr.gzf, 'I') : putc('I', corr.f);
              gz ? gzputc(corr.gzf, '\n') : putc('\n', corr.f);
            }
          } else {
            gz ? gzprintf(out.gzf, "%s", line)
              : fprintf(out.f, "%s", line);
            if (corrOpt)
              gz ? gzprintf(corr.gzf, "%s", line)
                : fprintf(corr.f, "%s", line);
          }
        } else if (wasteOpt)
          gz ? gzprintf(waste.gzf, "%s", line)
            : fprintf(waste.f, "%s", line);

    }
  }
  return count;
}

/* void getPos()
 * Determines the start-end positions for the primer search.
 */
void getPos(char* pos, int* start, int* end) {
  if (pos == NULL)
    return;

  char* beg = strtok(pos, CSV);
  if (beg == NULL)
    exit(error(beg, ERRINT));
  *start = getInt(beg);

  char* stop = strtok(NULL, END);
  *end = (stop != NULL ? getInt(stop) + 1 : *start + 1);
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

/* void openGZWrite()
 * Open a (possibly gzip compressed) file for writing.
 */
void openGZWrite(char* outFile, File* out, int gz) {
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
  } else
    out->f = openWrite(outFile);
}

/* FILE* openRead()
 * Opens a file for reading.
 */
FILE* openRead(char* inFile) {
  FILE* in = fopen(inFile, "r");
  if (in == NULL)
    exit(error(inFile, ERROPEN));
  return in;
}

/* void openFiles()
 * Opens the files to run the program.
 */
void openFiles(char* outFile, File* out,
    char* primFile, FILE** prim, char* inFile, File* in,
    char* logFile, FILE** log, char* bedFile, FILE** bed,
    char* wasteFile, File* waste,
    char* corrFile, File* corr, int gz) {
  // open required files
  *prim = openRead(primFile);
  if (gz) {
    // gzip compressed files
    in->gzf = gzopen(inFile, "r");
    if (in->gzf == NULL)
      exit(error(inFile, ERROPEN));
    openGZWrite(outFile, out, gz);
  } else {
    in->f = openRead(inFile);
    out->f = openWrite(outFile);
  }

  // open optional files
  if (bedFile != NULL)
    *bed = openRead(bedFile);
  if (logFile != NULL)
    *log = openWrite(logFile);
  if (wasteFile != NULL)
    openGZWrite(wasteFile, waste, gz);
  if (corrFile != NULL)
    openGZWrite(corrFile, corr, gz);
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
  else if (in == 'Y') out = 'R';
  else if (in == 'R') out = 'Y';
  else if (in == 'W') out = 'W';
  else if (in == 'S') out = 'S';
  else if (in == 'K') out = 'M';
  else if (in == 'M') out = 'K';
  else if (in == 'B') out = 'V';
  else if (in == 'V') out = 'B';
  else if (in == 'D') out = 'H';
  else if (in == 'H') out = 'D';
  else if (in == 'N') out = 'N';
  else exit(error("", ERRPRIM));
  return out;
}

/* char* revComp(char*)
 * Reverse-complements the given sequence.
 */
char* revComp(char* seq) {
  int i = strlen(seq) - 1;
  char* out = (char*) memalloc(2 + i);
  int j;
  for (j = 0; i > -1; j++) {
    char nuc = seq[i--];
    out[j] = rc(nuc);
  }
  out[j] = '\0';
  return out;
}

/* int loadSeqs(FILE*)
 * Loads the primers from the given file.
 */
int loadSeqs(FILE* prim) {

  Primer* prev = NULL;
  int count = 0;
  while (fgets(line, MAX_SIZE, prim) != NULL) {

    if (line[0] == '#')
      continue;

    // load name and sequence
    char* name = strtok(line, CSV);
    char* seq = strtok(NULL, CSV);
    char* rev = strtok(NULL, DEL);
    if (name == NULL || seq == NULL || rev == NULL) {
      error("", ERRPRIM);
      continue;
    }

    // check for duplicate
    for (Primer* pc = primo; pc != NULL; pc = pc->next)
      if (!strcmp(pc->name, name))
        exit(error(name, ERRPREP));

    // create primer
    Primer* p = (Primer*) memalloc(sizeof(Primer));
    p->name = (char*) memalloc(1 + strlen(name));
    p->fwd = (char*) memalloc(1 + strlen(seq));
    p->rev = (char*) memalloc(1 + strlen(rev));
    strcpy(p->name, name);
    strcpy(p->fwd, seq);
    strcpy(p->rev, rev);

    // save sequence rc's
    p->frc = revComp(p->fwd);
    p->rrc = revComp(p->rev);

    p->fcount = p->rcount = p->fcountr = p->rcountr = 0;
    p->next = NULL;
    if (primo == NULL)
      primo = p;
    else
      prev->next = p;
    prev = p;
    p->len = 0;
    p->fpos = -1;
    count++;
  }

  return count;
}

/* void getLengths()
 * Determine expected lengths of amplicons.
 */
void getLengths(FILE* bed) {
  while (fgets(line, MAX_SIZE, bed) != NULL) {
    if (line[0] == '#')
      continue;

    // load positions
    char* first = strtok(line, CSV);
    first = strtok(NULL, CSV);
    char* second = strtok(NULL, CSV);
    char* amp = strtok(NULL, DEL);
    if (first == NULL || second == NULL || amp == NULL) {
      error("", ERRBED);
      continue;
    }
    int firstPos = getInt(first);
    int secondPos = getInt(second);

    // find amplicon
    Primer* p;
    for (p = primo; p != NULL; p = p->next)
      if (!strcmp(p->name, amp))
        break;
    if (p == NULL || p->len)
      continue;

    // save length
    if (p->fpos != -1) {
      p->len = (p->fpos < firstPos ? firstPos - p->rpos :
        p->fpos - secondPos);
      if (p->len < 0) {
        error(amp, ERRBEDA);
        p->len = 0;
      }
    } else {
      p->fpos = firstPos;
      p->rpos = secondPos;
    }
  }
}

/* void getParams()
 * Parses the command line.
 */
void getParams(int argc, char** argv) {

  char* outFile = NULL, *inFile = NULL, *primFile = NULL,
    *bedFile = NULL, *fwdPos = NULL, *revPos = NULL,
    *bedPos = NULL, *logFile = NULL, *wasteFile = NULL,
    *corrFile = NULL;
  int misAllow = 0, revLen = 0, revMis = 0, revLMis = 0,
    revOpt = 0;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (!strcmp(argv[i], REVOPT))
      revOpt = 1;
    else if (i < argc - 1) {
      if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], INFILE))
        inFile = argv[++i];
      else if (!strcmp(argv[i], PRIMFILE))
        primFile = argv[++i];
      else if (!strcmp(argv[i], BEDFILE))
        bedFile = argv[++i];
      else if (!strcmp(argv[i], FWDPOS))
        fwdPos = argv[++i];
      else if (!strcmp(argv[i], REVPOS))
        revPos = argv[++i];
      else if (!strcmp(argv[i], BEDPOS))
        bedPos = argv[++i];
      else if (!strcmp(argv[i], LOGFILE))
        logFile = argv[++i];
      else if (!strcmp(argv[i], WASTEFILE))
        wasteFile = argv[++i];
      else if (!strcmp(argv[i], MISALLOW))
        misAllow = getInt(argv[++i]);
      else if (!strcmp(argv[i], REVMIS))
        revMis = getInt(argv[++i]);
      else if (!strcmp(argv[i], REVLENGTH))
        revLen = getInt(argv[++i]);
      else if (!strcmp(argv[i], REVLMIS))
        revLMis = getInt(argv[++i]);
      else if (!strcmp(argv[i], CORRFILE))
        corrFile = argv[++i];
      else
        exit(error(argv[i], ERRINVAL));
    } else
      exit(error(argv[i], ERRINVAL));
  }

  if (outFile == NULL || inFile == NULL || primFile == NULL)
    usage();
  int gz = 0;
  if (!strcmp(inFile + strlen(inFile) - strlen(GZEXT), GZEXT))
    gz = 1;

  // open files, load primer sequences
  File out, in, waste, corr;
  FILE* prim = NULL, *log = NULL, *bed = NULL;
  openFiles(outFile, &out, primFile, &prim, inFile, &in,
    logFile, &log, bedFile, &bed, wasteFile, &waste,
    corrFile, &corr, gz);
  int pr = loadSeqs(prim);

  // get start and end locations
  int fwdSt = 0, fwdEnd = 1, revSt = 0, revEnd = 1,
    bedSt = 0, bedEnd = 1;
  getPos(fwdPos, &fwdSt, &fwdEnd);
  getPos(revPos, &revSt, &revEnd);
  if (bed != NULL) {
    getLengths(bed);
    getPos(bedPos, &bedSt, &bedEnd);
  }

  // read file
  int match = 0, rcmatch = 0;  // counting variables
  int count = readFile(in, out, misAllow, &match, &rcmatch,
    fwdSt, fwdEnd, revSt, revEnd, bedSt, bedEnd,
    waste, wasteFile != NULL, revMis, revLen, revLMis,
    revOpt, corr, corrFile != NULL, gz);

  // print log output
  if (log != NULL) {
    fprintf(log, "Primer pairs: %d\nRead count: %d\n", pr, count);
    fprintf(log, "Primer matches: %d\nBoth primer matches: %d\n\n", match, rcmatch);
    fprintf(log, "Matches:\nPrimer\tFwd\tFwd-Both\tRev\tRev-Both\n");
    for (Primer* p = primo; p != NULL; p = p->next)
      fprintf(log, "%s\t%d\t%d\t%d\t%d\n", p->name, p->fcount, p->fcountr,
        p->rcount, p->rcountr);
  }

  // close files
  if ( (gz && (gzclose(in.gzf) != Z_OK || gzclose(out.gzf) != Z_OK ||
      (wasteFile != NULL && gzclose(waste.gzf) != Z_OK) ||
      (corrFile != NULL && gzclose(corr.gzf) != Z_OK) ) ) ||
      ( ! gz && (fclose(in.f) || fclose(out.f) ||
      (wasteFile != NULL && fclose(waste.f)) ||
      (corrFile != NULL && fclose(corr.f)) ) ) ||
      fclose(prim) || (log != NULL && fclose(log)) ||
      (bed != NULL && fclose(bed)) )
    exit(error("", ERRCLOSE));
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  line = (char*) memalloc(MAX_SIZE);
  hline = (char*) memalloc(MAX_SIZE);
  primo = NULL;
  getParams(argc, argv);
  freeMemory();
  return 0;
}
