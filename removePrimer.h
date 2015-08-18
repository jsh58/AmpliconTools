/*
  John Gaspar
  July 2014

  Header file for removePrimer.c.
*/

#define MAX_SIZE    1024    // maximum length for input line
#define CSV         ",\t"
#define DEL         ",\t\n"
#define END         "\0"

// output labels
#define FWD         " fwd"  // read matched a fwd primer
#define REV         " rev"  // read matched a rev primer
#define BOTH        " both" // read matched primers on both ends

// command-line options and parameters
#define HELP        "-h"
#define INFILE      "-i"
#define OUTFILE     "-o"
#define PRIMFILE    "-p"
#define FWDPOS      "-fp"
#define REVPOS      "-rp"
#define MISALLOW    "-ef"
#define REVMIS      "-er"
#define REVLENGTH   "-rl"
#define REVLMIS     "-el"
#define BEDFILE     "-b"
#define BEDPOS      "-bp"
#define REVOPT      "-rq"
#define LOGFILE     "-l"
#define WASTEFILE   "-w"
#define CORRFILE    "-c"

// error messages
#define ERROPEN     0
#define MERROPEN    "cannot open file for reading"
#define ERRCLOSE    1
#define MERRCLOSE   "cannot close file"
#define ERROPENW    2
#define MERROPENW   "cannot open file for writing"
#define ERRUNK      3
#define MERRUNK     "unknown file type (not fasta or fastq)"
#define ERRMEM      4
#define MERRMEM     "cannot allocate memory"
#define ERRSEQ      5
#define MERRSEQ     "cannot load sequence"
#define ERRPRIM     6
#define MERRPRIM    "cannot load primer sequence"
#define ERRPREP     7
#define MERRPREP    ": cannot repeat primer name"
#define ERRINT      8
#define MERRINT     ": cannot convert to int"
#define ERRBED      9
#define MERRBED     "cannot load value from BED file"
#define ERRBEDA     10
#define MERRBEDA    ": error determining length from BED file"
#define ERRINVAL    11
#define MERRINVAL   ": invalid parameter or usage"
#define DEFERR      "Unknown error"

typedef struct primer {
  char* name;
  char* fwd;
  char* rev;
  char* frc;
  char* rrc;
  int len;  // expected amplicon length
  int fpos;
  int rpos;
  int fcount;
  int rcount;
  int fcountr;
  int rcountr;
  struct primer* next;
} Primer;
