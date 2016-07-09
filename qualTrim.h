/*
  John Gaspar
  July 2014

  Header file for qualTrim.c.
*/

#define MAX_SIZE    1024   // maximum length for input line
#define OFFSET      33     // ASCII-based offset of quality scores
#define GZEXT       ".gz"  // file extension for gzip compression

// command-line parameters
#define HELP        "-h"
#define INFILE      "-i"
#define OUTFILE     "-o"
#define WINDOWLEN   "-l"   // length for sliding window of quality scores
#define WINDOWAVG   "-q"   // average quality score for sliding window
#define QUALAVG     "-t"   // average quality score for entire read
#define MINLEN      "-n"   // minimum length of a read
#define FIVEOPT     "-5"   // option to trim only at 5' end
#define THREEOPT    "-3"   // option to trim only at 3' end
#define VERBOSE     "-ve"  // option to print counts to stdout

// error messages
#define ERROPEN     0
#define MERROPEN    ": cannot open file for reading"
#define ERRCLOSE    1
#define MERRCLOSE   "Cannot close file"
#define ERROPENW    2
#define MERROPENW   ": cannot open file for writing"
#define ERRMEM      3
#define MERRMEM     "Cannot allocate memory"
#define ERRPARAM    4
#define MERRPARAM   ": unknown command-line parameter"
#define ERRSEQ      5
#define MERRSEQ     "Cannot load sequence"
#define ERRFLOAT    6
#define MERRFLOAT   ": cannot convert to float"
#define ERRINT      7
#define MERRINT     ": cannot convert to int"
#define DEFERR      "Unknown error"

union File {
  FILE* f;
  gzFile gzf;
};
