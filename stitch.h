/*
  John Gaspar
  April 2015

  Header file for stitch.c.
*/

#define MAX_SIZE    1024   // maximum length for input line
#define NOTMATCH    1.5f   // stitch failure

// command-line parameters
#define HELP        "-h"
#define FIRST       "-1"
#define SECOND      "-2"
#define OUTFILE     "-o"
#define UNFILE1     "-u1"
#define UNFILE2     "-u2"
#define LOGFILE     "-l"
#define OVERLAP     "-m"
#define MISMATCH    "-p"
#define DOVEOPT     "-d"
#define DOVEFILE    "-dl"
#define MAXOPT      "-n"
#define VERBOSE     "-ve"

// default parameter values
#define DEFOVER     20
#define DEFMISM     0.0f

// third parameter to copyStr()
#define FWD         0
#define RC          1
#define REV         2

// error messages
#define ERROPEN     0
#define MERROPEN    ": cannot open file for reading"
#define ERRCLOSE    1
#define MERRCLOSE   "Cannot close file"
#define ERROPENW    2
#define MERROPENW   ": cannot open file for writing"
#define ERRUNK      3
#define MERRUNK     "Unknown nucleotide"
#define ERRMEM      4
#define MERRMEM     "Cannot allocate memory"
#define ERRSEQ      5
#define MERRSEQ     "Cannot load sequence"
#define ERRQUAL     6
#define MERRQUAL    "Sequence/quality scores do not match"
#define ERRHEAD     7
#define MERRHEAD    ": not matched in input files"
#define ERRINT      8
#define MERRINT     ": cannot convert to int"
#define ERRFLOAT    9
#define MERRFLOAT   ": cannot convert to float"
#define ERRPARAM    10
#define MERRPARAM   ": unknown command-line parameter"
#define ERROVER     11
#define MERROVER    "Overlap must be greater than 0"
#define ERRMISM     12
#define MERRMISM    "Mismatch must be in [0,1)"
#define DEFERR      "Unknown error"
