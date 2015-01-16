:Begin:
:Function: readform_file
:Pattern: ReadForm[filename_String]
:Arguments: {filename}
:ArgumentTypes: {Manual}
:ReturnType: Manual
:End:

:Begin:
:Function: readform_exec
:Pattern: ReadForm[formcmd_String, filename_String]
:Arguments: {formcmd, filename}
:ArgumentTypes: {Manual}
:ReturnType: Manual
:End:

:Evaluate: _ReadForm := (Message[ReadForm::syntax]; Abort[])

:Begin:
:Function: readformdebug
:Pattern: ReadFormDebug[debug_Integer, filename_:""]
:Arguments: {debug, filename}
:ArgumentTypes: {Integer, Manual}
:ReturnType: Manual
:End:

:Evaluate: ReadForm::syntax = "Bad syntax."

:Evaluate: ReadForm::noopen = "Cannot open ``."

:Evaluate: ReadForm::nooutput =
  "Something went wrong, there was no output from FORM."

:Evaluate: ReadForm::formerror = "``"


/*
	ReadForm.tm
		reads FORM output back into Mathematica
		this file is part of FormCalc
		last modified 15 Feb 13 th

Note: FORM code must have
	1. #- (no listing),
	2. off stats,
	3. should produce output with print (not print +s).

Debug:
	bit 0 = DEB_IN  = listing of input
	bit 1 = DEB_OUT = transfer of results
	bit 2 = DEB_CHU = internal chains before OrderChain
	bit 3 = DEB_CHO = internal chains after OrderChain
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <assert.h>
#include <signal.h>
#include <sys/wait.h>
#include <sys/types.h>

#define TERMBUF 500000

#define DEB_IN 1
#define DEB_OUT 2
#define DEB_CHU 4
#define DEB_CHO 8

#define RED "\e[31m"
#define BLUE "\e[34m"
#define RESET "\e[0m"


/*

A term in the FORM output is organized into the TERM structure
in the following way:

 ____4_____     __3___     __2___     ___0____     _____1_____
/          \   /      \   /      \   /        \   /           \
SumOver(...) * Mat(...) * Den(...) * paveM(...) * ..... * (...)

Hierarchy of collecting:
4. SumOver, PowerOf
3. Mat
2. Den
1. [coefficient]
0. paveM, cutM

*/

#define LEVEL_LOOP 0
#define LEVEL_COEFF 1
#define LEVEL_DEN 2
#define LEVEL_MAT 3
#define LEVEL_SUMOVER 4
#define NLEVELS 5

typedef const int cint;
typedef const size_t csize_t;
typedef unsigned char byte;
typedef byte *string;
typedef MLCONST byte *cstring;

typedef struct {
  cstring name;
  int level;
} FUN;

static const FUN funtab[] = {
  {(cstring)"SumOver", LEVEL_SUMOVER},
  {(cstring)"PowerOf", LEVEL_SUMOVER},
  {(cstring)"Mat",     LEVEL_MAT},
  {(cstring)"Den",     LEVEL_DEN},
  {(cstring)"helM",    LEVEL_DEN},
  {(cstring)"paveM",   LEVEL_LOOP},
  {(cstring)"cutM",    LEVEL_LOOP},
  {(cstring)"A0i",     LEVEL_LOOP},
  {(cstring)"B0i",     LEVEL_LOOP},
  {(cstring)"C0i",     LEVEL_LOOP},
  {(cstring)"D0i",     LEVEL_LOOP},
  {(cstring)"E0i",     LEVEL_LOOP},
  {(cstring)"F0i",     LEVEL_LOOP},
  {(cstring)"Acut",    LEVEL_LOOP},
  {(cstring)"Bcut",    LEVEL_LOOP},
  {(cstring)"Ccut",    LEVEL_LOOP},
  {(cstring)"Dcut",    LEVEL_LOOP},
  {(cstring)"Ecut",    LEVEL_LOOP},
  {(cstring)"Fcut",    LEVEL_LOOP},
  {(cstring)"A0",      LEVEL_LOOP},
  {(cstring)"A00",     LEVEL_LOOP},
  {(cstring)"B0",      LEVEL_LOOP},
  {(cstring)"B1",      LEVEL_LOOP},
  {(cstring)"B00",     LEVEL_LOOP},
  {(cstring)"B11",     LEVEL_LOOP},
  {(cstring)"B001",    LEVEL_LOOP},
  {(cstring)"B111",    LEVEL_LOOP},
  {(cstring)"DB0",     LEVEL_LOOP},
  {(cstring)"DB1",     LEVEL_LOOP},
  {(cstring)"DB00",    LEVEL_LOOP},
  {(cstring)"C0",      LEVEL_LOOP},
  {(cstring)"D0",      LEVEL_LOOP},
  {(cstring)"E0",      LEVEL_LOOP},
  {(cstring)"F0",      LEVEL_LOOP}
};

static byte nothing[] = "";
static int debug = 0;
static FILE *stddeb;
static const char *levname[] = {"LOOP", "COEFF", "DEN", "MAT", "SUM"};

typedef struct term {
  struct term *next;
  int flags, nterms[NLEVELS];
  string f[NLEVELS];
  byte expr[];
} TERM;

#define FLAGS_LAST 1
#define FLAGS_COLL 2

/******************************************************************/

#define Atoi(s) atoi((const char *)(s))
#define Strlen(s) strlen((const char *)(s))
#define Strchr(s, c) (string)strchr((const char *)(s), (char)(c))
#define Strspn(s, sx) strspn((const char *)(s), (const char *)(sx))
#define Strstr(s1, s2) (string)strstr((const char *)(s1), (const char *)(s2))
#define Strcmp(s1, s2) strcmp((const char *)(s1), (const char *)(s2))
#define Strncmp(s1, s2, n) strncmp((const char *)(s1), (const char *)(s2), n)
#define Strncpy(s1, s2, n) strncpy((char *)(s1), (const char *)(s2), n)
#define MLPutStr(mlp, s) MLPutByteString(mlp, (string)s, Strlen(s))

/******************************************************************/

static void PrintChain(TERM *tp, const char *info)
{
  int n = 0, lev;
  while( tp ) {
    fprintf(stddeb, "\n" BLUE "%s term %d (%p):" RESET "\n", info, ++n, tp);
    for( lev = 0; lev < NLEVELS; ++lev )
      if( *tp->f[lev] )
        fprintf(stddeb, BLUE "  %s" RESET " |%s|\n", levname[lev], tp->f[lev]);
    fprintf(stddeb, BLUE "  flags: " RESET "subM=%d coll=%d last=%d\n",
      tp->flags >> 2, (tp->flags & FLAGS_COLL) >> 1, tp->flags & FLAGS_LAST);
    tp = tp->next;
  }
  fprintf(stddeb, "\n");
}

/******************************************************************/

static inline string MLString(MLINK mlp)
{
  cstring s;
  string d;
  int n;

  if( MLGetByteString(mlp, &s, &n, ' ') == 0 ) {
    MLClearError(mlp);
    MLNewPacket(mlp);
    return NULL;
  }
  d = malloc(n + 1);
  if( d ) {
    memcpy(d, s, n);
    d[n] = 0;
  }
  MLReleaseByteString(mlp, s, n);
  return d;
}

/******************************************************************/

static inline void MLSendPacket(MLINK mlp)
{
  MLEndPacket(mlp);
  while( MLNextPacket(mlp) != RETURNPKT )
    MLNewPacket(mlp);
}

/******************************************************************/

static inline void MLEmitMessage(MLINK mlp, const char *tag, cstring arg)
{
  MLPutFunction(mlp, "EvaluatePacket", 1);

  MLPutFunction(mlp, "Message", (arg) ? 2 : 1);
  MLPutFunction(mlp, "MessageName", 2);
  MLPutSymbol(mlp, "ReadForm");
  MLPutStr(mlp, tag);
  if( arg ) MLPutStr(mlp, arg);
  MLSendPacket(mlp);
  MLNewPacket(mlp);	/* discard returned Null */
}

/******************************************************************/

static inline int MLPutExpr(MLINK mlp, TERM *tp, cint lev)
{
  if( debug & DEB_OUT )
    fprintf(stddeb, RED "  %s" RESET " |%s|\n", levname[lev], tp->f[lev]);
  MLPutFunction(mlp, "ToExpression", 1);
  return MLPutStr(mlp, tp->f[lev]);
}

/******************************************************************/

static inline void InsertLor(string s)
{
  s -= 3;
  do s[3] = s[0]; while( *--s < 'A' );
  memcpy(s, "Lor[", 4);
}

/******************************************************************/

static inline void Shift(TERM *tp, csize_t len, cstring old, cstring new)
{
  if( new != old ) {
    string *f;
    for( f = tp->f; f < &tp->f[NLEVELS]; ++f )
      if( (size_t)(*f - old) < len ) *f += new - old;
  }
}

/******************************************************************/

static inline TERM *Resize(TERM *tp, csize_t len, cint newsize)
{
  TERM *new;
  assert(new = realloc(tp, newsize));
  Shift(new, len, tp->expr, new->expr);
  return new;
}

/******************************************************************/

static inline void MoveToEnd(TERM *tp, cint lev, cstring end)
{
  string s = tp->f[lev];
  csize_t len = Strlen(s);
  cstring begin = s + len + 1;
  csize_t rest = end - begin;
  string tmp;

  assert(tmp = malloc(len));
  memcpy(tmp, s, len);
  memmove(s, begin, rest);
  s += rest;
  *s++ = 0;
  memcpy(s, tmp, len);
  free(tmp);

  Shift(tp, rest, begin, tp->f[lev]);
  tp->f[lev] = s;
}

/******************************************************************/

static string PutFactor(string to, cstring from)
{
  if( *from ) {
    csize_t len = Strlen(from) + 1;
    memcpy(to, from, len);
    return to + len;
  }
  *to++ = '1';
  *to++ = 0;
  return to;
}

/******************************************************************/

static TERM *CollectLoop(TERM *termp, csize_t maxloopsize)
{
  TERM *old;

  do {
    TERM *tp;
    string s = termp->f[LEVEL_LOOP];
    old = termp;

    while( (tp = old->next) ) {
      int lev;
      for( lev = LEVEL_LOOP + 1; lev < NLEVELS; ++lev )
        if( Strcmp(termp->f[lev], tp->f[lev]) != 0 ) {
          old = tp;
          goto loop;
        }
      if( (termp->flags & FLAGS_COLL) == 0 ) {
        assert(termp->f[LEVEL_LOOP] = malloc(maxloopsize));
        s = PutFactor(termp->f[LEVEL_LOOP], s);
        termp->flags |= FLAGS_COLL;
      }
      s[-1] = '+';
      s = PutFactor(s, tp->f[LEVEL_LOOP]);
      old->next = tp->next;
      free(tp);
loop: ;
    }

    if( termp->flags & FLAGS_COLL )
      assert(termp->f[LEVEL_LOOP] =
        realloc(termp->f[LEVEL_LOOP], s - termp->f[LEVEL_LOOP]));

    old = termp;
  } while( (termp = termp->next) );

  return old;
}

/******************************************************************/

static TERM *OrderChain(TERM *t1p, cint level)
{
  TERM *old1;
  int *const nterms = &t1p->nterms[level];
  int c = 0;

  do {
    TERM *t2p, *old2, *ini = t1p;
    int c2 = 0;
    ++c;

    do {
      ++c2;
      old1 = t1p;
      t1p = old1->next;
      if( t1p == NULL ) goto next;
    } while( Strcmp(ini->f[level], t1p->f[level]) == 0 );

    t2p = t1p;

    do {
      old2 = t2p;
over:
      t2p = old2->next;
      if( t2p == NULL ) goto next;
    } while( Strcmp(ini->f[level], t2p->f[level]) != 0 );

    old1->next = t2p;
    old1 = t2p;
    old2->next = t2p->next;
    ++c2;
    goto over;

next:
    if( level > LEVEL_COEFF ) {
      old1->next = NULL;
      old1 = OrderChain(ini, level - 1);
    }
    else ini->nterms[level - 1] = c2;
  } while( (old1->next = t1p) );

  *nterms = c;
  return old1;
}

/******************************************************************/

static TERM *Transmit(TERM *tp, int level, int simp)
{
  cint nterms = tp->nterms[level];
  int term;

  if( level == LEVEL_SUMOVER ) MLPutFunction(stdlink, "List", nterms);
  else if( nterms > 1 ) MLPutFunction(stdlink, "Plus", nterms);

  for( term = 1; term <= nterms; ++term ) {
    int lev, ntimes = (*tp->f[level] != 0);

    if( debug & DEB_OUT )
      fprintf(stddeb, RED "term %d/%d of %s" RESET "\n",
        term, nterms, levname[level]);

    for( lev = level - 1; lev > LEVEL_LOOP; --lev ) {
      ++ntimes;
      if( tp->nterms[lev] > 1 ) goto sendit;
      if( *tp->f[lev] == 0 ) --ntimes;
    }
	/* OrderChain goes down only to LEVEL_COEFF, hence: */
    if( *tp->f[LEVEL_LOOP] ) ++ntimes;

sendit:
    switch( ntimes ) {
    case 0:
      MLPutInteger(stdlink, 1);
      break;

    default:
      MLPutFunction(stdlink, "Times", ntimes);
    case 1:
      if( *tp->f[level] ) MLPutExpr(stdlink, tp, level);
      for( lev = level - 1; lev > LEVEL_LOOP; --lev ) {
        if( simp && lev == LEVEL_MAT - 1 )
          MLPutFunction(stdlink, "FormMat", 1);
        if( tp->nterms[lev] > 1 ) {
          tp = Transmit(tp, lev, simp);
          goto loop;
        }
        if( *tp->f[lev] ) MLPutExpr(stdlink, tp, lev);
      }
      if( *tp->f[LEVEL_LOOP] ) MLPutExpr(stdlink, tp, LEVEL_LOOP);
      break;
    }
    tp = tp->next;
loop: ;
  }

  return tp;
}

/******************************************************************/

static void ReadForm(FILE *file, cstring errarg)
{
  TERM *termp = NULL, **exprp = NULL, *ep;
  TERM *anchor = NULL, **last = &anchor;
  byte br[64], lhs[64];
  int b = 0, thislev = 0;
  int nterms = 0, nexpr = 0, nsubM = 0;
  size_t maxloopsize = 0;
  cstring beg = NULL;
  string delim = NULL, di = NULL;
  enum { nfun = sizeof funtab/sizeof *funtab };
  int lineno = 0;
  const char *errtag = "noopen";
  byte line[256];
  byte errmsg[512];

  *lhs = 0;

  if( file ) for( ; ; ) {
    string errend = errmsg;
    cstring si;
    size_t termsize;
    int i;

    *errend = 0;

nextline:
    for( ; ; ) {
      string pos;

      *line = 0;
      si = (string)fgets((char *)line, sizeof line, file);

      if( si == NULL || MLAbort ) {
        fclose(file);

        if( MLAbort ) goto abort;

        if( errend > errmsg ) {
          errend[-1] = 0;	/* discard last \n */
          errtag = "formerror";
          errarg = errmsg;
          goto error;
        }

        if( nexpr == 0 ) goto error;

        /* successful exit */
        if( debug & DEB_OUT )
          fprintf(stddeb, "\n" RED "nexpr = %d  nsubM = %d" RESET "\n",
            nexpr, nsubM);
        MLPutFunction(stdlink, "FormExpr", nexpr - nsubM);
        if( nsubM ) MLPutFunction(stdlink, "CompoundExpression", nsubM + 1);
        for( termp = anchor; nexpr--; ) {
          cint n = termp->flags >> 2;
          if( n ) {
            if( debug & DEB_OUT )
              fprintf(stddeb, "\n" RED "%s[%d] = " RESET, lhs, n);
            MLPutFunction(stdlink, "Set", 2);
            MLPutFunction(stdlink, (const char *)lhs, 1);
            MLPutInteger(stdlink, n);
            if( Strncmp(termp->f[LEVEL_COEFF], "mulM", 4) == 0 )
              MLPutFunction(stdlink, "FormSub", 1);
            termp = Transmit(termp, LEVEL_MAT, 0);
          }
          else {
            if( debug & DEB_OUT )
              fprintf(stddeb, "\n" RED "expr = " RESET);
            termp = Transmit(termp, LEVEL_SUMOVER, 1);
          }
        }

        goto quit;
      }

      errtag = "nooutput";
      if( debug & DEB_IN ) fprintf(stddeb, "%06d %c%s",
        ++lineno,
        (exprp) ? (b ? '&' : ' ') : '*',
        (const char *)line);

      if( (pos = Strstr(si, "-->")) ||
          (pos = Strstr(si, "==>")) ||
          (pos = Strstr(si, "===")) ) {
        pos += 4;
        if( Strstr(errmsg, pos) == NULL ) {
          Strncpy(errend, pos, errmsg + sizeof errmsg - errend);
          errend += Strlen(errend);
        }
        continue;
      }

#define FinalizeTerm() \
  *di++ = 0; \
  if( termp->f[LEVEL_LOOP] ) \
    maxloopsize += Strlen(termp->f[LEVEL_LOOP]) + 2; \
  *last = termp = Resize(termp, di - (string)termp, di - (string)termp); \
  last = &termp->next

      if( exprp ) {
        if( *si >= ' ' ) break;  /* catch both \n and \r (on Windows) */
        if( di == termp->expr ) continue;
        FinalizeTerm();
      }
      else if( (pos = Strchr(si, '=')) == NULL ) continue;

      assert(termp = malloc(termsize = sizeof *termp + TERMBUF));
      termp->next = NULL;
      termp->flags = 0;
      memset(termp->nterms, 0, sizeof termp->nterms);
      for( i = 0; i < NLEVELS; ++i ) termp->f[i] = nothing;
      beg = delim = di = termp->f[thislev = LEVEL_COEFF] = termp->expr;
      ++nterms;

      if( exprp == NULL ) {
        string s = memchr(si, '(', pos - si);
        if( s ) {
          termp->flags = Atoi(s + 1) << 2;
          if( *lhs == 0 ) {
            si += Strspn(si, " ");
            memcpy(lhs, si, s - si);
            lhs[s-si] = 0;
          }
          ++nsubM;
        }
        si = pos + 1;
        exprp = last;
        ++nexpr;
        break;
      }
    }

    if( di > (string)termp + termsize - sizeof line ) {
      csize_t len = di - (string)termp;
      termp = Resize(termp, len, termsize += TERMBUF);
      di = (string)termp + len;
    }

    while( *si ) {
      byte c = *si++;
      if( c <= ' ' ) continue;

      switch( c ) {
      case '+':
      case '-':
      case '*':
        if( b == 0 ) delim = di + 1;
        break;

      case '(':
        if( b == 0 ) {
          int newlev = LEVEL_COEFF;
          *di = 0;
          if( *delim ) for( i = 0; i < nfun; ++i )
            if( Strcmp(delim, funtab[i].name) == 0 ) {
              newlev = funtab[i].level;
              break;
            }

          if( thislev != newlev ) {
            if( *termp->f[newlev] )
              MoveToEnd(termp, newlev, delim - 1);
            else {
              if( delim == beg ) termp->f[thislev] = nothing;
              else delim[-1] = 0;
              termp->f[newlev] = delim;
            }
            thislev = newlev;
          }
        }

        if( di == beg || Strchr("+-*/^,([", di[-1]) ) br[b++] = ')';
        else c = '[', br[b++] = ']';
        break;

      case ')':
        if( b > 0 ) c = br[--b];
        break;

      case '[':
        *di++ = '\\';
        break;

      case ';':
        FinalizeTerm();
        ep = *exprp;
        if( maxloopsize ) termp = CollectLoop(ep, maxloopsize);
        if( debug & DEB_CHU ) PrintChain(ep, "unordered");
        if( (ep->flags >> 2) ) ep->nterms[NLEVELS-2] = nterms;
        else {
          termp = OrderChain(ep, NLEVELS - 1);
          if( debug & DEB_CHO ) PrintChain(ep, "ordered");
        }
        termp->flags |= FLAGS_LAST;
        maxloopsize = nterms = 0;
        exprp = NULL;
        goto nextline;

      case '_':
        *di++ = c = '$';
        break;

      case '?':
        InsertLor(di++);
        c = ']';
        break;
      }

      *di++ = c;
    }
  }

error:
  MLEmitMessage(stdlink, errtag, errarg);

abort:
  MLPutFunction(stdlink, "Abort", 0);

quit:
  MLEndPacket(stdlink);

  while( anchor ) {
    TERM *tp = anchor;
    anchor = anchor->next;
    if( tp->flags & FLAGS_COLL ) free(tp->f[LEVEL_LOOP]);
    free(tp);
  }
}

/******************************************************************/

static void readform_file(void)
{
  string filename = MLString(stdlink);
  ReadForm(fopen((const char *)filename, "r"), filename);
  free(filename);
}

/******************************************************************/

static void readform_exec(void)
{
  int fd[2], status;
  pid_t pid;
  string argv[16], *argp = argv;
  string p = MLString(stdlink);

  while( (p = Strchr(*argp++ = p, '|')) ) *p++ = 0;
  argp[0] = MLString(stdlink);
  argp[1] = NULL;

  signal(SIGCHLD, SIG_IGN);
  assert(pipe(fd) != -1 && (pid = fork()) != -1);

  if( pid == 0 ) {
    usleep(500);
    close(fd[0]);
    dup2(fd[1], 1);
    dup2(fd[1], 2);
    close(fd[1]);
    exit(execvp((char *)argv[0], (char **)argv));
  }

  close(fd[1]);
  ReadForm(fdopen(fd[0], "r"), argv[0]);

  kill(pid, SIGKILL);
  wait(&status);

  free(argp[0]);
  free(argv[0]);
}

/******************************************************************/

static void readformdebug(cint deb)
{
  cstring filename = MLString(stdlink);
  debug = deb;

  stddeb = stderr;
  if( filename && *filename ) {
    stddeb = fopen((const char *)filename, "w");
    if( stddeb == NULL ) {
      MLEmitMessage(stdlink, "noopen", filename);
      MLPutSymbol(stdlink, "$Failed");
      MLEndPacket(stdlink);
      debug = 0;
      return;
    }
    setbuf(stddeb, NULL);
  }

  MLPutSymbol(stdlink, "True");
  MLEndPacket(stdlink);
}

/******************************************************************/

int main(int argc, char **argv)
{
  int fd;

	/* make sure a pipe will not overlap with 0, 1, 2 */
  do fd = open("/dev/null", O_WRONLY); while( fd <= 2 );
  close(fd);

  return MLMain(argc, argv);
}

