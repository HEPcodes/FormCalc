:Begin:
:Function: readform_file
:Pattern: ReadForm[filename_String, debug_Integer]
:Arguments: {filename, debug}
:ArgumentTypes: {String, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: readform_exec
:Pattern: ReadForm[formcmd_String, filename_String, debug_Integer]
:Arguments: {formcmd, filename, debug}
:ArgumentTypes: {String, String, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: clearabbr
:Pattern: ClearAbbr[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Evaluate: ReadForm::noopen = "Cannot open `1`."

:Evaluate: ReadForm::nooutput =
  "Something went wrong, there was no output from FORM."

:Evaluate: ReadForm::toomany =
  "Too many expressions.  Increase MAXEXPR in ReadForm.tm."

:Evaluate: ReadForm::formerror = "`1`"


/*
	ReadForm.tm
		reads FORM output back into Mathematica
		this file is part of FormCalc
		last modified 9 Jun 09 th

Note: FORM code must have
	1. #- (no listing),
	2. off stats,
	3. should produce output with print (not print +s).
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include <sys/wait.h>

#define MAXEXPR 5000
#define TERMBUF 500000
#define STRINGSIZE 32767

#define Allocate(p, n) \
  if( (p = malloc(n)) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }

#define Reallocate(p, n) \
  if( (p = realloc(p, n)) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }

#define DEBUG "\e[31m"
#define RESET "\e[0m\n"


/*

A term in the FORM output is organized into the TERM structure
in the following way:

 ____4_____     __3___     __2___     ___0____     _____1_____
/          \   /      \   /      \   /        \   /           \
SumOver(...) * Mat(...) * Den(...) * paveM(...) * ..... * (...)

Hierarchy of collecting:
4. SumOver
3. Mat
2. Den
1. [coefficient]
0. paveM

*/

#define LEVEL_PAVE 0
#define LEVEL_COEFF 1
#define LEVEL_DEN 2
#define LEVEL_MAT 3
#define LEVEL_SUMOVER 4
#define LEVELS 5

typedef const int cint;
typedef char *string;
typedef MLCONST char *cstring;

typedef struct {
  cstring name;
  int level;
} FUN;

static FUN funtab[] = {
  {"SumOver", LEVEL_SUMOVER},
  {"Mat",     LEVEL_MAT},
  {"Den",     LEVEL_DEN},
  {"paveM",   LEVEL_PAVE}
};

typedef struct term {
  struct term *last;
  string f[LEVELS];
  int nterms[LEVELS], coll;
} TERM;

typedef struct btree {
  struct btree *lt, *gt;
  string abbr;
  char expr[0];
} BTREE;

static char zero[] = "";
static BTREE *root = NULL;

/******************************************************************/

static inline void MLSendPacket(MLINK mlp)
{
  MLEndPacket(mlp);
  while( MLNextPacket(mlp) != RETURNPKT )
    MLNewPacket(mlp);
}

/******************************************************************/

static inline void MLEmitMessage(MLINK mlp, cstring tag, cstring arg)
{
  MLPutFunction(mlp, "EvaluatePacket", 1);

  MLPutFunction(mlp, "Message", (arg) ? 2 : 1);
  MLPutFunction(mlp, "MessageName", 2);
  MLPutSymbol(mlp, "ReadForm");
  MLPutString(mlp, tag);
  if( arg ) MLPutString(mlp, arg);
  MLSendPacket(mlp);
  MLNewPacket(mlp);	/* discard returned Null */
}

/******************************************************************/

static inline void InsertLor(string s)
{
  s -= 3;
  do s[3] = s[0]; while( *--s < 'A' );
  memcpy(s, "Lor[", 4);
}

/******************************************************************/

static string GetAbbr(string expr)
{
  BTREE *lp, **node = &root;
  cstring abbr;
  int exprlen, abbrlen;

  while( (lp = *node) ) {
    cint t = strcmp(expr, lp->expr);
    if( t == 0 ) return lp->abbr;
    node = (t < 0) ? &lp->lt : &lp->gt;
  }

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "FormEval", 1);
  MLPutString(stdlink, expr);
  MLSendPacket(stdlink);

  if( MLGetString(stdlink, &abbr) == 0 ) return expr;

  exprlen = strlen(expr);
  abbrlen = strlen(abbr);
  Allocate(lp, exprlen + abbrlen + 2 + sizeof(BTREE));
  lp->lt = lp->gt = NULL;
  *node = lp;
  memcpy(lp->expr, expr, exprlen);
  lp->abbr = lp->expr + exprlen + 1;
  memcpy(lp->abbr, abbr, abbrlen);
  lp->abbr[abbrlen] = 0;

  MLDisownString(stdlink, abbr);

  return lp->abbr;
}

/******************************************************************/

static TERM *Downsize(TERM *tp, cstring end)
{
  cstring begin = (cstring)tp;
  TERM *new = realloc(tp, end - begin);
  cint offset = begin - (cstring)new;

  if( offset ) {
    string *f;
    for( f = new->f; f < &new->f[LEVELS]; ++f )
      if( *f >= begin && *f <= end ) *f -= offset;
  }

  return new;
}

/******************************************************************/

static string PutFactor(string to, cstring from)
{
  if( *from ) return memccpy(to, from, 0, STRINGSIZE);
  *to++ = '1';
  *to++ = 0;
  return to;
}

/******************************************************************/

static void CollectPaVe(TERM *termp, cint maxpavesize)
{
  string s;

  do {
    TERM *old = termp, *tp;

    while( (tp = old->last) ) {
      int i;
      for( i = LEVEL_PAVE + 1; i < LEVELS; ++i )
        if( strcmp(termp->f[i], tp->f[i]) != 0 ) {
          old = tp;
          goto loop;
        }
      if( termp->coll == 0 ) {
        s = termp->f[LEVEL_PAVE];
        Allocate(termp->f[LEVEL_PAVE], maxpavesize);
        s = PutFactor(termp->f[LEVEL_PAVE], s);
        termp->coll = 1;
      }
      s[-1] = '+';
      s = PutFactor(s, tp->f[LEVEL_PAVE]);
      old->last = tp->last;
      free(tp);
loop: ;
    }

    if( termp->coll ) termp->f[LEVEL_PAVE] =
      realloc(termp->f[LEVEL_PAVE], s - termp->f[LEVEL_PAVE]);
  } while( (termp = termp->last) );
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
      t1p = old1->last;
      if( t1p == NULL ) goto next;
    } while( strcmp(ini->f[level], t1p->f[level]) == 0 );

    t2p = t1p;

    do {
      old2 = t2p;
over:
      t2p = old2->last;
      if( t2p == NULL ) goto next;
    } while( strcmp(ini->f[level], t2p->f[level]) != 0 );

    old1->last = t2p;
    old1 = t2p;
    old2->last = t2p->last;
    ++c2;
    goto over;

next:
    if( level > LEVEL_COEFF ) {
      old1->last = NULL;
      old1 = OrderChain(ini, level - 1);
    }
    else ini->nterms[level - 1] = c2;
  } while( (old1->last = t1p) );

  *nterms = c;
  return old1;
}

/******************************************************************/

static TERM *Transmit(TERM *tp, int level)
{
  int n = tp->nterms[level];

  if( level == LEVEL_SUMOVER ) MLPutFunction(stdlink, "List", n);
  else if( n > 1 ) MLPutFunction(stdlink, "Plus", n);

  while( n-- ) {
    int i, ntimes = (*tp->f[level] != 0);

    for( i = level - 1; i > LEVEL_PAVE; --i ) {
      ++ntimes;
      if( tp->nterms[i] > 1 ) goto sendit;
      if( *tp->f[i] == 0 ) --ntimes;
    }
	/* OrderChain goes down only to LEVEL_COEFF, hence: */
    if( *tp->f[LEVEL_PAVE] ) ++ntimes;

sendit:
    switch( ntimes ) {
    case 0:
      MLPutInteger(stdlink, 1);
      break;

    default:
      MLPutFunction(stdlink, "Times", ntimes);
    case 1:
      if( *tp->f[level] ) {
        MLPutFunction(stdlink, "ToExpression", 1);
        MLPutString(stdlink, tp->f[level]);
      }
      for( i = level - 1; i > LEVEL_PAVE; --i ) {
        if( tp->nterms[i] > 1 ) {
          tp = Transmit(tp, i);
          goto loop;
        }
        if( *tp->f[i] ) {
          MLPutFunction(stdlink, "ToExpression", 1);
          MLPutString(stdlink, tp->f[i]);
        }
      }
      if( *tp->f[LEVEL_PAVE] ) {
        MLPutFunction(stdlink, "ToExpression", 1);
        MLPutString(stdlink, tp->f[LEVEL_PAVE]);
      }
      break;
    }
    tp = tp->last;
loop: ;
  }

  return tp;
}

/******************************************************************/

static void ReadForm(FILE *file, cint debug)
{
  TERM *expressions[MAXEXPR], **exprp = expressions;
  TERM *termp = NULL, *tp;
  char br[64];
  int inexpr = 0, maxpavesize = 0, b = 0, thislev;
  enum { nfun = sizeof(funtab)/sizeof(FUN) };

  maxpavesize = 0;

  for( ; ; ) {
    char line[1024];
    string di, delim;
    cstring si, beg;
    char errmsg[512];
    string errend = errmsg;
    int i;

    *errend = 0;

nextline:
    for( ; ; ) {
      string pos;

      *line = 0;
      si = fgets(line, sizeof(line), file);
      if( debug ) fputs(line, stderr);

      if( si == NULL || MLAbort ) {
        int nexpr;
        TERM **ep;

        if( !feof(file) ) goto abort;

        if( errend > errmsg ) {
          errend[-1] = 0;	/* discard last \n */
          MLEmitMessage(stdlink, "formerror", errmsg);
          goto abort;
        }

        nexpr = (int)(exprp - expressions);
        if( nexpr == 0 ) {
          MLEmitMessage(stdlink, "nooutput", NULL);
          goto abort;
        }

        /* successful exit */
        MLPutFunction(stdlink, "FormExpr", nexpr);
        for( ep = expressions; ep < exprp; ++ep ) {
          OrderChain(*ep, LEVELS - 1);
          Transmit(*ep, LEVELS - 1);
        }
        goto quit;
      }

      if( (pos = strstr(si, "-->")) ||
          (pos = strstr(si, "==>")) ||
          (pos = strstr(si, "===")) ) {
        pos += 4;
        if( strstr(errmsg, pos) == NULL ) {
          strncpy(errend, pos, errmsg + sizeof(errmsg) - errend);
          errend += strlen(errend);
        }
        continue;
      }

      if( inexpr ) {
        if( *si != '\n' ) break;
        if( di == (cstring)tp + sizeof(TERM) ) continue;
        *di++ = 0;
        termp = Downsize(termp, di);
      }
      else if( (pos = strchr(si, '=')) == NULL ) continue;

      Allocate(tp, sizeof(TERM) + TERMBUF);
      beg = delim = di = (string)tp + sizeof(TERM);
      tp->last = termp;
      termp = tp;
      for( i = 0; i < LEVELS; ++i ) tp->f[i] = zero;
      tp->f[thislev = LEVEL_COEFF] = di;
      tp->coll = 0;

      if( !inexpr ) {
        inexpr = 1;
        si = pos + 1;
        break;
      }
    }

    while( *si ) {
      char c = *si++;
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
          for( i = 0; i < nfun; ++i )
            if( strcmp(delim, funtab[i].name) == 0 ) {
              newlev = funtab[i].level;
              break;
            }

          if( thislev != newlev ) {
            if( delim > beg ) delim[-1] = 0;
            switch( thislev ) {
            case LEVEL_MAT:
              termp->f[LEVEL_MAT] = GetAbbr(termp->f[LEVEL_MAT]);
              break;
            case LEVEL_PAVE:
              maxpavesize += (int)(delim - termp->f[LEVEL_PAVE]) + 2;
              break;
            }
            thislev = newlev;
            termp->f[thislev] = delim;
          }
        }

        if( di == beg || strchr("+-*/^,([", di[-1]) ) br[b++] = ')';
        else c = '[', br[b++] = ']';
        break;

      case ')':
        if( b > 0 ) c = br[--b];
        break;

      case '[':
        *di++ = '\\';
        break;

      case ';':
        *di++ = 0;
        *exprp++ = termp = Downsize(termp, di);
        if( exprp >= expressions + MAXEXPR ) {
          MLEmitMessage(stdlink, "toomany", NULL);
          goto abort;
        }
        if( maxpavesize ) CollectPaVe(termp, maxpavesize);
        termp = NULL;
        maxpavesize = inexpr = 0;
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

abort:
  MLPutFunction(stdlink, "Abort", 0);

quit:
  MLEndPacket(stdlink);

  while( exprp > expressions ) {
    TERM *last;
    for( tp = *--exprp; tp; tp = last ) {
      if( tp->coll ) free(tp->f[LEVEL_PAVE]);
      last = tp->last;
      free(tp);
    }
  }
}

/******************************************************************/

static void readform_file(cstring filename, cint debug)
{
  FILE *file = fopen(filename, "r");

  if( file ) {
    ReadForm(file, debug);
    fclose(file);
  }
  else {
    MLEmitMessage(stdlink, "noopen", filename);
    MLPutFunction(stdlink, "Abort", 0);
    MLEndPacket(stdlink);
  }
}

/******************************************************************/

static inline int writeall(cint h, cstring buf, long n)
{
  long w = 0;
  while( (w = write(h, buf, n -= w)) > 0 );
  if( w < 0 ) close(h);
  return w;
}

/******************************************************************/

static int ToMma(cint hw, string expr)
{
  int b = -10000;
  cstring result, r;
  string s;
  char c;

//fprintf(stderr, DEBUG "to mma (%d bytes)" RESET, strlen(expr));
//fprintf(stderr, "expr=|%s|\n", expr);
  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "FormEvalDecl", 1);
  MLPutString(stdlink, expr);
  MLSendPacket(stdlink);

  if( MLGetString(stdlink, &result) == 0 ) return 0;

//fprintf(stderr, DEBUG "from mma raw (%d bytes)" RESET, strlen(result));
//fprintf(stderr, "expr=|%s|\n", result);
  for( r = result, s = expr; (c = *r++); ) {
    if( c <= ' ' && b >= 0 ) continue;
    switch( c ) {
    case '$':
      if( *r == c ) ++r, c = '_';
      break;
    case '[':
    case '(':
      c = '(';
      ++b;
      break;
    case ']':
    case ')':
      c = ')';
      --b;
      break;
    case '\\':
      if( *r == '0' ) c = strtol(r, (string *)&r, 8);
      break;
    case ',':
      if( b ) break;
      *s++ = '\n';
    case '{':
//fprintf(stderr, DEBUG "from mma (%d bytes)" RESET, (int)(s - expr));
//*s=0; fprintf(stderr, "expr=|%s|\n", expr);
      *s++ = '\n';
      if( writeall(hw, expr, s - expr) < 0 ) {
        MLDisownString(stdlink, result);
        return 0;
      }
      s = expr;
      b = 0;
      continue;
    }
    *s++ = c;
  }

  sync();  /* another sync needed for Mac OS, again unclear why */

  MLDisownString(stdlink, result);
  return 1;
}

/******************************************************************/

static void *ExtIO(void *h)
{
  cint hw = ((int *)h)[1];
  cint hr = ((int *)h)[2];
  string expr;
  char br[64];
  enum { blocksize = 40960, linesize = 512, ahead = 32 };
  int size = 2*blocksize, b, w, n;

  if( (n = read(hr, br, sizeof br) - 1) < 0 ||
      writeall(hw, br, n + sprintf(br + n, ",%d\n", getpid())) < 0 )
    pthread_exit(NULL);

  Allocate(expr, size);

loop:
  w = 1;
  b = 0;

  for( ; ; ) {
    long r = w + ahead;
    string s;

    if( r + linesize > size ) {
      size += blocksize;
      Reallocate(expr, size);
    }

    sync();  /* needed for Mac OS, unclear why */

    n = read(hr, s = expr + r, size - r);
    if( n <= 0 ) {
      free(expr);
      pthread_exit(NULL);
    }

//expr[r+n] = 0;
//fprintf(stderr, "expr=|%s|\n", expr+r);

    do {
      char c = *s++;
      if( c <= ' ' ) continue;
      switch( c ) {
      case '#':
        expr[w++] = '0';
        expr[w++] = '}';
        expr[w] = 0;
        *expr = '{';
        if( ToMma(hw, expr) ) goto loop;
        free(expr);
        pthread_exit(NULL);
      case '?':
        InsertLor(expr + w++);
        c = ']';
        break;
      case '_':
        expr[w++] = c = '$';
        break;
      case '(':
        if( strchr("+-*/^,([{", expr[w-1]) ) br[b++] = ')';
        else c = '[', br[b++] = ']';
        break;
      case ')':
        if( b > 0 ) c = br[--b];
        break;
      }
      expr[w++] = c;
    } while( --n );
  }
}

/******************************************************************/

static void readform_exec(cstring cmd, cstring filename, cint debug)
{
  int h[6] = {-1, -1, -1, -1, -1, -1}, i;
  pthread_t tid;
  void *tj = NULL;
  pid_t pid = -1;
  FILE *file;
  char arg[32];
  cstring argv[5] = {cmd, filename, NULL, filename, NULL};

  if( pipe(&h[0]) != -1 &&
      pipe(&h[2]) != -1 &&
      pthread_create(&tid, NULL, ExtIO, h) == 0 ) {
    tj = (void *)1;
    sprintf(arg, "%d,%d", h[0], h[3]);
    argv[1] = "-pipe";
    argv[2] = arg;
  }

  signal(SIGCHLD, SIG_IGN);
  if( pipe(&h[4]) == -1 || (pid = fork()) == -1 ) goto abort;

  if( pid == 0 ) {
    if( h[1] != -1 ) close(h[1]);
    if( h[2] != -1 ) close(h[2]);
    close(h[4]);
    dup2(h[5], 1);
    close(h[5]);
    exit(execvp(cmd, (char *const *)argv));
  }

  if( h[0] != -1 ) close(h[0]);
  if( h[3] != -1 ) close(h[3]);
  close(h[5]);
  h[0] = h[3] = h[5] = -1;

  file = fdopen(h[4], "r");
  if( file ) {
    ReadForm(file, debug);
    fclose(file);
    h[4] = -1;
  }
  else {
abort:
    MLEmitMessage(stdlink, "noopen", cmd);
    MLPutFunction(stdlink, "Abort", 0);
    MLEndPacket(stdlink);
  }

  if( pid > 0 ) {
    kill(pid, SIGKILL);
    wait(&i);
  }

  for( i = 5; i >= 0; --i ) if( h[i] != -1 ) close(h[i]);

  if( tj ) pthread_join(tid, &tj);
}

/******************************************************************/

static void CutBranch(BTREE *node)
{
  if( node ) {
    CutBranch(node->lt);
    CutBranch(node->gt);
    free(node);
  }
}

/******************************************************************/

static void clearabbr(void)
{
  CutBranch(root);
  root = NULL;
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}

/******************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

