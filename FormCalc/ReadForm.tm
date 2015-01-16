:Begin:
:Function: readform
:Pattern: ReadForm[filename_String, debug_Integer]
:Arguments: {filename, debug}
:ArgumentTypes: {String, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: clearcache
:Pattern: ClearCache[]
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
		last modified 15 Jan 07 th

Note: FORM code must have
	1. #- (no listing),
	2. off stats,
	3. should produce output with print (not print +s).
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXEXPR 5000
#define TERMBUF 500000
#define STRINGSIZE 32767

#define Allocate(p, n) \
  if( (p = malloc(n)) == NULL ) { \
    fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
    exit(1); \
  }


/*

A term in the FORM output is organized into the TERM structure
in the following way:

 ____4_____     __3___     __2___     ___0___     _____1_____
/          \   /      \   /      \   /       \   /           \
SumOver(...) * Mat(...) * Den(...) * pave(...) * ..... * (...)

Hierarchy of collecting:
4. SumOver
3. Mat
2. Den
1. [coefficient]
0. pave

*/

#define LEVEL_PAVE 0
#define LEVEL_COEFF 1
#define LEVEL_DEN 2
#define LEVEL_MAT 3
#define LEVEL_SUMOVER 4
#define LEVELS 5

typedef char *string;
typedef MLCONST char *cstring;
typedef const int cint;

typedef struct {
  cstring name;
  int level;
} FUN;

static FUN funtab[] = {
  {"SumOver", LEVEL_SUMOVER},
  {"Mat",     LEVEL_MAT},
  {"Den",     LEVEL_DEN},
  {"pave",    LEVEL_PAVE}
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

#define MLPutString(mlp, s) \
  MLPutByteString(mlp, (unsigned MLCONST char *)s, strlen(s))

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
  int pkt;

  MLPutFunction(mlp, "EvaluatePacket", 1);

  MLPutFunction(mlp, "Message", 2);
  MLPutFunction(mlp, "MessageName", (arg) ? 2 : 1);
  MLPutSymbol(mlp, "ReadForm");
  MLPutString(mlp, tag);
  if( arg ) MLPutString(mlp, arg);
  MLSendPacket(mlp);
  MLNewPacket(mlp);	/* discard returned Null */
}

/******************************************************************/

static string GetAbbr(cstring expr)
{
  BTREE *lp, **node = &root;
  const unsigned char *abbr;
  long exprlen, abbrlen;
  int pkt;

  while( (lp = *node) ) {
    cint t = strcmp(expr, lp->expr);
    if( t == 0 ) return lp->abbr;
    node = (t < 0) ? &lp->lt : &lp->gt;
  }

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "ToString", 1);
  MLPutFunction(stdlink, "ToExpression", 1);
  MLPutString(stdlink, expr);
  MLSendPacket(stdlink);

  MLGetByteString(stdlink, &abbr, &abbrlen, 255);

  exprlen = strlen(expr);
  Allocate(lp, exprlen + abbrlen + 2 + sizeof(BTREE));
  lp->lt = lp->gt = NULL;
  *node = lp;
  strcpy(lp->expr, expr);
  lp->abbr = lp->expr + exprlen + 1;
  memcpy(lp->abbr, abbr, abbrlen);
  lp->abbr[abbrlen] = 0;

  MLDisownByteString(stdlink, abbr, abbrlen);

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
      *(s - 1) = '+';
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

static void OrderChain(TERM *t1p, cint level)
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
      OrderChain(ini, level - 1);
    }
    else ini->nterms[level - 1] = c2;
  } while( (old1->last = t1p) );

  *nterms = c;
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

void readform(cstring filename, cint debug)
{
  TERM *expressions[MAXEXPR], **exprp = expressions;
  TERM *termp = NULL, *tp;
  char brackets[20], *br = brackets;
  int inexpr = 0, maxpavesize = 0, thislev;

  FILE *file = (*filename == '!') ? popen(filename + 1, "r") :
                                    fopen(filename, "r");
  if( file == NULL ) {
    MLEmitMessage(stdlink, "noopen", filename);
    MLPutFunction(stdlink, "Abort", 0);
    MLEndPacket(stdlink);
    return;
  }

  maxpavesize = 0;

  for( ; ; ) {
    char line[1024];
    string di, ind, delim;
    cstring si, beg;
    char errmsg[512];
    string errend = errmsg;
    int i;

    *errend = 0;

nextline:
    for( ; ; ) {
      string pos;

      if( MLAbort ) goto abort;

      if( feof(file) ) {
        int nexpr;
        TERM **ep;

        if( errend > errmsg ) {
          *(errend - 1) = 0;	/* discard last \n */
          MLEmitMessage(stdlink, "formerror", errmsg);
          goto abort;
        }

        nexpr = (int)(exprp - expressions);
        if( nexpr == 0 ) {
          MLEmitMessage(stdlink, "nooutput", NULL);
          goto abort;
        }

        /* successful exit */
        MLPutFunction(stdlink, "List", nexpr);
        for( ep = expressions; ep < exprp; ++ep ) {
          OrderChain(*ep, LEVELS - 1);
          Transmit(*ep, LEVELS - 1);
        }
        goto quit;
      }

      *line = 0;
      si = fgets(line, sizeof(line), file);
      if( debug ) fputs(line, stderr);

      if( (pos = strstr(line, "-->")) ||
          (pos = strstr(line, "==>")) ||
          (pos = strstr(line, "===")) ) {
        pos += 4;
        if( strstr(errmsg, pos) == NULL ) {
          strncpy(errend, pos, errmsg + sizeof(errmsg) - errend);
          errend += strlen(errend);
        }
        continue;
      }

      if( inexpr ) {
        if( *line != '\n' ) break;
        if( di == (cstring)tp + sizeof(TERM) ) continue;
        *di++ = 0;
        termp = Downsize(termp, di);
      }
      else if( (pos = strchr(line, '=')) == NULL ) continue;

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

    for( ; *si; ++si ) {
      if( *si > ' ' ) switch( *si ) {
      case '+':
      case '-':
      case '*':
        *di++ = *si;
        if( br == brackets ) delim = di;
        break;

      case '(':
        if( br == brackets ) {
          int newlev = LEVEL_COEFF;
          *di = 0;
          for( i = 0; i < sizeof(funtab)/sizeof(FUN); ++i )
            if( strcmp(delim, funtab[i].name) == 0 ) {
              newlev = funtab[i].level;
              break;
            }

          if( thislev != newlev ) {
            if( delim > beg ) *(delim - 1) = 0;
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

        if( di == beg || strchr("+-*/^,([", *(di - 1)) )
          *di++ = '(', *br++ = ')';
        else *di++ = '[', *br++ = ']';
        break;

      case ')':
        *di++ = *--br;
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
        *di++ = 'J';
        break;

      case '?':
        ind = di - 2;
        do *(ind + 3) = *ind; while( *--ind != 'N' );
        *ind++ = 'L';
        *ind++ = 'o';
        *ind++ = 'r';
        *ind++ = '[';
        di += 2;
        *di++ = ']';
        break;

      default:
        *di++ = *si;
        break;
      }
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

  ((*filename == '!') ? pclose : fclose)(file);
}

/******************************************************************/

static void cutbranch(BTREE *node)
{
  if( node ) {
    cutbranch(node->lt);
    cutbranch(node->gt);
    free(node);
  }
}

/******************************************************************/

void clearcache(void)
{
  cutbranch(root);
  root = NULL;
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}

/******************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

