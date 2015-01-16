:Begin:
:Function: readform
:Pattern: ReadForm[filename_String]
:Arguments: {filename}
:ArgumentTypes: {String}
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
  "Too many expressions. Increase MAXEXPR in ReadForm.tm."

:Evaluate: ReadForm::formerror = "`1`"


/*
	ReadForm.tm
		reads FORM output back into Mathematica
		this file is part of FormCalc
		last modified 11 Dec 02 th
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXEXPR 5000
#define TERMBUF 100000
#define STRINGSIZE 32767

#ifndef MLCONST
#define MLCONST
#endif


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

typedef struct {
  char *name;
  int level;
} FUN;

FUN funtab[] = {
  {"SumOver", LEVEL_SUMOVER},
  {"Mat",     LEVEL_MAT},
  {"Den",     LEVEL_DEN},
  {"pave",    LEVEL_PAVE}
};

typedef struct term {
  struct term *last;
  char *f[LEVELS];
  int nterms[LEVELS], coll;
} TERM;

typedef struct btree {
  struct btree *lt, *gt;
  char *sym;
  char abb[0];
} BTREE;

char *tok;
TERM *termp, *old1;
int maxintsize;
char zero[] = "";
BTREE *root = NULL;


void report_error(MLCONST char *tag, const char *arg)
{
  int p;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Message", arg ? 2 : 1);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "ReadForm");
  MLPutString(stdlink, tag);
  if( arg ) MLPutString(stdlink, arg);
  MLEndPacket(stdlink);

  do {
    p = MLNextPacket(stdlink);
    MLNewPacket(stdlink);
  } while( p != RETURNPKT );
}


char *getabbr(char *s)
{
  BTREE *lp, **node = &root;
  MLCONST char *mmares;
  int p;

  while( (lp = *node) ) {
    p = strcmp(s, lp->abb);
    if( p == 0 ) return lp->sym;
    node = (p < 0) ? &lp->lt : &lp->gt;
  }

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "ToString", 1);
  MLPutFunction(stdlink, "ToExpression", 1);
  MLPutString(stdlink, s);
  MLEndPacket(stdlink);

  while( MLNextPacket(stdlink) != RETURNPKT )
    MLNewPacket(stdlink);
  MLGetString(stdlink, &mmares);

  p = strlen(s);
  lp = malloc(p + strlen(mmares) + 2 + sizeof(BTREE));
  lp->lt = lp->gt = NULL;
  *node = lp;
  strcpy(lp->abb, s);
  s = lp->sym = lp->abb + p + 1;
  strcpy(s, mmares);

  MLDisownString(stdlink, mmares);
  return s;
}


TERM *downsize(TERM *tp, char *end)
{
  TERM *new;
  char **f;
  int off;

  new = realloc(tp, end - (char *)tp);
  if( (off = (char *)tp - (char *)new) ) {
    for( f = new->f; f < &new->f[LEVELS]; ++f )
      if( *f >= (char *)tp && *f <= end ) *f -= off;
  }
  return new;
}


char *putfac(char *to, char *from)
{
  if( *from ) return memccpy(to, from, 0, STRINGSIZE);
  *to++ = '1';
  *to++ = 0;
  return to;
}


void collect_pave()
{
  TERM *tp, *old;
  char *s;
  int i;

  do {
    for( old = termp; (tp = old->last); ) {
      for( i = LEVEL_PAVE + 1; i < LEVELS; ++i )
        if( strcmp(termp->f[i], tp->f[i]) ) {
          old = tp;
          goto loop;
        }
      if( termp->coll == 0 ) {
        s = termp->f[LEVEL_PAVE];
        s = putfac(termp->f[LEVEL_PAVE] = malloc(maxintsize), s);
        termp->coll = 1;
      }
      *(s - 1) = '+';
      s = putfac(s, tp->f[LEVEL_PAVE]);
      old->last = tp->last;
      free(tp);
loop: ;
    }
    if( termp->coll ) termp->f[LEVEL_PAVE] =
      realloc(termp->f[LEVEL_PAVE], s - termp->f[LEVEL_PAVE]);
  } while( (termp = termp->last) );
}


void orderchain(TERM *t1p, int level)
{
  TERM *t2p, *old2, *ini;
  int c = 0, c2, *nterms = &t1p->nterms[level];

  do {
    ++c;
    c2 = 0;
    ini = t1p;
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
    } while( strcmp(ini->f[level], t2p->f[level]) );
    old1->last = t2p;
    old1 = t2p;
    old2->last = t2p->last;
    ++c2;
    goto over;
next:
    if( level > LEVEL_COEFF ) {
      old1->last = NULL;
      orderchain(ini, level - 1);
    }
    else ini->nterms[level - 1] = c2;
  } while( (old1->last = t1p) );
  *nterms = c;
}


TERM *transmit(TERM *tp, int level)
{
  int n = tp->nterms[level], ntimes, i;

  if( level == LEVEL_SUMOVER ) MLPutFunction(stdlink, "List", n);
  else if( n > 1 ) MLPutFunction(stdlink, "Plus", n);

  while( n-- ) {
    ntimes = *tp->f[level] != 0;
    for( i = level - 1; i > LEVEL_PAVE; --i ) {
      ++ntimes;
      if( tp->nterms[i] > 1 ) goto sendit;
      if( *tp->f[i] == 0 ) --ntimes;
    }
	/* orderchain goes down only to LEVEL_COEFF, hence: */
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
          tp = transmit(tp, i);
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
    }
    tp = tp->last;
loop: ;
  }
  return tp;
}


void readform(const char *filename)
{
  FILE *file;
  char line[1024], *si, *di, *ind, *delim, *beg, **pp;
  char *er, errmsg[512], *erp = errmsg;
  char brackets[20], *br = brackets;
  int inexpr = 0, thislev, newlev;
  TERM *expressions[MAXEXPR], **exprp = expressions, **ep;
  TERM *tp, *last;
  FUN *funp;

  file = (*filename == '!') ? popen(filename + 1, "r") : fopen(filename, "r");
  if( file == NULL ) {
    report_error("noopen", filename);
    MLPutFunction(stdlink, "Abort", 0);
    MLEndPacket(stdlink);
    return;
  }

  termp = NULL;
  maxintsize = 0;

  for( ; ; ) {

nextline:
    do {
      if( MLAbort ) goto abort;
      if( feof(file) ) {
        if( erp > errmsg ) {
          *(erp - 1) = 0;	/* discard last \n */
          report_error("formerror", errmsg);
          goto abort;
        }
        inexpr = (int)(exprp - expressions);
        if( inexpr == 0 ) {
          report_error("nooutput", NULL);
          goto abort;
        }
        MLPutFunction(stdlink, "List", inexpr);
        for( ep = expressions; ep < exprp; ++ep ) {
          orderchain(*ep, LEVELS - 1);
          transmit(*ep, LEVELS - 1);
        }
        goto quit;
      }
      *line = 0;
      si = fgets(line, sizeof(line), file);
      if( (er = strstr(line, "-->")) ||
          (er = strstr(line, "==>")) ||
          (er = strstr(line, "===")) ) {
        er += 4;
        *erp = 0;
        if( !strstr(errmsg, er) &&
            (int)(erp - errmsg) + strlen(er) < sizeof(errmsg) )
          erp = memccpy(erp, er, '\n', sizeof(errmsg));
      }
      if( inexpr && *line == '\n' ) {
        *di++ = 0;
        termp = downsize(termp, di);
        goto newterm;
      }
      if( !inexpr && (er = strchr(line, '=')) ) {
newterm:
        tp = malloc(sizeof(TERM) + TERMBUF);
        beg = delim = di = (char *)tp + sizeof(TERM);
        tp->last = termp;
        termp = tp;
        for( pp = tp->f; pp < tp->f + LEVELS; ++pp ) *pp = zero;
        tp->f[thislev = LEVEL_COEFF] = di;
        tp->coll = 0;
        if( inexpr ) goto nextline;
        inexpr = 1;
        si = er + 1;
        break;
      }
    } while( !inexpr || erp > errmsg );

    for( ; *si; ++si)
      if( *si > ' ' ) switch(*si) {
      case '+':
      case '-':
      case '*':
        *di++ = *si;
        if( br == brackets ) delim = di;
        break;

      case '(':
        if( br == brackets ) {
          *di = 0;
          newlev = LEVEL_COEFF;
          for( funp = funtab;
               funp < &funtab[sizeof(funtab)/sizeof(FUN)];
               ++funp )
            if( strcmp(delim, funp->name) == 0 ) {
              newlev = funp->level;
              break;
            }
          if( thislev != newlev ) {
            if( delim > beg ) *(delim - 1) = 0;
            switch( thislev ) {
            case LEVEL_MAT:
              termp->f[LEVEL_MAT] = getabbr(termp->f[LEVEL_MAT]);
              break;
            case LEVEL_PAVE:
              maxintsize += (int)(delim - termp->f[LEVEL_PAVE]) + 2;
            }
            termp->f[thislev = newlev] = delim;
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
        *exprp++ = termp = downsize(termp, di);
        if( exprp >= expressions + MAXEXPR ) {
          report_error("toomany", NULL);
          goto abort;
        }
        if( maxintsize ) collect_pave();
        termp = NULL;
        maxintsize = inexpr = 0;
        goto nextline;

      case '_':
        if( *(di - 1) == 'i' ) *(di - 1) = 'I';
        else *di++ = '$';
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

abort:
  MLPutFunction(stdlink, "Abort", 0);
quit:
  MLEndPacket(stdlink);
  while( exprp > expressions )
    for( tp = *--exprp; tp; tp = last ) {
      if( tp->coll ) free(tp->f[LEVEL_PAVE]);
      last = tp->last;
      free(tp);
    }
  ((*filename == '!') ? pclose : fclose)(file);
}


void cutbranch(BTREE *node)
{
  if( node ) {
    cutbranch(node->lt);
    cutbranch(node->gt);
    free(node);
  }
}

void clearcache(void)
{
  cutbranch(root);
  root = NULL;
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}


main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

