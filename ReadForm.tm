:Begin:
:Function:	readform
:Pattern:	ReadForm[filename_String]
:Arguments:	{filename}
:ArgumentTypes:	{String}
:ReturnType:	Manual
:End:

:Begin:
:Function:	powercountingfor
:Pattern:	PowerCountingFor[mom___]
:Arguments:	{mom}
:ArgumentTypes:	{Manual}
:ReturnType:	Manual
:End:

:Evaluate:	ReadForm::notfound = "File not found."
:Evaluate:	ReadForm::nooutput =
		"Something went wrong, there was no output from FORM."
:Evaluate:	ReadForm::unknown =
		"Unimplemented higher tensor functions found."
:Evaluate:	ReadForm::toomany =
		"Too many expressions. Increase MAXEXPR in ReadForm.tm."
:Evaluate:	ReadForm::formerror = "`1`"


/*
	ReadForm.tm
		reads FORM output back into Mathematica
		this file is part of FormCalc
		last modified 12 Jan 00 th

Note: This is the fancy version which performs some simplifications
      while sending the FORM output to Mma. Also in other respects
      (e.g. power counting of momenta), it is pretty much adapted to
      FormCalc. For a simple FORM -> Mma translator see FormGet.tm.
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char copyleft[] =
  "@(#) ReadForm utility for FormCalc, 12 Jan 00 Thomas Hahn";

#define MAXEXPR 2000
#define STRINGSIZE 32767

#ifndef MLCONST
#define MLCONST
#endif


/*

A term in the FORM output is organized into the MATHOBJ structure
in the following way:

 ____5_____     __4___     __3___     _____0______     _2_
/          \   /      \   /      \   /            \   /   \
SumOver(...) * Mat(...) * DEN(...) * A0/B0/...(...) * ..... * (
                                                 _
    + 3/16*EL^2 * abb(...) [ * $Scale^n ]         |
      \_     _/   \_   __/ \_         __/         | 1
        coeff       abb      scale = n           _|
    + ...
)

The lines inside the parentheses are referred to as "factors", the whole
term as "expression".

Hierarchy of collecting:
5. SumOver
4. Mat
3. DEN
2. [coefficient]
1. [factor]
0. [integral]

*/

#define LEVELS 6
#define LEVEL_INTEGRAL 0
#define LEVEL_FACTOR 1
#define LEVEL_COEFFICIENT 2
#define LEVEL_MAT 4

typedef struct {
  char *func;
  int level;
} COLL;

COLL collecttable[] = {
  {"SumOver", 5},
  {"Mat", LEVEL_MAT},
  {"DEN", 3},
  {"A0", LEVEL_INTEGRAL},
  {"B0m", LEVEL_INTEGRAL},
  {"B1m", LEVEL_INTEGRAL},
  {"B00m", LEVEL_INTEGRAL},
  {"B11m", LEVEL_INTEGRAL},
  {"pave3", LEVEL_INTEGRAL},
  {"pave4", LEVEL_INTEGRAL},
  {"pave5", LEVEL_INTEGRAL}
};


typedef struct exprterm {
  struct exprterm *last;
  char *f[LEVELS];
  int nterms[LEVELS], sign, coll;
} EXPRTERM;

typedef struct facterm {
  struct facterm *last;
  char *abb;
  int sign, scale;
  char coeff[8];
} FACTERM;

typedef struct btree {
  struct btree *lt, *gt;
  char *sym;
  int scale;
  char abb[8];
} BTREE;

/* strict ANSI C doesn't allow incomplete arrays in structs, so the
   "incomplete" arrays actually have 8 bytes (8 to avoid unnecessary
   padding), and to compensate we need */
#define sizeofstruct(x) (sizeof(x) - 8)


char *tok;
FACTERM *fp;
EXPRTERM *ep, *old1;
BTREE *abb_root = NULL, *abbsum_root = NULL, *mat_root = NULL;
int infac, maxabbsize, maxintsize;

static char zero[] = "";


char *pctab[50], **pctabend = pctab;
char pcstore[256];


void report_error(MLCONST char *tag, const char *arg)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", arg ? 2 : 1);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "ReadForm");
  MLPutString(stdlink, tag);
  if(arg) MLPutString(stdlink, arg);
  MLPutFunction(stdlink, "Abort", 0);
  MLEndPacket(stdlink);
}


char *getfac(char *s)
{
  int bracket = 0;

  for(tok = s; ; ++tok)
    switch(*tok) {
    case '[':
    case '(':
      ++bracket;
      break;
    case ')':
    case ']':
      --bracket;
      break;
    case '*':
      if(bracket) continue;
      *tok++ = 0;
    case 0:
      return s;
    }
}


char *putfac(char *to, char *from)
{
  if(*from) return memccpy(to, from, 0, STRINGSIZE);
  *to++ = '1';
  *to++ = 0;
  return to;
}


#define nonalpha(c) ((unsigned char)((c | 0x20) - 'a') > 25)

int powercount(char *s)
{
  char **pct, *p, *t;
  int c = 0;

  for(pct = pctab; pct < pctabend; ++pct)
    for(p = s; p = strstr(p, *pct); ) {
      t = p - 1;
      p += strlen(*pct);
      if(nonalpha(*t) && nonalpha(*p) &&
          !(*t == '[' && memcmp(t - 6, "Spinor", 6))) ++c;
    }
  return c;
}


char *getabbr(char *head, char *s, BTREE **root, int *count)
{
  MLCONST char *mmares;
  BTREE *lp;
  int p;

  while(lp = *root) {
    p = strcmp(s, lp->abb);
    if(p == 0) {
      if(count) *count = lp->scale;
      return lp->sym;
    }
    root = p < 0 ? &lp->lt : &lp->gt;
  }

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "ToString", 1);
  MLPutFunction(stdlink, head, 1);
  MLPutString(stdlink, s);
  MLEndPacket(stdlink);

  while((p = MLNextPacket(stdlink)) && p != RETURNPKT)
    MLNewPacket(stdlink);
  MLGetString(stdlink, &mmares);

  p = strlen(s);
  lp = malloc(sizeofstruct(BTREE) + 4 + p + strlen(mmares));
			/* 4 = strlen("()\0\0") */
  lp->lt = lp->gt = NULL;
  *root = lp;
  strcpy(lp->abb, s);
  if(count) *count = lp->scale = powercount(s);
  s = lp->sym = lp->abb + p + 1;
  if(strpbrk(mmares, "+-")) sprintf(s, "(%s)", mmares);
  else strcpy(s, mmares);

  MLDisownString(stdlink, mmares);
  return s;
}


void fac_chopup(char *si)
{
  char *pabb;
  int size;
  FACTERM *mp;

  size = strlen(si) + 1 + sizeofstruct(FACTERM);
  if(pabb = strstr(si, "abb[")) {
    getfac(pabb + 3);
    size -= (int)(tok - pabb);
  }
  mp = malloc(size);
  mp->last = fp;
  fp = mp;

  mp->sign = 1;
  switch(*si) {
  case '-':
    mp->sign = -1;
  case '+':
    ++si;
  }

  if(pabb) {
    mp->abb = getabbr("FromForm", pabb, &abb_root, &mp->scale);
    maxabbsize += strlen(mp->abb) + 1;
    if(size = pabb - si) {
      if(*tok == 0) --size;		/* remove trailing "*" */
      memcpy(mp->coeff, si, size);
    }
    strcpy(mp->coeff + size, tok);
  }
  else {
    mp->scale = 0;
    mp->abb = zero;
    strcpy(mp->coeff, si);
    maxabbsize += 2;
  }
}


void sum_up_fac()
{
  FACTERM *f1p = fp, *f2p, *old;
  char *pool, *s;
  int overallsign, maxsumsize = 5;	/* 5 = strlen("o2[]\0") */

  infac = 0;
  pool = malloc(maxabbsize + 9);	/* 9 = strlen("abbsum[]\0") */
  do {
    s = NULL;
    for(old = f1p; f2p = old->last; )
      if(strcmp(f1p->coeff, f2p->coeff) ||
         f1p->scale != f2p->scale) old = f2p;
      else {
        if(s == NULL) {
          strcpy(s = pool, "abbsum[");
          s = putfac(s + 7, f1p->abb);
        }
        *(s - 1) = f1p->sign*f2p->sign < 0 ? '-' : '+';
        s = putfac(s, f2p->abb);
        old->last = f2p->last;
        free(f2p);
      }
    if(s) {
      *(s - 1) = ']';
      *s = 0;
      f1p->abb = getabbr("ToExpression", pool, &abbsum_root, NULL);
    }
    maxsumsize += strlen(f1p->coeff) + strlen(f1p->abb) + 16;
			/* 16 = strlen("+o1[]**$Scale^nn") */
  } while(f1p = f1p->last);
  free(pool);

  ep->sign *= overallsign = fp->sign;
  ep->f[LEVEL_FACTOR] = s = malloc(maxsumsize);
  *s++ = 'o';
  *s++ = '2';
  *s++ = '[';
  do {
    *s++ = overallsign*fp->sign < 0 ? '-' : '+';
    *s++ = 'o';
    *s++ = '1';
    *s++ = '[';
    s = putfac(s, fp->coeff);
    for(old = fp; f2p = old->last; )
      if(strcmp(fp->abb, f2p->abb)) old = f2p;
      else {
        *(s - 1) = fp->sign*f2p->sign < 0 ? '-' : '+';
        s = putfac(s, f2p->coeff);
        old->last = f2p->last;
        free(f2p);
      }
    *(s - 1) = ']';
    if(fp->scale) {
      s += sprintf(s, "*$Scale^%d", fp->scale);
      if(fp->scale == 1) s -= 2;
    }
    if(*fp->abb) {
      *s++ = '*';
      s = memccpy(s, fp->abb, 0, STRINGSIZE) - 1;
    }
    old = fp->last;
    free(fp);
  } while(fp = old);
  *s++ = ']';
  *s++ = 0;
  ep->f[LEVEL_FACTOR] =
    realloc(ep->f[LEVEL_FACTOR], s - ep->f[LEVEL_FACTOR]);
}


void expr_chopup(char *si)
{
  char *pfunc, *pfunc0, *pcoeff, *pcoeff0, *di, *c, **pp, s[20];
  COLL *cp;
  int size, scale = 0, thislev = -1;
  EXPRTERM *mp;

  size = strlen(si) + sizeof(EXPRTERM) + 1;
  mp = malloc(size);
  mp->last = ep;
  ep = mp;

  mp->coll = 0;
  mp->sign = 1;
  switch(*si) {
  case '-':
    mp->sign = -1;
  case '+':
    ++si;
  }
  for(pp = mp->f; pp < mp->f + LEVELS; ++pp) *pp = zero;
  pcoeff = pcoeff0 = (char *)mp + size;
  pfunc = pfunc0 = (char *)mp + sizeof(EXPRTERM);

  for(di = getfac(si); *di; di = getfac(tok)) {
    if(c = strchr(di, '[')) {
      *c = 0;
      for(cp = collecttable;
          cp < &collecttable[sizeof(collecttable)/sizeof(COLL)];
          ++cp)
        if(strcmp(di, cp->func) == 0) {
          *c = '[';
          if(cp->level == LEVEL_MAT)	/* note: if there's more than one
					   Mat function per line, we're
					   in trouble here */
            mp->f[thislev = LEVEL_MAT] =
              getabbr("FromForm", di, &mat_root, &scale);
          else {
            if(thislev == cp->level) *(pfunc - 1) = '*';
            else mp->f[thislev = cp->level] = pfunc;
            pfunc = memccpy(pfunc, di, 0, STRINGSIZE);
          }
          goto loop;
        }
      *c = '[';
    }
    thislev = -1;
    pcoeff -= tok - di + (*tok == 0);
    *((char *)memccpy(pcoeff, di, 0, STRINGSIZE) - 1) = '*';
loop: ;
  }

  if(scale) {
    if(scale == 1) memcpy(pcoeff -= 7, "$Scale*", 7);
    else {
      size = sprintf(s, "$Scale^%d*", scale);
      memcpy(pcoeff -= size, s, size);
    }
  }
  if(pcoeff != pcoeff0) {
    *(pcoeff0 - 1) = 0;
    mp->f[LEVEL_COEFFICIENT] = pcoeff;
  }
  maxintsize += strlen(mp->f[LEVEL_INTEGRAL]) + 2;
}


void collect_integrals()
{
  EXPRTERM *f2p, *old;
  char *s;
  int i;

  do {
    for(old = ep; f2p = old->last; ) {
      for(i = LEVEL_INTEGRAL + 1; i < LEVELS; ++i)
        if(strcmp(ep->f[i], f2p->f[i])) {
          old = f2p;
          goto loop;
        }
      if(ep->coll == 0) {
        s = ep->f[LEVEL_INTEGRAL];
        s = putfac(ep->f[LEVEL_INTEGRAL] = malloc(maxintsize), s);
        ep->coll = 1;
      }
      *(s - 1) = ep->sign*f2p->sign < 0 ? '-' : '+';
      s = putfac(s, f2p->f[LEVEL_INTEGRAL]);
      old->last = f2p->last;
      free(f2p);
loop: ;
    }
    if(ep->coll) ep->f[LEVEL_INTEGRAL] =
      realloc(ep->f[LEVEL_INTEGRAL], s - ep->f[LEVEL_INTEGRAL]);
  } while(ep = ep->last);
}


void orderchain(EXPRTERM *e1p, int level)
{
  EXPRTERM *e2p, *old2, *ini;
  int c = 0, c2, *nterms = &e1p->nterms[level];

  do {
    ++c;
    c2 = 0;
    ini = e1p;
    do {
      ++c2;
      old1 = e1p;
      e1p = old1->last;
      if(e1p == NULL) goto next;
    } while(strcmp(ini->f[level], e1p->f[level]) == 0);
    e2p = e1p;
    do {
      old2 = e2p;
over:
      e2p = old2->last;
      if(e2p == NULL) goto next;
    } while(strcmp(ini->f[level], e2p->f[level]));
    old1->last = e2p;
    old1 = e2p;
    old2->last = e2p->last;
    ++c2;
    goto over;
next:
    if(level > LEVEL_COEFFICIENT) {
      old1->last = NULL;
      orderchain(ini, level - 1);
    }
    else ini->nterms[level - 1] = c2;
  } while(old1->last = e1p);
  *nterms = c;
}


void transmit(EXPRTERM **xp, int level)
{
  EXPRTERM *ep = *xp, *last;
  int n = ep->nterms[level], ntimes, i;

  if(n > 1) MLPutFunction(stdlink, "Plus", n);
  while(n--) {
    ntimes = *ep->f[level] != 0;
    for(i = level - 1; i > LEVEL_INTEGRAL; --i) {
      ++ntimes;
      if(ep->nterms[i] > 1) goto sendit;
      if(*ep->f[i] == 0) --ntimes;
    }
	/* orderchain goes down only to LEVEL_FACTOR, hence: */
    if(*ep->f[LEVEL_INTEGRAL]) ++ntimes;
    if(ep->sign < 0) ++ntimes;
sendit:
    switch(ntimes) {
    case 0:
      MLPutInteger(stdlink, 1);
      break;
    default:
      MLPutFunction(stdlink, "Times", ntimes);
    case 1:
      if(*ep->f[level]) {
        MLPutFunction(stdlink, "ToExpression", 1);
        MLPutString(stdlink, ep->f[level]);
      }
      for(i = level - 1; i > LEVEL_INTEGRAL; --i) {
        if(ep->nterms[i] > 1) {
          transmit(&ep, i);
          goto loop;
        }
        if(*ep->f[i]) {
          MLPutFunction(stdlink, "ToExpression", 1);
          MLPutString(stdlink, ep->f[i]);
        }
      }
      if(*ep->f[LEVEL_INTEGRAL]) {
        MLPutFunction(stdlink, "ToExpression", 1);
        MLPutString(stdlink, ep->f[LEVEL_INTEGRAL]);
      }
      if(ep->sign < 0) MLPutInteger(stdlink, -1);
    }
    if(ep->coll) free(ep->f[LEVEL_INTEGRAL]);
    free(ep->f[LEVEL_FACTOR]);
    last = ep->last;
    free(ep);
    ep = last;
loop: ;
  }
  *xp = ep;
}


void expr_as_fac(char *expr)
{
  EXPRTERM *np;
  char **pp;

  if(!infac) {
    maxabbsize = 0;
    np = malloc(sizeof(EXPRTERM));
    for(pp = np->f; pp < np->f + LEVELS; ++pp) *pp = zero;
    np->sign = 1;
    np->coll = 0;
    np->last = ep;
    ep = np;
    fp = NULL;
    infac = 1;
  }
  fac_chopup(expr);
}


void readform(const char *filename)
{
  FILE *file;
  char physicalline[280], logicalline[2048];
  char *si, *di = logicalline;
  char *er, errmsg[512], *erp = errmsg;
  int brackets[16], *br = brackets;
  int inexpr = 0;
  EXPRTERM *delim[MAXEXPR], **dp = delim, **dp2;

  file = *filename == '!' ?
    popen(filename + 1, "r") : fopen(filename, "r");
  if(file == NULL) {
    report_error("notfound", NULL);
    return;
  }

  ep = NULL;
  maxintsize = infac = 0;

  for( ; ; ) {
nextline:
    do {
      if(feof(file)) {
        if(erp == errmsg) {
          (*filename == '!' ? pclose : fclose)(file);
          inexpr = (int)(dp - delim);
          if(inexpr == 0) report_error("nooutput", NULL);
          else {
            if(inexpr != 1) MLPutFunction(stdlink, "List", inexpr);
            for(dp2 = delim; dp2 < dp; ++dp2) {
              orderchain(*dp2, LEVELS - 1);
              transmit(dp2, LEVELS - 1);
            }
            MLEndPacket(stdlink);
          }
          return;
        }
        *(erp - 1) = 0;		/* discard last \n */
        report_error("formerror", errmsg);
        return;
      }
      *physicalline = 0;
      fgets(si = physicalline, sizeof(physicalline), file);
      er = physicalline - 2;
      if(*physicalline > ' ' || (er = strchr(physicalline, '>'))) {
        er += 2;
        *erp = 0;
        if(!strstr(errmsg, er) &&
             (int)(erp - errmsg) + strlen(er) < sizeof(errmsg))
          erp = memccpy(erp, er, '\n', sizeof(errmsg));
      }
      if(!inexpr && (er = strchr(physicalline, '='))) {
        si = er + 1;
        inexpr = 1;
      }
    } while(!inexpr || erp > errmsg);
    while(*si == ' ') ++si;
    if(strstr(si, "unknown")) {
      report_error("unknown", NULL);
      return;
    }
    for( ; *si; ++si)
      if(*si > ' ') switch(*si) {
      case '_':
        if(*(di - 1) == 'i') *(di - 1) = 'I';
        else *di++ = '$';
        break;
      case '?':
        break;
      case '*':
        *di++ = br == brackets ? '*' : ' ';
        break;
      case '(':
        if(*(di - 1) == '*') {
          *(di - 1) = 0;
          fp = NULL;
          maxabbsize = 0;
          expr_chopup(di = logicalline);
          infac = 1;
        }
        else {
          *br = *(si - 1) < '0' || *(si - 1) == '^';
          *di++ = *br++ ? '(' : '[';
        }
        break;
      case ')':
        if(infac && br == brackets) {
          *di = 0;
          fac_chopup(di = logicalline);
          sum_up_fac();
        }
        else *di++ = *--br ? ')' : ']';
        break;
      case ';':
        if(di > logicalline) {
          *di = 0;
          expr_as_fac(di = logicalline);
        }
        if(infac) sum_up_fac();
        *dp++ = ep;
        if(dp >= delim + MAXEXPR) {
          report_error("toomany", NULL);
          return;
        }
        collect_integrals();
        maxintsize = inexpr = 0;
        goto nextline;
      case '+':
      case '-':
        if(di > logicalline && *(di - 1) != '^' && br == brackets) {
          *di = 0;
          expr_as_fac(di = logicalline);
        }
      default:
        *di++ = *si;
        break;
      }
  }
}


void powercountingfor(void)
{
  long argc;
  const char *arg;
  char *pcs = pcstore, **pct = pctab;

  if(MLReady(stdlink)) {
    MLGetFunction(stdlink, &arg, &argc);
    MLDisownSymbol(stdlink, arg);
    while(argc--) {
      MLGetSymbol(stdlink, &arg);
      pcs = memccpy(*pct++ = pcs, arg, 0, STRINGSIZE);
      MLDisownSymbol(stdlink, arg);
    }
    pctabend = pct;
    pct = pctab;
  }
  MLPutFunction(stdlink, "List", (int)(pctabend - pctab));
  while(pct < pctabend) MLPutSymbol(stdlink, *pct++);
  MLEndPacket(stdlink);
}


main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

