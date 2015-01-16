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
		last modified 13 Jun 99 th

Note: This is the fancy version which performs some simplifications
      while sending the FORM output to Mma. Also in other respects
      (e.g. power counting of momenta), it is pretty much adapted to
      FormCalc. For a simple FORM -> Mma translator see FormGet.tm.
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXEXPR 2000
#define ITEMSIZE 32767

#ifndef MLCONST
#define MLCONST
#endif

#define FMEOFF ((char *)&dummy.fme - (char *)&dummy)
#define FACOFF ((char *)&dummy.fac - (char *)&dummy)
#define STR(p, off) *(char **)((char *)p + off)

typedef struct mathobj {
  struct mathobj *last;
  char *sme, *fme, *fac, *fac2;
  int sign, coll, nfac, nterms;
} MATHOBJ;

char *tok;
MATHOBJ *fp, *ep, dummy, *old1;
int infac, maxsmesize, maxintsize;

static char zero[] = "";

char *pctab[50] = { "p1", "p2", "p3", "k1", "k2", "k3" };
char **pctabend = &pctab[6];
char pcstore[256];


void report_error(MLCONST char *tag, const char *arg)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", arg ? 2 : 1);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "ReadForm");
  MLPutString(stdlink, tag);
  if(arg) MLPutString(stdlink, arg);
  MLPutSymbol(stdlink, "$Failed");
  MLEndPacket(stdlink);
}


char *getfac(char *s)
{
  int bracket = 0;

  for(tok = s; ; ++tok)
    switch(*tok) {
    case '(':
    case '[':
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


int powercount(char *s)
{
  char **pct;

  for(pct = pctab; pct < pctabend; ++pct)
    if(strcmp(s, *pct) == 0) return 1;
  return 0;
}


int chopup(char *si, MATHOBJ **running, char *wrap)
{
  char *psme, *pfac, *pfme, *pend, *di, *p, *sp, *c, *c2, h;
  int size, sc = 0;
  MATHOBJ *mp;

  size = strlen(si) + sizeof(MATHOBJ) + 1;
  if(*wrap) size += strlen(wrap) + 12;
  mp = malloc(size);
  mp->last = *running;
  *running = mp;

  mp->coll = 0;
  mp->sign = 1;
  switch(*si) {
  case '-':
    mp->sign = -1;
  case '+':
    ++si;
  }
  mp->fac2 = mp->fme = zero;
  pfme = NULL;
  pfac = pend = (char *)mp + size;
  mp->sme = (char *)mp + sizeof(MATHOBJ);
  p = psme = memccpy(mp->sme, wrap, 0, 128) - 1;

  for(di = getfac(si); *di; di = getfac(tok)) {
	/* dot products, eps tensors & one-loop functions -> psme */
    if(*wrap == 0 && strncmp(di, "fme", 3) == 0) {
      pfme = di;
      continue;
    }
    if(c = strchr(di, '[')) {	/* function */
      if(*di == 'e' && *(di + 1) == '$') {	/* eps tensor */
        sp = di;
        if(*wrap)	/* count powers of momenta */
          do {
            for(c2 = ++c; *c2 != ',' && *c2 != ']'; ++c2) ;
            h = *c2;
            *c2 = 0;
            sc += powercount(c);
            *(c = c2) = h;
          } while(h != ']');
      }
      else {
        h = *--c;
        *c = 0;
        sp = strstr("pave,B0,B1,B00,B11,A", di);
        *c = h;
      }
    }
    else {
      sp = strchr(di, '.');	/* dot product */
      if(sp && *wrap) {
        *sp = 0;
        sc += powercount(di) + powercount(sp + 1);
        *sp = '.';
      }
    }
    if(sp) {
      psme = memccpy(psme, di, 0, 256);
      *(psme - 1) = '*';
    }
    else {
      pfac -= tok - di;
      if(*tok == 0) --pfac;
      *((char *)memccpy(pfac, di, 0, 256) - 1) = '*';
    }
  }
  *(pend - 1) = 0;
  mp->fac = pfac - (pfac == pend);

  if(psme == p) mp->sme = zero;
  else {
    --psme;
    if(*wrap) {
      *psme++ = ']';
      if(sc) {
        *psme = 0;
        sprintf(mp->fac2 = ++psme, "$O1ME^%d", sc);
        psme += strlen(psme);
        if(sc == 1) psme -= 2;
      }
    }
    *psme++ = 0;
  }
  if(pfme) strcpy(mp->fme = psme, pfme);
  return *mp->sme ? psme - mp->sme : 2;
}


char *putfac(char *to, char *from)
{
  if(*from) return memccpy(to, from, 0, ITEMSIZE);
  *to++ = '1';
  *to++ = 0;
  return to;
}


void sum_up_fac2()
{
  MATHOBJ *f1p = fp, *f2p, *old;
  char *pool, *s, *s2;
  int overallsign, maxsumsize = maxsmesize + 5;

  infac = 0;
  pool = s = malloc(maxsmesize);
  do {
    for(old = f1p; f2p = old->last; )
      if(strcmp(f1p->fac, f2p->fac) ||
         strcmp(f1p->fac2, f2p->fac2)) old = f2p;
      else {
        if(f1p->coll == 0) {
          s2 = f1p->sme;
          s = putfac(f1p->sme = s, s2);
          f1p->coll = 1;
        }
        *(s - 1) = f1p->sign*f2p->sign < 0 ? '-' : '+';
        s = putfac(s, f2p->sme);
        old->last = f2p->last;
        free(f2p);
      }
    maxsumsize += strlen(f1p->fac) + 30;
  } while(f1p = f1p->last);

  ep->sign *= overallsign = fp->sign;
  ep->fac2 = s = malloc(maxsumsize);
  *s++ = 'o';
  *s++ = '2';
  *s++ = '[';
  do {
    *s++ = overallsign*fp->sign < 0 ? '-' : '+';
    *s++ = 'o';
    *s++ = '1';
    *s++ = '[';
    s = putfac(s, fp->fac);
    for(old = fp; f2p = old->last; )
      if(strcmp(fp->sme, f2p->sme)) old = f2p;
      else {
        *(s - 1) = fp->sign*f2p->sign < 0 ? '-' : '+';
        s = putfac(s, f2p->fac);
        old->last = f2p->last;
        free(f2p);
      }
    *(s - 1) = ']';
    if(*fp->fac2) {
      *s++ = '*';
      s = memccpy(s, fp->fac2, 0, ITEMSIZE) - 1;
    }
    if(*fp->sme) {
      *s++ = '*';  
      if(fp->coll) {
        strcpy(s, "smeplus[");
        s += 8;
      }
      s = memccpy(s, fp->sme, 0, ITEMSIZE) - 1;
      if(fp->coll) *s++ = ']';
    }
    old = fp->last;
    free(fp);
  } while(fp = old);
  free(pool);
  *s++ = ']';
  *s++ = 0;
  ep->fac2 = realloc(ep->fac2, s - ep->fac2);
}


void collect()
{
  MATHOBJ *f2p, *old;
  char *s;

  do {
    for(old = ep; f2p = old->last; )
      if(strcmp(ep->fac, f2p->fac) ||
         strcmp(ep->fac2, f2p->fac2) ||
         strcmp(ep->fme, f2p->fme)) old = f2p;
      else {
        if(ep->coll == 0) {
          s = ep->sme;
          s = putfac(ep->sme = malloc(maxintsize), s);
          ep->coll = 1;
        }
        *(s - 1) = ep->sign*f2p->sign < 0 ? '-' : '+';
        s = putfac(s, f2p->sme);
        old->last = f2p->last;
        free(f2p);
      }
    if(ep->coll) ep->sme = realloc(ep->sme, s - ep->sme);
  } while(ep = ep->last);
}


int orderchain(MATHOBJ *e1p, int off)
{
  MATHOBJ *e2p, *old2, *ini;
  int c = 0, c2;

  do {
    ++c;
    c2 = 0;
    ini = e1p;
    do {
      ++c2;
      old1 = e1p;
      e1p = old1->last;
      if(e1p == NULL) goto next;
    } while(strcmp(STR(ini, off), STR(e1p, off)) == 0);
    e2p = e1p;
    do {
      old2 = e2p;
over:
      e2p = old2->last;
      if(e2p == NULL) goto next;
    } while(strcmp(STR(ini, off), STR(e2p, off)));
    old1->last = e2p;
    old1 = e2p;
    old2->last = e2p->last;
    ++c2;
    goto over;
next:
    if(off == FACOFF) ini->nterms = c2;
    else {
      old1->last = NULL;
      ini->nfac = orderchain(ini, FACOFF);
    }
  } while(old1->last = e1p);
  return c;
}


void transmit(MATHOBJ *ep)
{
  MATHOBJ *last;
  int nfme, nfac, ntimes, nplus, overallsign, s;

  nfme = orderchain(ep, FMEOFF);
  if(nfme > 1) MLPutFunction(stdlink, "Plus", nfme);
  while(nfme--) {
    if(*ep->fme) {
      MLPutFunction(stdlink, "Times", 2);
      MLPutFunction(stdlink, "ToExpression", 1);
      MLPutString(stdlink, ep->fme);
    }
    nfac = ep->nfac;
    if(nfac > 1) MLPutFunction(stdlink, "Plus", nfac);
    while(nfac--) {
      overallsign = ep->sign;
      nplus = ep->nterms;
      ntimes = ((1 - overallsign) >> 1) + (*ep->fac != 0);
      if(s = nplus > 1) ++ntimes;
      else ntimes += (*ep->sme != 0) + (*ep->fac2 != 0);
      if(ntimes == 0) {
        MLPutInteger(stdlink, 1);
        goto noterm;
      }
      if(ntimes > 1) MLPutFunction(stdlink, "Times", ntimes);
      if(overallsign < 0) MLPutInteger(stdlink, -1);
      if(*ep->fac) {
        MLPutFunction(stdlink, "ToExpression", 1);
        MLPutString(stdlink, ep->fac);
      }
      if(s) MLPutFunction(stdlink, "Plus", nplus);
      while(nplus--) {
        if(s) {
          s = overallsign*ep->sign;
          ntimes = ((1 - s) >> 1) + (*ep->sme != 0) + (*ep->fac2 != 0);
          if(ntimes == 0) {
            MLPutInteger(stdlink, 1);
            goto noterm;
          }
          if(ntimes > 1) MLPutFunction(stdlink, "Times", ntimes);
          if(s < 0) MLPutInteger(stdlink, -1);
        }
        if(*ep->sme) {
          MLPutFunction(stdlink, "ToExpression", 1);
          MLPutString(stdlink, ep->sme);
        }
        if(*ep->fac2) {
          MLPutFunction(stdlink, "ToExpression", 1);
          MLPutString(stdlink, ep->fac2);
        }
noterm:
        last = ep->last;
        free(ep->fac2);
        if(ep->coll) free(ep->sme);
        free(ep);
        ep = last;
      } /* sme*fac2 */
    } /* fac */
  } /* fme */
}


void expr_as_fac(char *expr)
{
  if(!infac) {
    maxsmesize = 0;
    fp = malloc(sizeof(MATHOBJ));
    fp->sme = fp->fac = fp->fme = zero;
    fp->sign = 1;
    fp->coll = 0;
    fp->last = ep;
    ep = fp;
    fp = NULL;
    infac = 1;
  }
  maxsmesize += chopup(expr, &fp, "sme[");
}


void readform(const char *filename)
{
  FILE *file;
  char physicalline[280], logicalline[2048];
  char *si, *di = logicalline;
  char *er, errmsg[512], *erp = errmsg;
  int brackets[16], *br = brackets;
  int inexpr = 0;
  MATHOBJ *delim[MAXEXPR], **dp = delim, **dp2;

  file = (*filename == '!' ?
    popen(filename + 1, "r") : fopen(filename, "r"));
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
            for(dp2 = delim; dp2 < dp; ++dp2) transmit(*dp2);
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
          maxsmesize = 0;
          maxintsize += chopup(di = logicalline, &ep, "");
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
          maxsmesize += chopup(di = logicalline, &fp, "sme[");
          sum_up_fac2();
        }
        else *di++ = *--br ? ')' : ']';
        break;
      case ';':
        if(di > logicalline) {
          *di = 0;
          expr_as_fac(di = logicalline);
        }
        if(infac) sum_up_fac2();
        *dp++ = ep;
        if(dp >= delim + MAXEXPR) {
          report_error("toomany", NULL);
          return;
        }
        collect();
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
      pcs = memccpy(*pct++ = pcs, arg, 0, 128);
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

