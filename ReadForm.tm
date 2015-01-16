:Begin:
:Function:	readform
:Pattern:	ReadForm[filename_String]
:Arguments:	{filename}
:ArgumentTypes:	{String}
:ReturnType:	Manual
:End:

:Evaluate:	ReadForm::notfound = "File not found."
:Evaluate:	ReadForm::unknown =
		"Unimplemented higher tensor functions found."
:Evaluate:	ReadForm::formerror = "`1`"

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#define MAXFACS 1500
#define MAXEXPR 1200
#define EXPRSPACESIZE 1500000
#define FACSIZE 128
#define ITEMSIZE 32767
#define FACSPACESIZE (2 * MAXFACS * FACSIZE)

#ifndef MLCONST
#define MLCONST
#endif

typedef struct mathobj {
  int sign, coll;
  char *sme, *fac2;
  union {
    char *fac;
    struct mathobj *next;
  } f;
} MATHOBJ;

char theinput[256], errmsg[512], *tok;
char exprspace[EXPRSPACESIZE];
char facspace[FACSPACESIZE];
MATHOBJ expr[MAXEXPR], facs[MAXFACS];

#ifdef DEBUG
int chkEXPRSPACESIZE = 0, chkMAXEXPR;
int chkFACSPACESIZE = 0, chkMAXFACS = 0;
#define BOUND(var_, expr_, desc_, max_) \
  if(expr_ > var_) \
    fprintf(stderr, desc_ " = %d(%d)\n", var_ = expr_, max_)
#else
#define BOUND(var_, expr_, desc_, max_)
#endif

void report_error(MLCONST char *tag, int n)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", n);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "ReadForm");
  MLPutString(stdlink, tag);
  if(n == 2) MLPutString(stdlink, errmsg);
  MLPutSymbol(stdlink, "$Failed");
  MLEndPacket(stdlink);
}

int getline(FILE *fi)
{
  char *si, *di, s[256];
  static int brackets[16], *br = brackets, brcount = 0;
  int errlen = 0, nu;

  do {
    do {
      if(feof(fi)) {
        if(errlen == 0) return 0;
        *(errmsg + errlen - 1) = 0;
        report_error("formerror", 2);
        return -1;
      }
      *s = 0;
      fgets(si = s, sizeof(s), fi);
      di = s - 2;
      if(*s > ' ' || (di = strchr(s, '>'))) {
        if((nu = errlen + strlen(di += 2)) < sizeof(errmsg)) {
          strcpy(errmsg + errlen, di);
          errlen = nu;
        }
      }
    } while(errlen);
    for(di = theinput; *si; ++si)
      if(*si > ' ') switch(*si) {
      case '_':
        *di++ = '$';
        break;
      case '(':
        if(*br++ = (si == s || *(si - 1) < '0' || *(si - 1) == '^'))
          *di++ = '(';
        else {
          *di++ = '[';
          ++brcount;
        }
        break;
      case '*':
        *di++ = (brcount == 0 ? '*' : ' ');
        break;
      case ')':
        *di++ = (br == brackets || *--br ? ')' : (--brcount, ']'));
      case ';':
      case '?':
        break;
      default:
        *di++ = *si;
        break;
      }
  } while(di == theinput);
  *di = 0;
  if(strstr(theinput, "unknown")) {
    report_error("unknown", 1);
    return -1;
  }
  return 1;
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

void chopup(char *si, MATHOBJ *mp, char *wrap, char **put)
{
  char *psme, *pfac, *di, *p, **sp, *c;
  int sc;

  mp->coll = 0;
  mp->sign = 1;
  switch(*si) {
  case '-':
    mp->sign = -1;
  case '+':
    ++si;
  }
  mp->f.fac = pfac = *put;
  mp->sme = pfac + strlen(si);
  p = psme = memccpy(mp->sme, wrap, 0, 128) - 1;
  for(di = getfac(si); *di; di = getfac(tok)) {
	/* dot products, eps tensors & one-loop functions -> psme */
    if(strchr(di, '.')) sp = &psme;
    else if(!(c = strchr(di, '['))) sp = &pfac;
    else {
      *c = 0;
      sp = strstr("e$,fme,pave4,pave3,B0m,B1m,B00m,B11m,A0", di) ?
        &psme : &pfac;
      *c = '[';
    }
    *sp = memccpy(*sp, di, 0, 256);
    *(*sp - 1) = '*';
  }
  if(pfac > *put) --pfac;
  *(mp->fac2 = pfac) = 0;
  if(psme == p) mp->sme = psme = pfac;
  else {
    --psme;
    if(*wrap) {
      *psme++ = ']';
      for(sc = 0; p < psme; ++p)
        if(strncmp(p, "Spinor", 6) == 0) p += 8;
        else switch(*p) {
          case 'e':
            ++p;
            break;
          case 'k':
          case 'p':
            if(*(p + 1) >= '0' && *(p + 1) <= '9') ++sc;
          }
      if(sc) {
        *psme = 0;
        sprintf(mp->fac2 = ++psme, "$O1ME^%d", sc);
        psme += strlen(psme);
        if(sc == 1) psme -= 2;
      }
    }
    *psme = 0;
  }
  *put = psme + 1;
}

char *putfac(char *to, char *from)
{
  if(*from) return memccpy(to, from, 0, ITEMSIZE);
  *to++ = '1';
  *to++ = 0;
  return to;
}

void collect(MATHOBJ *f1p, MATHOBJ *last, char *s)
{
  MATHOBJ *f2p;
  char *s2;

  for( ; f1p < last; ++f1p) {
    if(f1p->coll < 0) continue;
    for(f2p = f1p + 1; f2p < last; ++f2p)
      if(f2p->coll >= 0 && strcmp(f1p->f.fac, f2p->f.fac) == 0 &&
          strcmp(f1p->fac2, f2p->fac2) == 0) {
        if(!f1p->coll) {
          s2 = f1p->sme;
          s = putfac(f1p->sme = s, s2);
        }
        *(s - 1) = (f1p->sign * f2p->sign < 0 ? '-' : '+');
        s = putfac(s, f2p->sme);
        f2p->coll = -1;
        f1p->coll = 1;
      }
  }
  BOUND(chkFACSPACESIZE, (int)(s - facspace),
    "FACSPACESIZE", FACSPACESIZE);
}

void gather(MATHOBJ *first, MATHOBJ *last, char **put)
{
  MATHOBJ *f1p, *f2p;
  char *s = *put;

  *s++ = 'o';
  *s++ = '2';
  *s++ = '[';
  for(f1p = first; f1p < last; ++f1p) {
    if(f1p->coll < 0) continue;
    *s++ = (first->sign * f1p->sign < 0 ? '-' : '+');
    *s++ = 'o';
    *s++ = '1';
    *s++ = '[';
    s = putfac(s, f1p->f.fac);
    for(f2p = f1p+1; f2p < last; ++f2p)
      if(f2p->coll >= 0 && strcmp(f1p->sme, f2p->sme) == 0 &&
          strcmp(f1p->fac2, f2p->fac2) == 0) {
        *(s - 1) = (f1p->sign * f2p->sign < 0 ? '-' : '+');
        s = putfac(s, f2p->f.fac);
        f2p->coll = -1;
      }
    *(s - 1) = ']';
    if(*f1p->fac2) {
      *s++ = '*';
      s = memccpy(s, f1p->fac2, 0, ITEMSIZE) - 1;
    }
    if(*f1p->sme) {
      *s++ = '*';  
      if(f1p->coll) {
        strcpy(s, "smeplus[");
        s += 8;
      }
      s = memccpy(s, f1p->sme, 0, ITEMSIZE) - 1;
      if(f1p->coll) *s++ = ']';
    }
  }
  *s++ = ']';
  *s++ = 0;
  *put = s;
  BOUND(chkEXPRSPACESIZE, (int)(s - exprspace),
    "EXPRSPACESIZE", EXPRSPACESIZE);
}

void readform(const char *filename)
{
  FILE *file;
  char *si, *di, *esp = exprspace, *fsp;
  int infac = 0, l, s;
  MATHOBJ *ep = expr, *fp, *xp, *x2p, **pp;

  file = (*filename == '!' ?
    popen(filename + 1, "r") : fopen(filename, "r"));
  if(file == NULL) {
    report_error("notfound", 1);
    return;
  }
  do {
    switch(getline(file)) {
    case 0:
      MLPutInteger(stdlink, 0);
      MLEndPacket(stdlink);
    case -1:
      return;
    }
  } while((si = strchr(theinput, '=')) == NULL);
  if(*++si) goto thisline;

  for( ; ; ) {
    l = getline(file);
    if(l == 0) {
      if(infac) goto finishfac;
      break;
    }
    if(l == -1) return;
    si = theinput;
thisline:
    if(infac) {
      if(*si == ')' || *(di = si + strlen(si) - 1) == '(') {
        infac = 0;
finishfac:
        BOUND(chkMAXFACS, (int)(fp - facs), "MAXFACS", MAXFACS);
        collect(facs, fp, fsp);
        ep->fac2 = esp;
        gather(facs, fp, &esp);
        ep->sign *= facs->sign;
        ++ep;
#ifdef DEBUG
	if((int)(ep - expr) % 100 == 0)
          fprintf(stderr, "MAXEXPR = %d(%d)\n", (int)(ep - expr), MAXEXPR);
#endif
        if(l == 0) break;
        if(*si != ')') goto isexpr;
        continue;
      }
isfac:
      chopup(si, fp++, "sme[", &fsp);
    }
    else {
      di = si + strlen(si) - 1;
isexpr:
      ep->sign = infac = 1;
      fp = facs;
      fsp = facspace;
      if(*di != '(') {
        ep->sme = ep->f.fac = esp;
        ep->coll = 0;
        *esp++ = 0;
        goto isfac;
      }
      *(di - 1) = 0;
      chopup(si, ep, "", &esp);
    }
  }
  (*filename == '!' ? pclose : fclose)(file);
  BOUND(chkMAXEXPR, (int)(ep - expr), "MAXEXPR", MAXEXPR);
  collect(expr, ep, facspace);
  for(l = 0, xp = expr; xp < ep; ++xp) {
    if(xp->coll < 0) continue;
    ++l;
    xp->coll = 1;
    if(*xp->f.fac == 0) continue;
    for(pp = &xp->f.next, x2p = xp + 1; x2p < ep; ++x2p)
      if(x2p->coll >= 0 && strcmp(x2p->f.fac, (char *)*pp) == 0) {
        ++xp->coll;
        x2p->coll = -1;
        *pp = x2p;
        pp = &x2p->f.next;
      }
  }
  if(l > 1) MLPutFunction(stdlink, "Plus", l);
  for(xp = expr; xp < ep; ++xp) {
    if(xp->coll < 0) continue;
    l = (1 - xp->sign) >> 1;
    if(infac = (xp->coll > 1)) l += 2;
    else l += (*xp->f.fac != 0) + (*xp->sme != 0) + (*xp->fac2 != 0);
    if(l > 1) MLPutFunction(stdlink, "Times", l);
    if(xp->sign < 0) MLPutInteger(stdlink, -1);
    if(infac) MLPutFunction(stdlink, "Plus", xp->coll);
    for(x2p = xp; ; x2p = x2p->f.next) {
      s = (1 - xp->sign * x2p->sign) >> 1;
      if(infac) {
        l = s + (*x2p->sme != 0) + (*x2p->fac2 != 0);
        if(l > 1) MLPutFunction(stdlink, "Times", l);
      }
      if(s) MLPutInteger(stdlink, -1);
      if(*x2p->sme) {
        MLPutFunction(stdlink, "ToExpression", 1);
        MLPutString(stdlink, x2p->sme);
      }
      if(*x2p->fac2) {
        MLPutFunction(stdlink, "ToExpression", 1);
        MLPutString(stdlink, x2p->fac2);
      }
      if(--xp->coll == 0) break;
    }
    if(*x2p->f.fac) {
      MLPutFunction(stdlink, "ToExpression", 1);
      MLPutString(stdlink, x2p->f.fac);
    }
  }
  MLEndPacket(stdlink);
}

main(int argc, char **argv)
{
  return MLMain(argc, argv);
}
