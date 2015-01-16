/*
r8.c
replaces all real constants by Fortran-style double precision
numbers (1.234D0) in Mma FortranForm output
this file is part of FormCalc
last modified 28 Feb 00 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

main(int argc, char **argv)
{
  char s[200], next[200], buf[200], *p, *si, *di, *num;
  int gotnext = 0, didexp;

  while(!feof(stdin)) {
    if(gotnext) {
      strcpy(s, next);
      gotnext = 0;
    }
    else {
      *s = 0;
      fgets(s, sizeof(s), stdin);
    }
    p = s + strlen(s) - 1;
    *p-- = 0;
    if(*p == '\\') {
      fgets(next, sizeof(next), stdin);
      di = next + 6 + strspn(next + 6, " \t");
      si = di + strcspn(di, " */()");
      memcpy(p, di, (int)(si - di));
      *(p + (int)(si - di)) = 0;
      strcpy(di, si);
      gotnext = 1;
    }
    p = s;
    while(p = strpbrk(p, "0123456789")) {
      if(*(p - 1) >= 'A') {
        p += strspn(p, "0123456789");
        continue;
      }
      strtod(num = p, &p);
      if(memchr(num, '.', p - num)) {
        strcpy(buf, p);
        didexp = 0;
        for(di = si = num; si < p; ++si)
          if(*si != 'e') *di++ = *si;
          else {
            while(di > num && *(di - 1) == '0') --di;
            *di++ = 'D';
            didexp = 1;
          }
        if(!didexp) {
          while(di > num && *(di - 1) == '0') --di;
          *di++ = 'D';
          *di++ = '0';
        }
        strcpy(p = di, buf);
      }
    }
    if(*s != '#') {
      for(si = di = s; ; ++si) {
        if(*si != '"') *di++ = *si;
        if(*si == 0) break;
/*        if(*si++ == '=') while(*(si + 1) == ' ') ++si; */
      }
      if(*s >= 'A') putchar('\t');
    }
    puts(s);
  } /* eof */
}

