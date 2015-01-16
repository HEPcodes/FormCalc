/*
	r8.c
		replaces all real constants by Fortran-style
		double precision numbers (1.234D0) in Mma
		FortranForm output
		last modified 16 Mar 99 th
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
      gets(s);
    }
    p = s + strlen(s) - 1;
    if(*p == '\\') {
      gets(next);
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
      for(si = di = s; ; ) {
        if(*si != '"') *di++ = *si;
        if(*si == 0) break;
        if(*si++ == '=') while(*(si + 1) == ' ') ++si;
      }
    }
    puts(s);
  } /* eof */
}

