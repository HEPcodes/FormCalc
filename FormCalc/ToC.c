/*
	ToC.c
		post-processes Mathematica's CForm output,
		removes " and continuation characters
		this file is part of FormCalc
		last modified 21 May 12 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
  char line[256], next[256], *s, *p;

  for( *next = 0;
       (*next) ? strcpy(line, next), *next = 0, (char *)1 :
                 fgets(line, sizeof line, stdin);
       puts(line) ) {
    char *eol = line + strlen(s = line) - 1;
    if( *eol != '\n' ) *++eol = '\n';
    if( eol[-1] == '\\' ) {
      if( fgets(next, sizeof next, stdin) == NULL ) break;
      else {
        char *p = next + strspn(next, " \t");
        int n = strcspn(p, " */()");
        memcpy(--eol, p, n);
        eol += n;
        memmove(p, p + n, strlen(p + n) + 1);
      }
    }
    *eol = 0;

    if( *s == '!' ) {
      *s++ = ' ';
      memmove(s, s - 1, eol - s + 2);
    }
    else if( *s != '#' && (p = strchr(s, '"')) ) {
      char *d = p;
      do {
        ++p;
        if( *p != '"' ) *d++ = *p;
      } while( *p );
    }
  } /* eof */

  return 0;
}

