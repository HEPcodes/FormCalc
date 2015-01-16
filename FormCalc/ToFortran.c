/*
	ToFortran.c
		post-processes Mathematica's FortranForm output
		- replaces all real constants by Fortran-style
		  double precision numbers (1.234D0),
		- removes " and indents lines not starting with #
		  (i.e. does not touch preprocessor statements)
		- replaces the continuation character (if any) by &
		this file is part of FormCalc
		last modified 14 Feb 12 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
  static const char signdigits[] = "+-0123456789";
  static const char *digits = signdigits + 2;
  char line[256], next[256], *s;

  for( *next = 0;
       (*next) ? strcpy(line, next), *next = 0, (char *)1 :
                 fgets(line, sizeof line, stdin);
       puts(s) ) {
    char *si, *di, *p;
    char *eol = line + strlen(s = line) - 1;
    if( *eol != '\n' ) *++eol = '\n';
    if( eol[-1] == '\\' ) {
      if( fgets(next, sizeof next, stdin) == NULL ) break;
      else {
        char *p = next + 6 + strspn(next + 6, " \t");
        int n = strcspn(p, " */()");
        memcpy(--eol, p, n);
        eol += n;
        memmove(p, p + n, strlen(p + n) + 1);
      }
    }
    *eol = 0;

    if( *s == '*' ) continue;

    si = s;
    while( (p = strpbrk(si, digits)) ) {
      char term;

      si = p + strspn(p, digits);
      if( p[-1] >= 'A' ) continue;  /* belongs to variable name */

      term = *si;
      if( term == '.' ) {
        if( *++si >= 'A' ) continue;
        si += strspn(si, digits);
      }

      if( (*si++ & 0xde) != 'D' ) {
        if( term != '.' ) continue;  /* is an integer */
        for( di = eol += 2; di > si; --di ) *di = di[-2];
        *si = '0';
      }
      si[-1] = 'D';
      si += strspn(si, signdigits);
    }

    if( *s == '!' ) *s = '\t';
    else if( *s != '#' ) {
      char *p = s;
      if( strncmp(p, "     ", 5) == 0 && *(p += 5) != ' ' ) *p++ = '&';
      if( (p = strchr(p, '"')) ) {
        char *d = p;
        do {
          ++p;
          if( *p != '"' ) *d++ = *p;
        } while( *p );
      }
      if( *s >= 'A' ) putchar('\t');
    }
  } /* eof */

  return 0;
}

