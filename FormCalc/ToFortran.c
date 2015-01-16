/*
	ToFortran.c
		replaces all real constants by Fortran-style double
		precision numbers (1.234D0) in Mma FortranForm output
		this file is part of FormCalc
		last modified 6 Jun 05 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
  static const char signdigits[] = "+-0123456789";
  static const char *digits = signdigits + 2;

  char s[200], next[200], term, *eol, *pos, *si, *di;
  int gotnext = 0, n;

  while( !feof(stdin) ) {
    if( gotnext ) {
      strcpy(s, next);
      gotnext = 0;
    }
    else {
      *s = 0;
      fgets(s, sizeof(s), stdin);
    }
    eol = s + strlen(s) - 1;
    if( *eol != '\n' ) *++eol = '\n';
    if( *(eol - 1) == '\\' ) {
      fgets(next, sizeof(next), stdin);
      di = next + 6 + strspn(next + 6, " \t");
      n = strcspn(di, " */()");
      memcpy(--eol, di, n);
      eol += n;
      strcpy(di, di + n);
      gotnext = 1;
    }
    *eol = 0;

    si = s;
    while( (pos = strpbrk(si, digits)) ) {
      si = pos + strspn(pos, digits);
      if( *(pos - 1) >= 'A' ) continue;  /* belongs to variable name */

      term = *si;
      if( term == '.' ) si += strspn(++si, digits);

      if( (*si++ & 0xde) != 'D' ) {
        if( term != '.' ) continue;  /* is an integer */
        for( di = eol += 2; di > si; --di ) *di = *(di - 2);
        *si = '0';
      }
      *(si - 1) = 'D';
      si += strspn(si, signdigits);
    }

    if( *s != '#' ) {
      for( si = di = s; ; ++si ) {
        if( *si != '"' ) *di++ = *si;
        if( *si == 0 ) break;
      }
      if( *s >= 'A' ) putchar('\t');
    }
    puts(s);
  } /* eof */
}

