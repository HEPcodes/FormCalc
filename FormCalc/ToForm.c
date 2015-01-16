/*
	ToForm.c
		rearranges Mma's InputForm output to yield
		acceptable FORM input
		this file is part of FormCalc
		last modified 7 Apr 06 th
*/

#include <stdio.h>
#include <string.h>

int main()
{
  char in[2048], out[2048], *s, *d;

  while( !feof(stdin) ) {
    *in = 0;
    fgets(in, sizeof(in), stdin);

    if( *in == '#' ) {
      fputs(in, stdout);
      continue;
    }

    for( s = in, d = out; *s; ++s ) {
      if( *s == '\\' ) fgets(s, sizeof(in) - (int)(s - in), stdin);
      switch( *s ) {
      case '"':
        break;
      case '[':
        *d++ = '(';
        break;
      case ']':
        *d++ = ')';
        break;
      case ' ':
        if( *(s + 1) == '.' || *(s - 1) == '.' ) break;
        if( (int)(s - in) > 75 ) {
          *d++ = '\n';
          *d = 0;
          fputs(d = out, stdout);
          s = strcpy(in, s);
        }
        *d++ = *s;
        break;
      case '*':
      case '=':
        if( *(s - 1) == *s ) break;
      default:
        *d++ = *s;
      }
    }

    *d = 0;
    fputs(out, stdout);
  }
}

