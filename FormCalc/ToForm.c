/*
	ToForm.c
		rearranges Mathematica's InputForm output to
		yield acceptable FORM input
		this file is part of FormCalc
		last modified 28 Oct 16 th
*/

#include <stdio.h>
#include <string.h>

int main() {
  char in[1024], out[1024];
  int q = 1, m = -1;

  while( fgets(in, sizeof in, stdin) ) {
    char *s, *d, closing = ')';

    if( in[0] == '#' ) {
      if( in[1] == '#' ) {
        while( fgets(in, sizeof in, stdin) ) {
          if( in[0] == '#' && in[1] == '#' ) break;
          fputs(in, stdout);
        }
        continue;
      }
      m = 0xdf;
    }

    for( s = in, d = out; *s; ++s ) {
more:
      switch( *s & m ) {
      case '\\':
        if( s[1] == '\n' ) {
          fgets(s, in + sizeof in - s, stdin);
          goto more;
        }
        break;
      case '"':
        break;
      case '[':
        *d++ = (s[-1] == '\\') ? (closing = ']', '[') : '(';
        break;
      case ']':
        *d++ = closing;
        closing = ')';
        break;
      case '{':
        *d++ = 'L';
        *d++ = 'i';
        *d++ = 's';
        *d++ = 't';
        *d++ = '(';
        break;
      case '}':
        *d++ = ')';
        break;
      case '$':
        *d++ = (s[-1] > '@' || (unsigned)(s[-1] - '0') < 10) ? 'S' : '$';
        break;
      case ' ':
        if( s[1] == '.' || s[-1] == '.' ) break;
        if( s - in > 75 ) {
          *d++ = '\n';
          *d = 0;
          fputs(d = out, stdout);
          s = strcpy(in, s);
        }
        *d++ = *s;
        break;
      case '*':
      case '=':
        if( s[-1] == *s ) break;
      default:
        q ^= (*s == '"');
        *d++ = *s;
      }
    }

	/* line break not allowed inside a dot product (rare): */
    if( d > out + 2 && d[-2] == '.' ) {
      fgets(in, sizeof in, stdin);
      s = in + strspn(in, " \t");
      --d;
      goto more;
    }

    *d = 0;
    fputs(out, stdout);

    m |= -q;
  }

  return 0;
}

