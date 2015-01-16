/*
	ToForm.c
		rearranges Mma's InputForm output to yield
		acceptable FORM input
		this file is part of FormCalc
		last modified 20 Feb 01 th
*/

#include <stdio.h>

main()
{
  char line[2048], *s, *d;

  while(!feof(stdin)) {
    *line = 0;
    fgets(line, sizeof(line), stdin);
    if(*line != '#') {
      for(s = d = line; *s; ++s)
        switch(*s) {
        case '"':
          break;
        case '[':
          *d++ = '(';
          break;
        case ']':
          *d++ = ')';
          break;
        case ' ':
          if(*(s + 1) != '.' && *(s - 1) != '.') *d++ = *s;
          break;
        case '\\':
          if(*(s - 1) == '*') --d;
          break;
        case '*':
        case '=':
          if(*(s - 1) == *s) break;
        default:
          *d++ = *s;
        }
      *d = 0;
    }
    fputs(line, stdout);
  }
}

