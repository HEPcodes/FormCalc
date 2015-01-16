/*
	restack.c
		restacks data from num.F-output for easier access
		e.g. with gnuplot
		this file is part of FormCalc
		last modified 4 Jun 01 th

Syntax: restack var1 var2 var3... < infile > outfile

For example, if you have the following source file

# MH = 300.00000   MT = 175.00000
data1

# MH = 300.00000   MT = 150.00000
data2

# MH = 350.00000   MT = 175.00000
data3

# MH = 350.00000   MT = 150.00000
data4

the output from running "restack MH" will be

# MT = 175.00000
data1 300.00000
data3 350.00000

# MT = 150.00000
data2 300.00000
data4 350.00000

*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct dataline {
  struct dataline *next;
  char text[1024];
} DATALINE;

typedef struct commentline {
  struct commentline *next;
  struct dataline *lines;
  char varpick[512], varother[512];
} COMMENTLINE;


void putlines(COMMENTLINE *cp)
{
  DATALINE *dp, *dnext;

  for(dp = cp->lines; dp; dp = dnext) {
    printf("%s%s\n", dp->text, cp->varpick);
    dnext = dp->next;
    free(dp);
  }
}


void *memget(size_t amount)
{
  void *mem = malloc(amount);
  if(mem == NULL) {
    fputs("Out of memory. Exiting.\n", stderr);
    exit(2);
  }
  return mem;
}


main(int argc, char **argv)
{
  COMMENTLINE *cstart = NULL, *ccur, *clast;
  DATALINE *dcur, *dlast;
  char *d, *p, *v, s[1024];
  int c, found[25];

  if(argc < 2) {
    fprintf(stderr, "Usage: %s var1 var2 ... < infile > outfile\n"
      "  restacks data from num.F-output for easier\n"
      "  access e.g. with gnuplot.\n", argv[0]);
    exit(1);
  }

  memset(found, 0, sizeof(found));

  while(!feof(stdin)) {
    *s = 0;
    fgets(s, sizeof(s), stdin);
    if(*s == 0) break;
    *(s + strlen(s) - 1) = 0;
    v = s + strspn(s, " ");
    if(*v == 0) continue;
    d = s;
    do
      if(*v != ' ' || *(d - 1) != ' ') *d++ = *v;
    while(*v++);
    if(*s == '#') {
      if(!strchr(s, '=')) continue;
      *(cstart ? &clast->next : &cstart) =
        ccur = memget(sizeof(COMMENTLINE));
      clast = ccur;
      ccur->lines = NULL;
      ccur->next = NULL;
      d = ccur->varpick;
      for(c = 1; c < argc; ++c)
        if(p = strstr(s, argv[c])) {
          found[c] = 1;
          v = p + strlen(argv[c]);
          v += strspn(v, " =");
          *d++ = ' ';
          while(*v && *v != ' ') *d++ = *v++;
          strcpy(p, v);
        }
      *d = 0;
      strcpy(ccur->varother, s);
      dlast = NULL;
    }
    else {
      if(!cstart) {
        cstart = clast = ccur = memget(sizeof(COMMENTLINE));
        ccur->next = NULL;
        *ccur->varpick = *ccur->varother = 0;
        dlast = NULL;
      }
      *(dlast ? &dlast->next : &ccur->lines) =
        dcur = memget(sizeof(DATALINE));
      dlast = dcur;
      dcur->next = NULL;
      strcpy(dcur->text, s);
    }
  }

  while(cstart) {
    puts(cstart->varother);
    putlines(clast = cstart);
    for(ccur = cstart->next; ccur; ccur = clast->next) {
      if(strcmp(cstart->varother, ccur->varother)) clast = ccur;
      else {
        putlines(ccur);
        clast->next = ccur->next;
        free(ccur);
      }
    }
    clast = cstart->next;
    free(cstart);
    cstart = clast;
    putc('\n', stdout);
  }

  for(c = 1; c < argc; ++c)
    if(!found[c]) fprintf(stderr,
      "Warning: found no parameter \"%s\"\n", argv[c]);
}

