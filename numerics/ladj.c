/*
   ladj.c     2aug96 th
   adjusts the legend in gnuplot's pslatex output so that
   the labels don't stick together.

   usage: ladj [spacing] < infile > outfile

   where <spacing> is optional and defaults to 60

   [12aug97 th] patched up rel labels + gpl-3.6 adaptation
   [25sep97 th] patched for multiple legends
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

typedef struct data {
  struct data *next;
  char s[256];
} DATA;

int main(int argc, char **argv)
{
  char t[255], *pos, *p2;
  int addto = 0, distance = 60;
  int i, x, y, line = 0, lasty, lastx = 0, dx, dy;
  DATA *start, *da, *da2, *r, *on = NULL;

  if(argc > 1) distance = atoi(argv[1]);

  da = da2 = start = malloc(sizeof(DATA));
  while(!feof(stdin)) {
    da2->next = da = malloc(sizeof(DATA));
    *(da->s) = 0;
    gets(da->s);
    if(strncmp(da->s, "\\put(", 5) == 0 &&
       (pos = strstr(da->s, "\\makebox(0,0)[")) &&
       (*(pos += 14) == 'r' || *pos == 'l')) {
      sscanf(da->s + 5, "%d,%d", &x, &y);
      if((i = (x == lastx && y == lasty + 100)) && !on)
        on = da2, line = 0;
      lastx = x;
      lasty = y;
    }
    else i = 0;
    ++line;
    if(on && !i) {
      addto = distance*line;
      for(da2 = on; da2 != da; addto -= distance, da2 = da2->next) {
        x = strtod(strcpy(t, da2->s) + 5, &pos);
        y = strtod(++pos, &p2);
        *pos = 0;
        sprintf(da2->s, "%s%d%s", t, y - addto, p2);
        sprintf(t, "%d M", y);
        i = strlen(t);
        for(r = start->next; r != on; r = r->next)
          if(strncmp(t, pos = r->s + strlen(r->s) - i, i) == 0) {
            if(abs(x - atoi(r->s)) > 350) continue;
            sprintf(pos, "%d M", y - addto);
            pos = (r = r->next)->next->s;
            if(*(pos + strlen(pos) - 1) == 'R') {
              sscanf(pos, "%d %d", &dx, &dy);
              sprintf(pos, "%d %d R", dx, dy + addto);
            }
            break;
          }
      }
      on = NULL;
    }
    da2 = da;
  }

  while(start != da) {
    puts((da2 = start->next)->s);
    free(start);
    start = da2;
  }
}

