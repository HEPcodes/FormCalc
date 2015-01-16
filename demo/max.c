#include <stdlib.h>
#include <stdio.h>
#include <math.h>

main(int argc, char *argv[])
{
  FILE *f1;
  int i;
  double a, b, c, max = 0.;

  if(argc < 2) f1 = stdin;
  else {
    if((f1 = fopen(argv[1], "r")) == NULL) {
      fprintf(stderr, "%s not found\n", argv[1]);
      exit(1);
    }
  }
  for(i = 0; ; ++i) {
    if(fscanf(f1, "%lf %lf %lf", &a, &b, &c) == EOF) break;
    if(b > max) max = b;
    if(c > max) max = c;
  }
  fclose(f1);
/*
  c = .1 * pow(10., floor(log10(max)));
  b = ceil(max/c)*c;
  if((b - max)/b < .03) b += c;
*/
  i = (int)floor(log10(max));
  if(abs(i) <= 3) i = 0;
  printf("%d\n", i);
}

