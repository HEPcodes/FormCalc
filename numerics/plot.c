/*
	plot.c
		make plots of files produced by num.F
		last modified 6 Jul 99 th

This program works in two steps:

First, it generates a script, called "file.gpl" if "file" was the original
file name, which contains an invocation of gnuplot with a certain set of
parameter settings. These settings depend on the name "file", since num.F
encodes the parameters of a calculation in the file name, e.g.
run-diff.pol=UUUU.E=00500, which means: this file contains a differential
cross-section at a CMS energy of 500 GeV, all particles unpolarized.

Second, it executes the newly created .gpl script to produce the actual
EPS file. This is again not done directly from gnuplot. Instead, "pslatex"
mode is used which produces a LaTeX file in which the graphics part is
embedded in a PostScript \special. The advantage is that all labels are
genuine LaTeX labels. This file then is run through LaTeX and dvips -E to
produce the final .eps figure.

Clearly, this is quite some effort, and only justified with the hindsight
that after this program has terminated the .gpl file can be edited for
fine tuning. When editing the .gpl script, note that TeX metacharacters
are replaced as & -> \ and ? -> $ to avoid conflict with gnuplot 3.6.

Options:
  -born		use only cols 1:2 instead of 1:2,1:3
  -o outname	write output to outname
  -tot		use conventions for total cs plot
  -diff		use conventions for diff cs plot
  -lny		use log y scale regardless of diff/tot cs
  -nolny	same for linear scale

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>


typedef struct {
  enum PLOTTYPE {NONE, DIFF, TOT} type;
  int Efrom, Eto;
  char pol[40], file[128];
} PARA;


void analyse_filename(char *s, PARA *d)
{
  char pol[20], *p, *dpol;

  strcpy(d->file, s);
  if(p = strstr(s, "diff.")) {
    d->type = DIFF;
    sscanf(p + 5, "pol=%4s.E=%5d", pol, &d->Efrom);
  }
  else if(p = strstr(s, "tot.")) {
    d->type = TOT;
    sscanf(p + 4, "pol=%4s.E=%5d-%5d", pol, &d->Efrom, &d->Eto);
  }
  else {
    d->type = NONE;
    *d->pol = 0;
    return;
  }
  for(p = pol, dpol = d->pol; p < &pol[4]; ++p)
    if(*p == '-') {
      *dpol++ = '?';
      *dpol++ = *p;
      *dpol++ = '?';
    }
    else *dpol++ = *p;
  *dpol = 0;
}


static char *unit(char *var, int num)
{
  static char s[40];

  if(num % 1000 == 0)
    sprintf(s, "?%s = %d?~TeV, ", var, num/1000);
  else sprintf(s, "?%s = %d?~GeV, ", var, num);
  return s;
}


void max(char *file, double *max)
{
  FILE *tmp;
  char s[200];
  double born, loop;
  int m;

  tmp = fopen(file, "r");
  if(tmp == NULL) {
    fprintf(stderr, "%s doesn't exist.\n", file);
    exit(1);
  }
  while(fgets(s, sizeof(s), tmp)) {
    if(sscanf(s, "%*lf %lf %lf", &born, &loop) == 2) {
      if(born > *max) *max = born;
      if(loop > *max) *max = loop;
    }
  }
  fclose(tmp);
}


main(int argc, char **argv)
{
  PARA filez[20], *fp = filez, *fp2;
  int flag_lny = -1, flag_o = 0, flag_born = 0;
  enum PLOTTYPE thetype = NONE;
  char **argp = argv, outname[128], title[256];
  char theunit[30], thru[30];
  FILE *out;
  double ymax;
  int Efrom, Eto, m, difEfrom, start;
  char *difpol;

  *outname = 0;
  while(--argc) {
    if(**++argp == '-') {
      if(strcmp(*argp, "-lny") == 0) flag_lny = 1;
      else if(strcmp(*argp, "-nolny") == 0) flag_lny = 0;
      else if(strcmp(*argp, "-o") == 0) flag_o = 1;
      else if(strcmp(*argp, "-born") == 0) flag_born = 1;
      else if(strcmp(*argp, "-diff") == 0) thetype = DIFF;
      else if(strcmp(*argp, "-tot") == 0) thetype = TOT;
      else fprintf(stderr, "ignored: %s\n", *argp);
    }
    else {
      if(flag_o) {
        strcpy(outname, *argp);
        flag_o = 0;
      }
      else {
        if(strstr(*argp, ".eps") || strstr(*argp, ".gpl"))
          fprintf(stderr, "ignored: %s\n", *argp);
        else analyse_filename(*argp, fp++);
      }
    }
  }
  if(fp == filez) {
    fprintf(stderr, "usage: %s [-lny] [-nolny] [-tot] [-diff] "
      "[-born] [-o outputfile] files...\n", *argv);
    exit(1);
  }

  if(thetype == NONE) {
    thetype = filez->type;
    if(flag_lny == -1) flag_lny = thetype == TOT;
  }

  if(*outname == 0) {
    strcpy(outname, filez->file);
    strcat(outname, ".gpl");
  }
  out = fopen(outname, "w");
  if(out == NULL) {
    fprintf(stderr, "cannot open %s\n", outname);
    exit(1);
  }

  difEfrom = -1;
  difpol = NULL;
  start = 1;
  for(fp2 = filez; fp2 < fp; ++fp2)
    if(fp2->type != NONE) {
      if(start) {
        if(fp2->type == DIFF) difEfrom = fp2->Efrom;
        difpol = fp2->pol;
        start = 0;
      }
      else {
        if(fp2->type == DIFF && difEfrom != fp2->Efrom) difEfrom = -1;
        if(difpol && strcmp(difpol, fp2->pol)) difpol = NULL;
      }
    }
  *title = 0;
  if(!start) {
    if(difEfrom > 0) strcat(title, unit("&sqrt{s}", difEfrom));
    if(difpol) strcat(title, difpol);
    else if((m = strlen(title)) > 0) *(title + m - 2) = 0;
  }
  fprintf(out,
    "#!/bin/sh\n\n"
    "# This plotting script was generated automatically by plot.\n"
    "# If you change or add labels, remember to replace \\ by & and $ by ?\n"
    "# in LaTeX text to avoid conflict with gnuplot's metacharacters.\n\n"
    "trap \"rm -f v$$.*\" 0 1 2 3 9 15\n\n"
    "gnuplot << _EOF_\n\n"
    "# ----- The gnuplot commands start here -----\n\n"
    "set term pslatex norotate\n"
    "set output \"v$$.plot\"\n"
    "set key\n"
    "set title \"%s\" 0,-.5\n"
    "set lmargin 11\n"
    "set rmargin 15\n",
    title);

  *thru = 0;
  strcpy(theunit, "&mathrm{pb}");
  if(!flag_lny) {
    ymax = 0.;
    for(fp2 = filez; fp2 < fp; ++fp2)
      if(fp2->type != NONE) max(fp2->file, &ymax);
    if(ymax != 0. && abs(m = (int)floor(log10(ymax))) > 3) {
      sprintf(theunit, "10^{%d}&,&mathrm{pb}", m);
      sprintf(thru, " thru (x*1e%d)", -m);
    }
  }

  if(thetype == TOT) {
    Efrom = 100;
    Eto = 10000;
    for(fp2 = filez; fp2 < fp; ++fp2)
      if(fp2->type == TOT) {
        Efrom = fp2->Efrom;
        Eto = fp2->Eto;
        break;
      }
    fprintf(out,
      "set size 1,1.3\n"
      "set logscale x%s\n"
      "set xtics (200,500,1000,2000,5000,10000,20000,50000)\n"
      "set xrange [%d:%d]\n"
      "set xlabel \"?&sqrt s/?GeV\" 32,1.5\n"
      "set ylabel \"?&dfrac{&stot}{%s}?\" 0,14\n",
      flag_lny ? "y" : "", Efrom, Eto, theunit);
  }
  else {
    if(flag_lny) fprintf(out, "set logscale y\n");
    fprintf(out,
      "set size .8,1\n"
      "set xtics (\"&small ?0^&circ?\" 0, \"\" pi/4, \\\n"
      "  \"&small ?90^&circ?\" pi/2, \"\" 3*pi/4, \\\n"
      "  \"&small ?180^&circ?\" pi)\n"
      "set xrange [0:pi]\n"
      "set xlabel \"?&theta?\" 21,1.5\n"
      "set ylabel \"?&dfrac{&dsdO}{%s}?\" 0,9\n",
      theunit);
  }
  fprintf(out,
    "set format y \"?%s?\"\n"
    "plot \\\n",
    flag_lny ? "10^{%L}" : "%g");

  for(start = 1, fp2 = filez; ; ++start) {
    printf("> %s\n", fp2->file);
    *title = 0;
    if(fp2->type == DIFF && difEfrom < 0) 
      strcat(title, unit("&sqrt{s}", fp2->Efrom));
    if(!difpol) strcat(title, fp2->pol);
    else if((m = strlen(title)) > 0) *(title + m - 2) = 0;

    if(flag_born)
      fprintf(out,
        "  \"%s\"%s u 1:2 \\\n"
        "    t \"&small %s\" w l %d",
        fp2->file, thru, title, start);
    else {
      if(*title) strcat(title, ", ");
      fprintf(out,
        "  \"%s\"%s u 1:2 \\\n"
        "    t \"&small %sBorn\" w l %d,\\\n"
        "  \"\"%s u 1:3 \\\n"
        "    t \"&small full\" w l %d",
        fp2->file, thru, title, start + 3, thru, start);
    }
    if(++fp2 == fp) break;
    fprintf(out, ",\\\n");
  }

  fprintf(out,
    "\n\n"
    "# ----- The gnuplot commands end here -----\n\n"
    "_EOF_\n\n"
    "cat << _EOF_ > v$$.tex\n"
    "\\\\documentclass[11pt]{article}\n"
    "\\\\oddsidemargin=0pt\n"
    "\\\\evensidemargin=0pt\n"
    "\\\\parindent=0pt\n"
    "\\\\pagestyle{empty}\n"
    "\\\\def\\\\dfrac#1#2{{\\\\displaystyle{#1\\\\over #2}}}\n"
    "\\\\def\\\\d{\\\\mathrm{d}}\n"
    "\\\\def\\\\dsdO{\\\\d\\\\sigma/\\\\d\\\\Omega}\n"
    "\\\\def\\\\stot{\\\\sigma_{\\\\mathrm{tot}}}\n\n"
    "%% it's a good idea to use PostScript fonts here since the figures\n"
    "%% will likely be expanded or shrunk to fit the final size:\n\n"
    "\\\\renewcommand{\\\\rmdefault}{ppl}\n"
    "\\\\DeclareSymbolFont{operators}{OT1}{pplcm}{m}{n}\n"
    "\\\\DeclareSymbolFont{letters}{OML}{pplcm}{m}{it}\n"
    "\\\\DeclareSymbolFont{largesymbols}{OMX}{psycm}{m}{n}\n"
    "\\\\DeclareSymbolFont{bold}{OT1}{ppl}{bx}{n}\n"
    "\\\\DeclareSymbolFont{italic}{OT1}{ppl}{m}{it}\n"
    "\\\\DeclareMathAlphabet{\\\\mathrm}{OT1}{ppl}{m}{n}\n"
    "\\\\DeclareMathAlphabet{\\\\mathbf}{OT1}{ppl}{bx}{n}\n"
    "\\\\DeclareMathAlphabet{\\\\mathit}{OT1}{ppl}{m}{it}\n\n"
    "\\\\begin{document}\n"
    "_EOF_\n\n"
    "sed -e 's/&/\\\\/g' -e 's/?/$/g' -e '/endinput/d' v$$.plot | \\\n"
    "  ladj_$HOSTTYPE 40 >> v$$.tex\n\n"
    "echo '\\end{document}' >> v$$.tex\n\n"
    "epsfile=\"`basename $0 .gpl`.eps\"\n"
    "latex v$$.tex\n"
    "dvips -E -o $epsfile v$$.dvi\n"
    "gv -magstep 3 $epsfile &\n\n");

  fchmod(fileno(out), 0755);
  fclose(out);

  system(outname);
}

