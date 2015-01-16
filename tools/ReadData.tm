:Begin:
:Function:	readdata
:Pattern:	ReadData[filename_String, setn_Integer:1, para_:Para, data_:Data]
:Arguments:	{filename, setn, ToString[para], ToString[data]}
:ArgumentTypes:	{String, Integer, String, String}
:ReturnType:	Manual
:End:

:Evaluate:	DataRow = List
:Evaluate:	Attributes[Set0] = {HoldFirst}
:Evaluate:	Set0[x__] := (Set[x]; 0)
:Evaluate:	ReadData::noopen = "Cannot open `1`."

/*
	ReadData.tm
		reads data files produced by num.F into Mathematica
		this file is part of FormCalc
		last modified 16 Oct 01 th
*/

#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char copyleft[] =
  "@(#) ReadData utility for FormCalc, 16 Oct 01 Thomas Hahn";

#ifndef MLCONST
#define MLCONST
#endif


typedef struct {
  double val;
  char name[20];
} VAR;

VAR vars[20];
int varcount = 0;


void parse_para(char *s)
{
  int n;
  char *s2;

  for(s2 = s; *s2; ++s2)
    switch(*s2) {
    case '_':
      *s2 = '$';
      break;
    case '(':
      *s2 = '[';
      break;
    case ']':
      *s2 = ']';
      break;
    }

  do {
    sscanf(s, " %[^=]=%lg %n",
      vars[varcount].name, &vars[varcount].val, &n);
    ++varcount;
  } while(*(s += n));
}


void xmit_para(const char *head, int setn)
{
  int p, i;

  MLPutFunction(stdlink, "EvaluatePacket", 1);

  MLPutFunction(stdlink, "Set0", 2);
  MLPutFunction(stdlink, head, 1);
  MLPutInteger(stdlink, setn);
  MLPutFunction(stdlink, "List", varcount);
  for(p = 0; p < varcount; ++p) {
    MLPutFunction(stdlink, "Rule", 2);
    MLPutFunction(stdlink, "ToExpression", 1);
    MLPutString(stdlink, vars[p].name);
    MLPutDouble(stdlink, vars[p].val);
  }
  varcount = 0;

  MLEndPacket(stdlink);

  while((p = MLNextPacket(stdlink)) && p != RETURNPKT)
    MLNewPacket(stdlink);
  MLNewPacket(stdlink);
}


void xmit_data(char **begin, char ***end, const char *head, int setn)
{
  double dd[8];
  int i, c, n, p;

  n = (int)(*end - begin);
  *end = begin;

  MLPutFunction(stdlink, "EvaluatePacket", 1);

  MLPutFunction(stdlink, "Set0", 2);
  MLPutFunction(stdlink, head, 1);
  MLPutInteger(stdlink, setn);
  MLPutFunction(stdlink, "List", n);
  while(n--) {
    c = sscanf(*begin, "%lg %lg %lg %lg %lg %lg %lg %lg",
      &dd[0], &dd[1], &dd[2], &dd[3], &dd[4], &dd[5], &dd[6], &dd[7]);
    free(*begin++);
    MLPutFunction(stdlink, "DataRow", c);
    for(i = 0; i < c; ++i) MLPutDouble(stdlink, dd[i]);
  }

  MLEndPacket(stdlink);

  while((p = MLNextPacket(stdlink)) && p != RETURNPKT)
    MLNewPacket(stdlink);
  MLNewPacket(stdlink);
}


void readdata(const char *filename, int setn,
  const char *para, const char *data)
{
  FILE *file;
  char line[512], *s, *dataline[2048], **stor = dataline;
  int xmit = 0, p;

  file = fopen(filename, "r");
  if(file == NULL) {
    MLPutFunction(stdlink, "CompoundExpression", 2);
    MLPutFunction(stdlink, "Message", 2);
    MLPutFunction(stdlink, "MessageName", 2);
    MLPutSymbol(stdlink, "ReadData");
    MLPutString(stdlink, "noopen");
    MLPutString(stdlink, filename);
    MLPutSymbol(stdlink, "$Failed");
    return;
  }

  while(!feof(file)) {
    *line = 0;
    fgets(line, sizeof(line), file);
    s = line + strspn(line, " \t");
    if(*s >= '+' && *s <= '9') {
      if(varcount) xmit_para(para, setn);
      xmit = 1;
      *stor++ = strdup(s);
    }
    else {
      if(strchr(s, '=')) parse_para(s + strspn(s, "# \t"));
      if(xmit) {
        xmit = 0;
        xmit_data(dataline, &stor, data, setn++);
      }
    }
  }

  if(xmit) xmit_data(dataline, &stor, data, setn++);

  fclose(file);

  MLPutInteger(stdlink, setn - 1);
}


main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

