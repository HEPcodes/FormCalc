<< NumPrep.m

(* IMPORTANT: you must choose here the correct compiler flags for two
   common options in Fortran:
   a) a flag to make all real constants double precision
      (i.e. interpret 1.234 as 1.234D0)
   b) a flag to ignore the 72-column restriction of Fortran-77.

   Unfortunately, the naming of these flags varies widely across different
   computer systems. Choices for some systems are listed below in the
   $F77 = ... line. If you don't find your system among them, please look
   up the appropriate options in your f77 man page and insert them there.

   If your f77 does not support these switches, you need to replace "f77"
   by "../f77c" below; f77c is a script which performs the same task as
   f77 except that it calls first f2c and then cc. (If you need f2c, go
   to ftp://ftp.netlib.org/f2c.)
*)


f77::undef = "Warning: I don't know the correct f77 compiler switches for
your system. I'll use f2c and proceed with fingers crossed. If I run into
problems, please update demo_num.m."

$F77 = Switch[ $SystemID,
  "DEC-AXP",	"f77 -r8 -extend_source -O",
  "HP-RISC",	"fort77 -R8 +es +U77 -O2",
  "Linux",	"../f77c -r8 -f -O",
  _,		Message[f77::undef]; "../f77c -r8 -f -O" ]

$LoopToolsDir = "$(HOME)/LoopTools/"	(* specify in gmake fashion *)

MakeCode[ "fc_output/", "fortran/", "*.m" ]

