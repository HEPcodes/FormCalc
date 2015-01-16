(*
   IMPORTANT:
   Before using this program, you must adapt the NumPrep program to your
   Fortran compiler (NumPrep is located in the numerics/ subdirectory of
   the FormCalc tree). To this end, edit the NumPrep.m code, search for
   the string "compiler flags" and follow the instructions given in the
   comments there.
*)


<< ../../numerics/NumPrep.m

MakeCode[ "fc", "fortran",
  ProcessH -> "process.h_zzzz" ]

