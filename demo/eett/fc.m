<< ../../FormCalc.m

	(* if the calculation takes too long on your computer,
	   comment out the next line *)
o1 = o2 = Simplify

Small[ME] = Small[ME2] = 0;

time1 = SessionTime[]

counter = << fa/counter.amp

ProcessFile["fa/born.amp", "fc/born"]

ProcessFile[
  Join[ << fa/self.amp,
        Select[counter, DiagramType[#] == 2 &] ],
  "fc/self" ]

ProcessFile[
  Join[ << fa/vert.amp,
        Select[counter, DiagramType[#] == 1 &] ],
  "fc/vert" ]

	(* in ee -> tt, there are no box counter terms *)
ProcessFile["fa/box.amp", "fc/box"]

	(* calculate unpolarized matrix elements only *)
Hel[_] = 0;
HelicityME[All, << fc/bornF.m] >> fc/mat;
Abbreviations[] >> fc/abbr

Print["time used: ", SessionTime[] - time1];

