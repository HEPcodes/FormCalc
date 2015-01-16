<< ../../FormCalc.m

$Scale = Sqrt[S]

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

	(* in ZZ -> ZZ, there are no box counter terms *)
ProcessFile[ "fa/box.amp", "fc/box", Classification -> Tough ]

Abbreviations[] >> fc/abbr;
OptimizeAbbreviations["fc/abbr"]

Print["time used: ", SessionTime[] - time1];

