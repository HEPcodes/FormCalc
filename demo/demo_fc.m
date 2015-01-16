<< ../FormCalc.m

$O1ME = Sqrt[S]

SetOptions[ProcessFile, CancelEps -> True]

counter = << "fa_output/zzzz.counter.amp"

ProcessFile["fa_output/zzzz.born.amp", "fc_output/born"]

ProcessFile[
  Join[ << "fa_output/zzzz.self.amp",
        Select[counter, DiagramType[#] == 2 &] ],
  "fc_output/self" ]

ProcessFile[
  Join[ << "fa_output/zzzz.vert.amp",
        Select[counter, DiagramType[#] == 1 &] ],
  "fc_output/vert" ]

ProcessFile[ "fa_output/zzzz.box.amp",
  "fc_output/box", ToughStuff -> True ]

Abbreviations[] >> "fc_output/abbr";
OptimizeAbbreviations["fc_output/abbr"]

