<< ../FormCalc.m

	(* use only the transverse part of the self-energies: *)
Pair[k[i_], k[j_]] = K2;
Pair[e[i_], k[j_]] = 0;
Pair[e[i_], e[j_]] = 1

$OnShell = False

MZ2/: CW^2 MZ2 = MW2;
MW2/: C2 MW2 = MZ2

ProcessFile["rc/self.aa", "rc/saa"];
ProcessFile["rc/self.az", "rc/saz"];
ProcessFile["rc/self.hh", "rc/shh"];
ProcessFile["rc/self.ww", "rc/sww"];
ProcessFile["rc/self.pp", "rc/spp"];
ProcessFile["rc/self.zz", "rc/szz"];
ProcessFile["rc/self.cc", "rc/scc"]

abr = Select[Abbreviations[], FreeQ[#, HoldForm]&]

insabr[file_] :=
Block[ {amp = Get[file]},
  Print[file];
  Put[amp /. abr, file] ]

insabr/@ FileNames["*.m", "rc"];

