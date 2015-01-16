(*
	para.m
	  contains the constants and parameters of the Standard Model
	  to be used by NumPrep.m. Parameters that are not specified here
	  will appear as variables in the Fortran code.

Note: in spite of Mma's vigour to truncate digits to `allowable'
  precision it has proved valuable to enter the constants as rationals,
  i.e. numbers with infinite precision, even though this is nowhere
  justified on physical grounds.
*)

(* this will be the working precision;
   if you're pressed on time, NN=N will speed things up a little *)
NN = N[#, 35]&

GeV = 1;
MeV = 1/1000 GeV;
Alpha = 10^7/1370359895;
a2 = Alpha^2;
EL = Sqrt[4*Pi*Alpha]

(* boson masses *)
MW = 80401/1000 GeV;
MW2 = MW^2;
GammaW = 207/100 GeV;
MZ = 91188/1000 GeV;
MZ2 = MZ^2;
GammaZ = 2491/1000 GeV;
MLA = 0

(* fermion masses *)
fME = 51099906 10^-8 MeV;
fMD = 47 MeV;
fMU = 47 MeV;
fMM = 105658389 10^-6 MeV;
fMS = 150 MeV;
fMC = 155/100 GeV;
fML = 17771/10 MeV;
fMB = 45/10 GeV;
fMT = 1755/10 GeV

SW = Sqrt[(MZ2 - MW2)/MZ2];
CW = MW/MZ;
C4 = CW^-4;
C2 = CW^-2;
S4 = SW^-4;
S2 = SW^-2

BARN = 10^-28;
PICO = 10^-12;
HBARC2 = 38937966 10^-39;  (* \hbar c^2 *)
kinfactor = HBARC2/(PICO BARN)/(64 Pi^2 S)

