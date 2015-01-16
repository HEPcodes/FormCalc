(*
	para.m
		contains the constants and parameters of the Standard
		Model to be used by NumPrep.m
		last modified 2 Mar 99 th

Notes:

- Parameters that are not specified here will appear as variables in the
  Fortran code.

- In spite of Mma's vigour to truncate digits to `allowable' precision it
  has proved valuable to enter the constants as rationals, i.e. numbers
  with infinite precision, even though this is nowhere justified on
  physical grounds.

- IMPORTANT: if you change things here, be sure to re-run everything that
  depends on it -- not only the main code but also e.g. the
  renormalization constants.

*)


GeV = 1;
MeV = 1/1000 GeV;
Alpha = 10^7/1370359895;
Gmue = 116639 10^-10 GeV^-2;
a2 = Alpha^2;
EL = Sqrt[4 Pi Alpha]

(* boson masses *)
MW2 = (MW = 80356/1000 GeV)^2;
GammaW = 207/100 GeV;
MZ2 = (MZ = 911863/10000 GeV)^2;
GammaZ = 2491/1000 GeV;
MLA = mla = 0

(* fermion masses *)
ME2 = (ME = 51099907 10^-8 MeV)^2;
MU2 = (MU = 464/10 MeV)^2;
MD2 = (MD = 465/10 MeV)^2;
MM2 = (MM = 105658389 10^-6 MeV)^2;
MC2 = (MC = 150/100 GeV)^2;
MS2 = (MS = 150 MeV)^2;
ML2 = (ML = 1777 MeV)^2;
(* MT2 = (MT = 175 GeV)^2; *)
MB2 = (MB = 45/10 GeV)^2

SW = Sqrt[(MZ2 - MW2)/MZ2];
CW = MW/MZ;
C4 = CW^-4;
C2 = CW^-2;
S4 = SW^-4;
S2 = SW^-2

