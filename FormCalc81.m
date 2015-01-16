(*

This is FormCalc, Version 8.1
Copyright by Thomas Hahn 1996-2013
last modified 27 Mar 13 by Thomas Hahn

Release notes:

FormCalc is free software, but is not in the public domain.
Instead it is covered by the GNU Lesser General Public License.
In plain English this means:

1. We don't promise that this software works.
   (But if you find any bugs, please let us know!)

2. You can use this software for whatever you want.
   You don't have to pay us.

3. You may not pretend that you wrote this software.
   If you use it in a program, you must acknowledge
   somewhere in your publication that you've used
   our code.

If you're a lawyer, you can find the legal stuff at
http://gnu.org/licenses/lgpl.html.

The user guide for this program can be found at
http://feynarts.de/formcalc.

If you find any bugs, or want to make suggestions, or
just write fan mail, address it to:
	Thomas Hahn
	Max Planck Institute for Physics
	Foehringer Ring 6
	D-80805 Munich, Germany
	e-mail: hahn@feynarts.de

Have fun!

*)

Print[""];
Print["FormCalc 8.1"];
Print["by Thomas Hahn"];
Print["last revised 27 Mar 13"]


(* symbols from FeynArts *)

BeginPackage["FeynArts`"]

{ FeynAmp, FeynAmpList, Process, GraphID,
  Generic, Classes, Particles, S, F, U, V, T,
  Insertions, G, Mass, GaugeXi, VertexFunction,
  PropagatorDenominator, FeynAmpDenominator,
  FourMomentum, Internal, External, TheMass,
  Index, IndexDelta, IndexEps, IndexSum, SumOver,
  MatrixTrace, FermionChain, NonCommutative,
  CreateTopologies, ExcludeTopologies, Tadpoles,
  InitializeModel, $Model, Model, GenericModel, Reinitialize,
  InsertFields, InsertionLevel, ExcludeParticles,
  ExcludeFieldPoints, LastSelections, Restrictions,
  CreateFeynAmp, Truncated, RenConst, ProcessName,
  Paint, DiagramGrouping }

EndPackage[]


(* symbols from LoopTools *)

BeginPackage["LoopTools`"]

PaVeIntegral = A0i | B0i | C0i | D0i | E0i | F0i |
  A0 | A00 | B0 | B1 | B00 | B11 | B001 | B111 |
  DB0 | DB1 | DB00 | DB11

{ Nbb, bb0, bb1, bb00, bb11, bb001, bb111, dbb0, dbb1, dbb00, dbb11,
  Ncc, cc0, cc00, cc001, cc002, cc0000, cc0011, cc0022, cc0012,
  Ndd, dd0, dd0000, dd00001, dd00002, dd00003,
  Nee, ee0,
  Nff, ff0 }

A0::usage =
"A0[m] is the one-point one-loop scalar integral.  m is the mass
squared."

A00::usage =
"A00[m] is the one-point tensor coefficient of g_{mu nu}.  m is the mass
squared."

B0i::usage =
"B0i[id, p, m1, m2] is the generic two-point one-loop integral which
includes scalar and tensor coefficients as well as their derivatives
with respect to p, specified by id.  For example, B0i[bb0, ...] is the
scalar function B_0, B0i[bb11, ...] the tensor coefficient function
B_{11} etc.  p is the external momentum squared and m1 and m2 are the
masses squared."

Bput::usage =
"Bput[p, m1, m2] computes all two-point coefficients in LoopTools."

C0i::usage =
"C0i[id, p1, p2, p1p2, m1, m2, m3] is the generic three-point one-loop
integral which includes both scalar and tensor coefficients, specified
by id.  For example, C0i[cc0, ...] is the scalar function C_0,
C0i[cc112, ...] the tensor coefficient function C_{112} etc.  p1, p2,
and p1p2 are the external momenta squared and m1, m2, m3 are the masses
squared."

Cput::usage =
"Cput[p1, p2, p1p2, m1, m2, m3] computes all three-point coefficients in
LoopTools."

D0i::usage =
"D0i[id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the generic
four-point one-loop integral which includes both scalar and tensor
coefficients, specified by id.  For example, D0i[dd0, ...] is the scalar
function D_0, D0i[dd1233, ...] the tensor function D_{1233} etc. 
p1...p4 are the external momenta squared, p1p2 and p2p3 are the squares
of external momenta (1+2) and (2+3), respectively, and m1...m4 are the
masses squared."

Dput::usage =
"Dput[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] computes all
four-point coefficients in LoopTools."

E0i::usage =
"E0i[id, p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3,
m4, m5] is the generic five-point one-loop integral which includes both
scalar and tensor coefficients, specified by id.  For example,
E0i[ee0, ...] is the scalar function E_0, E0i[ee1244, ...] the tensor
function E_{1244} etc.  p1...p5 are the external momenta squared,
p1p2...p5p1 are the squares of external momenta (1+2)...(5+1),
respectively, and m1...m5 are the masses squared."

Eput::usage =
"Eput[p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4,
m5] computes all five-point coefficients in LoopTools."

F0i::usage =
"F0i[id, p1, p2, p3, p4, p5, p6, p1p2, p2p3, p3p4, p4p5, p5p6, p6p1,
p1p2p3, p2p3p4, p3p4p5, m1, m2, m3, m4, m5, m6] is the generic six-point
one-loop integral which includes both scalar and tensor coefficients,
specified by id.  For example, F0i[ff0, ...] is the scalar function F_0,
F0i[ff1244, ...] the tensor function F_{1244} etc.  p1...p6 are the
external momenta squared, p1p2...p6p1 are the squares of external
momenta (1+2)...(6+1), respectively, p1p2p3...p3p4p5 are the external
momenta (1+2+3)...(3+4+5) squared, and m1...m6 are the masses squared."

Fput::usage =
"Fput[p1, p2, p3, p4, p5, p6, p1p2, p2p3, p3p4, p4p5, p5p6, p6p1,
p1p2p3, p2p3p4, p3p4p5, m1, m2, m3, m4, m5, m6] computes all six-point
coefficients in a future version of LoopTools."


(* compatibility functions *)

B0::usage =
"B0[p, m1, m2] is the two-point one-loop scalar integral."

B1::usage =
"B1[p, m1, m2] is the coefficient of k_mu in the two-point one-loop
tensor integral B_mu."

B00::usage =
"B00[p, m1, m2] is the coefficient of g_{mu nu} in the two-point
one-loop tensor integral B_{mu nu}."

B11::usage =
"B11[p, m1, m2] is the coefficient of k_mu k_nu in the two-point
one-loop tensor integral B_{mu nu}."

B001::usage =
"B001[p, m1, m2] is the coefficient of g_{mu nu} k_rho in the two-point
one-loop tensor integral B_{mu nu rho}."

B111::usage =
"B111[p, m1, m2] is the coefficient of k_mu k_nu k_rho in the two-point
one-loop tensor integral B_{mu nu rho}."

DB0::usage =
"DB0[p, m1, m2] is the derivative of B0[p, m1, m2] with respect to p."

DB1::usage =
"DB1[p, m1, m2] is the derivative of B1[p, m1, m2] with respect to p."

DB00::usage =
"DB00[p, m1, m2] is the derivative of B00[p, m1, m2] with respect to p."

DB11::usage =
"DB11[p, m1, m2] is the derivative of B11[p, m1, m2] with respect to p."

C0::usage =
"C0[p1, p2, p1p2, m1, m2, m3] is the three-point scalar one-loop
integral."

D0::usage =
"D0[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the four-point scalar
one-loop integral."

E0::usage =
"E0[p1, p2, p3, p4, p5, p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4,
m5] is the five-point scalar one-loop integral."

F0::usage =
"F0[p1, p2, p3, p4, p5, p6, p1p2, p2p3, p3p4, p4p5, p5p6, p6p1, p1p2p3,
p2p3p4, p3p4p5, m1, m2, m3, m4, m5, m6] is the six-point scalar one-loop
integral."

PaVe::usage =
"PaVe[ind, {pi}, {mi}] is the generalized Passarino-Veltman function
used by FeynCalc.  It is converted to B0i, C0i, D0i, E0i, or F0i in
FormCalc."

ToOldBRules::usage =
"ToOldBRules is a list of rules for converting two-point functions to
the old (LoopTools 2.1) conventions."

ToNewBRules::usage =
"ToNewBRules is a list of rules for converting two-point functions to
the new (LoopTools 2.2) conventions."

EndPackage[]


(* OPP symbols *)

BeginPackage["OPP`"]

CutIntegral = Acut | Bcut | Ccut | Dcut | Ecut | Fcut

CutMasters = Amas | Bmas | Cmas | Dmas | Emas | Fmas

{Mbb, Mcc, Mdd, Mee, Mff}

Bcut::usage =
"Bcut[b0args][hel, r, num4, numE, k, m1, m2] is the two-point OPP
integral with numerator functions num4 (4-dim) and numE (D-4-dim), where
num4 contains r powers of the integration momentum.  For vanishing hel
the coefficient is known to be zero so the integral need not be
evaluated.  b0args are the same arguments as for B0, k is the external
momentum and m1, m2 are the masses squared."

Ccut::usage =
"Ccut[c0args][rank, num, numtilde, k1, k2, m1, m2, m3] is the
three-point OPP integral with numerator functions num4 (4-dim) and numE
(D-4-dim), where num4 contains r powers of the integration momentum.
For vanishing hel the coefficient is known to be zero so the integral
need not be evaluated.  c0args are the same arguments as for C0,
k1, k2 are the external momenta and m1..m3 are the masses squared."

Dcut::usage =
"Dcut[d0args][rank, num, numtilde, k1, k2, k3, m1, m2, m3, m4] is the 
four-point OPP integral with numerator functions num4 (4-dim) and numE 
(D-4-dim), where num4 contains r powers of the integration momentum.
For vanishing hel the coefficient is known to be zero so the integral 
need not be evaluated.  d0args are the same arguments as for D0,
k1..k3 are the external momenta and m1..m4 are the masses squared."

Ecut::usage =
"Ecut[e0args][rank, num, numtilde, k1, k2, k3, k4, m1, m2, m3, m4, m5]
is the five-point OPP integral with numerator functions num4 (4-dim)
and numE (D-4-dim), where num4 contains r powers of the integration
momentum.  For vanishing hel the coefficient is known to be zero so
the integral need not be evaluated.  e0args are the same arguments as
for E0, k1..k4 are the external momenta and m1..m5 are the masses
squared."

Fcut::usage =
"Fcut[f0args][rank, num, numtilde, k1, k2, k3, k4, k5, m1, m2, m3, m4,
m5, m6] is the six-point OPP integral with numerator functions num4
(4-dim) and numE (D-4-dim), where num4 contains r powers of the
integration momentum.  For vanishing hel the coefficient is known to
be zero so the integral need not be evaluated.  f0args are the same
arguments as for F0, k1..k5 are the external momenta and m1..m6 are
the masses squared."

MuTildeSq::usage =
"MuTildeSq represents mu-tilde squared in the OPP calculation."

EndPackage[]


BeginPackage["Form`"]

{ q1, intM, paveM, cutM, extM, numM, qfM, qcM, dm4M,
  dotM, abbM, fermM, sunM, helM, powM, addM, mulM, subM,
  d$$, e$$, i$$, dummy$$, g5M, g6M, g7M, dirM, sM }

EndPackage[]


(* symbols from the model files live in Global` *)

{ DiracMatrix, DiracSlash, ChiralityProjector,
  DiracSpinor, MajoranaSpinor, DiracObject,
  PolarizationVector, PolarizationTensor,
  MetricTensor, LeviCivita, FourVector, ScalarProduct,
  Lorentz, Lorentz4, EpsilonScalar,
  SUNT, SUNF, SUNTSum, SUNEps, Colour, Gluon }

TreeCoupling::usage =
"TreeCoupling[from -> to] calculates the tree-level contribution to
the process from -> to.  TreeCoupling[..., opt] specifies InsertFields
options to be used in the computation."

SelfEnergy::usage =
"SelfEnergy[from -> to, mass] calculates the self-energy with incoming
particle from and outgoing particle to, taken at k^2 = mass^2. 
SelfEnergy[f] calculates the self-energy of particle f on its mass
shell.  SelfEnergy[..., opt] specifies InsertFields options to be used
in the computation."

DSelfEnergy::usage =
"DSelfEnergy[from -> to, mass] calculates the derivative with respect to
k^2 of the self-energy with incoming particle from and outgoing particle
to, taken at k^2 = mass^2.  DSelfEnergy[f] calculates the self-energy of
particle f on its mass shell.  DSelfEnergy[..., opt] specifies
InsertFields options to be used in the computation."

ReTilde::usage =
"ReTilde[expr] takes the real part of loop integrals occurring in expr."

ImTilde::usage =
"ImTilde[expr] takes the imaginary part of loop integrals occurring in
expr."

LVectorCoeff::usage =
"LVectorCoeff[expr] returns the coefficient of DiracChain[6, k[1]]
(= k1slash omega_-) in expr."

RVectorCoeff::usage =
"RVectorCoeff[expr] returns the coefficient of DiracChain[7, k[1]]
(= k1slash omega_+) in expr."

LScalarCoeff::usage =
"LScalarCoeff[expr] returns the coefficient of DiracChain[7] (= omega_-)
in expr."

RScalarCoeff::usage =
"RScalarCoeff[expr] returns the coefficient of DiracChain[6] (= omega_+)
in expr."

SEPart::usage =
"SEPart[p, se] returns part p of self-energy se.  It is applied during
the calculation of a renormalization constant, where p is one of
LVectorCoeff, RVectorCoeff, LScalarCoeff, RScalarCoeff for fermion
self-energies, and Identity otherwise."

MassRC::usage =
"MassRC[f] computes the one-loop mass renormalization constant of field
f.  MassRC[f1, f2] computes the one-loop mass renormalization constant
for the f1-f2 transition.  For fermions the output is a list
{left-handed RC, right-handed RC}.  MassRC[..., opt] specifies
InsertFields options to be used in the computation."

FieldRC::usage =
"FieldRC[f] computes the one-loop field renormalization constant of
field f.  FieldRC[f1, f2] computes the one-loop field renormalization
constant for the f1-f2 transition.  FieldRC[f1, f2, c] subtracts c from
the self-energy entering into the calculation.  For fermions the output
is a list {left-handed RC, right-handed RC}.  FieldRC[..., opt]
specifies InsertFields options to be used in the computation."

TadpoleRC::usage =
"TadpoleRC[f] computes the one-loop tadpole renormalization constant of
field f.  TadpoleRC[..., opt] specifies InsertFields options to be used
in the computation."

WidthRC::usage =
"WidthRC[f] computes the one-loop width of field f.  WidthRC[..., opt]
specifies InsertFields options to be used in the computation."


BeginPackage["FormCalc`",
  {"FeynArts`", "LoopTools`", "OPP`", "Form`", "Global`"}]

(* some internal symbols must be visible for FORM/ReadForm *)

{ SUNSum, ReadForm, ReadFormDebug, FormExpr }

(* some internal symbols made visible for debugging *)

LoopIntegral = Join[PaVeIntegral, Blank/@ CutIntegral]

{ FormKins, FormProcs, KinFunc, InvSum, PairRules,
  CurrentProc, LastAmps, DenList, DenMatch, FormSetup }


(* symbols appearing in the output *)

Amp::usage =
"Amp[proc][expr1, expr2, ...] is the result of the calculation of
diagrams of the process proc.  The result is divided into parts expr1,
expr2, ..., such that index sums (marked by SumOver) apply to the whole
of each part."

Den::usage =
"Den[p, m] stands for 1/(p - m).  Note that in contrast to
PropagatorDenominator, p and m are the momentum and mass *squared*.
Den[p, m, d] is the denominator raised to the power d."

Num::usage =
"Num[expr] contains the numerator of a loop integral as a function of
the loop momentum.  This representation is used when calculating loop
integrals by the OPP packages."

DiracChain::usage =
"DiracChain[objs] is a chain of Dirac matrices contracted with the given
objects.  The integers 1, 5, 6, and 7 appearing as first argument denote
1, gamma_5, (1 + gamma_5)/2, and (1 - gamma_5)/2, respectively."

WeylChain::usage =
"WeylChain[objs] is a chain of sigma matrices contracted with the given
objects.  The integers 6, 7 respectively denote upper and lower indices
at the given position, and -1 stands for epsilon, the spinor metric."

Spinor::usage =
"Spinor[p, m, s] is a spinor with momentum p and mass m, i.e. a solution
of the Dirac equation (pslash + s m) Spinor[p, m, s] = 0.  On screen,
particle spinors (s = 1) are printed as u[p, m], antiparticle spinors
(s = -1) as v[p, m].
Inside a WeylChain, Spinor denotes a 2-component Weyl spinor and has two
additional arguments: Spinor[p, m, s, d, e].  An undotted spinor has
d = 1, a dotted d = 2; e indicates contraction with the spinor metric.  
Whether it corresponds to the upper or lower half of the 4-component
Dirac spinor is determined by the index convention of the WeylChain
(fixed by arguments 6 or 7), propagated to the position of the spinor."

e::usage =
"e[i] is the ith polarization vector."

ec::usage =
"ec[i] is the conjugate of the ith polarization vector."

z::usage =
"z[i] is the ith polarization vector in D - 4 dimensions."

zc::usage =
"zc[i] is the conjugate of the ith polarization vector in D - 4
dimensions."

eT::usage =
"eT[i] is the ith polarization tensor."

eTc::usage =
"eTc[i] is the conjugate of the ith polarization tensor."

k::usage =
"k[i] is the ith momentum."

nul::usage =
"nul is the zero vector.  It is kept only inside IGram to trace singular
Gram determinants."

SUNN::usage =
"SUNN specifies the N in SU(N), i.e. the number of colours."

S::usage =
"S is the Mandelstam variable s.  If k1 and k2 are the incoming momenta,
S = (k1 + k2)^2 = S12."

T::usage =
"T is the Mandelstam variable t.  If k1 denotes the first incoming and
k3 the first outgoing momentum, T = (k1 - k3)^2 = T13."

U::usage =
"U is the Mandelstam variable u.  If k2 denotes the first incoming and
k3 the second outgoing momentum, U = (k2 - k3)^2 = T23."

Dminus4::usage =
"Dminus4 represents the difference D - 4."

Dminus4Eps::usage =
"Dminus4Eps represents the Dminus4 arising from the contraction of
Levi-Civita tensors in HelicityME and PolarizationSum.  Choosing it
different from Dminus4 is not particularly consistent and used for
checking results only."

Dminus4Eps = Dminus4

IGram::usage =
"IGram[expr] contains the denominator arising from the reduction of a
tensor-coefficient function.  It is equivalent to 1/expr but is kept
separate for further simplification."

PowerOf::usage =
"PowerOf is an auxiliary function with the help of which it is easy
to split the amplitude into monomials in the couplings constants.
To this end, subject the FeynArts amplitude to the substitution
g -> g PowerOf[g] for all coupling constants g.  Each element
the Amp returned by CalcFeynAmp contains a single PowerOf only,
from which the power in the coupling constants can be read off."


(* ToFeynAmp, DeclareProcess, CalcFeynAmp and their options *)

ToFeynAmp::usage =
"ToFeynAmp[amps] converts hand-typed amplitudes into (approximate)
FeynArts conventions such that the amplitude can be processed by
DeclareProcess and CalcFeynAmp."

DeclareProcess::usage =
"DeclareProcess[amps] sets up internal definitions for the computation
of the amplitudes amps."

OnShell::usage =
"OnShell is an option of DeclareProcess.  It specifies whether FORM
should put the external particles on their mass shell, i.e. apply
ki^2 = mi^2."

Invariants::usage =
"Invariants is an option of DeclareProcess.  It specifies whether FORM
should introduce kinematical invariants, like the Mandelstam variables
for a 2 -> 2 process."

Transverse::usage =
"Transverse is an option of DeclareProcess.  It specifies whether FORM
should apply the transversality relations for polarization vectors
(ei.ki = 0)."

Normalized::usage =
"Normalized is an option of DeclareProcess.  It specifies whether FORM
should apply the normalization of polarization vectors (ei.ei^* = -1)."

InvSimplify::usage =
"InvSimplify is an option of DeclareProcess.  It specifies whether FORM
should try to simplify combinations of invariants as much as possible."

Antisymmetrize::usage =
"Antisymmetrize is an option of DeclareProcess.  It specifies whether
Dirac chains are antisymmetrized."

FormAmp::usage =
"FormAmp[proc][amps] contains a preprocessed form of the FeynArts
amplitudes amps for process proc, to be calculated by CalcFeynAmp."

CalcFeynAmp::usage =
"CalcFeynAmp[amps] calculates the Feynman amplitudes given in amps.  The
resulting expression is broken up into categories which are returned in
an Amp object."

CalcLevel::usage =
"CalcLevel is an option of CalcFeynAmp.  It specifies the level (Classes
or Particles) at which to calculate the amplitudes. Automatic takes
Classes level, if available, otherwise Particles."

Dimension::usage =
"Dimension is an option of CalcFeynAmp, HelicityME, and PolarizationSum. 
It specifies the space-time dimension in which to perform the
calculation and can take the values D, where dimensional regularization
is used, and 4, where constrained differential renormalization is used. 
The latter method is equivalent to dimensional reduction at the one-loop
level.  Dimension -> 0 retains the Dminus4 terms, i.e. does not emit
local (rational) terms."

MomElim::usage =
"MomElim is an option of CalcFeynAmp.  It controls in which way
momentum conservation is used to eliminate momenta.  False performs no
elimination, an integer between 1 and the number of legs substitutes
the specified momentum in favour of the others, and Automatic tries all
substitutions and chooses the one resulting in the fewest terms."

DotExpand::usage =
"DotExpand is an option of CalcFeynAmp.  It controls whether the terms
collected for momentum elimination are expanded again.  This prevents
kinematical invariants from appearing in the abbreviations but typically
leads to poorer simplification of the amplitude."

NoCostly::usage =
"NoCostly is an option of CalcFeynAmp.  Useful for 'tough' amplitudes,
NoCostly -> True turns off potentially time-consuming transformations
in FORM."

FermionChains::usage =
"FermionChains is an option of CalcFeynAmp.  It can take the three
values Chiral, VA, and Weyl, which specify how fermion chains are
returned by CalcFeynAmp.  Chiral and VA both select (4-dimensional)
Dirac chains, where the chiral decomposition is taken for Chiral and the
vector/axial-vector decomposition for VA.  Weyl selects (2-dimensional)
Weyl chains."

Chiral::usage =
"Chiral is a possible value for the FermionChains option of CalcFeynAmp. 
It instructs CalcFeynAmp to return fermion chains as left- and
right-handed (4-dimensional) Dirac chains."

VA::usage =
"VA is a possible value for the FermionChains option of CalcFeynAmp. 
It instructs CalcFeynAmp to return fermion chains as vector and
axial-vector parts of (4-dimensional) Dirac chains."

Weyl::usage =
"Weyl is a possible value for the FermionChains option of CalcFeynAmp. 
It instructs CalcFeynAmp to return fermion chains as (2-dimensional)
Weyl chains."

FermionOrder::usage =
"FermionOrder is an option of CalcFeynAmp.  It determines the ordering
of Dirac spinor chains in the output, i.e. requires FermionChains ->
Chiral or VA.  Possible values are None, Fierz, Automatic, Colour, or
an explicit ordering, e.g. {2, 1, 4, 3}.  None applies no reordering. 
Fierz applies the Fierz identities twice, thus simplifying the chains
but keeping the original order.  Colour applies the ordering of the
external colour indices (after simplification) to the spinors. 
Automatic chooses a lexicographical ordering."

Fierz::usage =
"Fierz is a possible value for the FermionOrder option of CalcFeynAmp. 
It instructs CalcFeynAmp to apply the Fierz identities twice, thus
simplifying the chains but keeping the original spinor order."

Colour::usage =
"Colour is a possible value for the FermionOrder option of CalcFeynAmp. 
It instructs CalcFeynAmp to bring the spinors into the same order as the
external colour indices, i.e. \"Fermion flow follows colour flow\"."

Evanescent::usage =
"Evanescent is an option of CalcFeynAmp.  It introduces factors
of the form Evanescent[original operator, Fierzed operator]
with the help of which one can detect problems due to the application
of the Fierz identities."

InsertionPolicy::usage =
"InsertionPolicy is an option of CalcFeynAmp.  Begin specifies that the
insertions shall be applied at the beginning of the FORM code (this
ensures that all loop integrals are fully symmetrized).  Default applies
them after simplifying the generic amplitudes (this is fastest).  A
positive integer does the same, except that insertions with a LeafCount
larger than that integer are inserted only after the amplitude comes
back from FORM (this is a workaround for the rare cases where the FORM
code aborts due to very long insertions)."

SortDen::usage =
"SortDen is an option of CalcFeynAmp.  It determines whether the
denominators of loop integrals shall be sorted.  This is usually done to
reduce the number of loop integrals appearing in an amplitude."

CombineDen::usage =
"CombineDen is an option of CalcFeynAmp.  It determines whether to
combine OPP integrals with common denominators, as in:
N2/(D0 D1) + N3/(D0 D1 D2) -> (N2 D2 + N3)/(D0 D1 D2).
True/False turns combining integrals on/off, Automatic combines integrals
only if the rank of the combined numerator is not larger than the sum
of the individual numerator ranks (which is faster)."

PaVeReduce::usage =
"PaVeReduce is an option of CalcFeynAmp.  False retains the one-loop
tensor-coefficient functions.  Raw reduces them to scalar integrals but
keeps the Gram determinants in the denominator in terms of dot products. 
True simplifies the Gram determinants using invariants."

SimplifyQ2::usage =
"SimplifyQ2 is an option of CalcFeynAmp.  It controls simplification
of terms involving the integration momentum squared.  If set to True,
powers of q1.q1 in the numerator are cancelled by a denominator, except
for OPP integrals, where conversely lower-N integrals are put on a
common denominator with higher-N integrals to reduce OPP calls, as in:
N2/(D0 D1) + N3/(D0 D1 D2) -> (N2 D2 + N3)/(D0 D1 D2)."

OPP::usage =
"OPP is an option of CalcFeynAmp.  It specifies an integer N starting
from which an N-point function is treated with OPP methods.  For
example, OPP -> 4 means that A, B, C functions are reduced with
Passarino-Veltman and D and up with OPP.
A negative N indicates that the rational terms for the OPP integrals
shall be added analytically whereas otherwise their computation is left
to the OPP package."

OPPQSlash::usage =
"OPPQSlash is an option of CalcFeynAmp.  The integration momentum q1
starts life as a D-dimensional object.  In OPP, any q1.q1 surviving
SimplifyQ2 are substituted by q1.q1 - MuTildeSq, after which q1 is
considered 4-dimensional.  The dimensionful scale MuTildeSq enables the
OPP libraries to reconstruct the rational terms.
OPPQSlash extends this treatment to the q1-slash on external fermion
chains, i.e. also substitutes q1-slash -> q1-slash + I gamma_5 MuTilde,
where odd powers of MuTilde are eventually set to zero."

Gamma5Test::usage =
"Gamma5Test is an option of CalcFeynAmp.  If set to True, each gamma_5
is multiplied by (1 + Gamma5Test (D - 4)) and the dependence on
Gamma5Test in the final result should vanish."

Gamma5ToEps::usage =
"Gamma5ToEps is an option of CalcFeynAmp.  It substitutes gamma_5 by
1/24 eps(mu, nu, ro, si) gamma(mu) gamma(nu) gamma(ro) gamma(si) in
fermion traces."

NoExpand::usage =
"NoExpand is an option of CalcFeynAmp.  NoExpand -> {sym1, sym2, ...}
specifies that sums containing any of sym1, sym2, ... are not expanded
during the FORM calculation."

NoBracket::usage =
"NoBracket is an option of CalcFeynAmp.  NoBracket -> {sym1, sym2, ...} 
specifies that sym1, sym2, ... are not collected inside a multiplication 
bracket during the FORM calculation."

MomRules::usage =
"MomRules is an option of CalcFeynAmp.  It specifies a set of rules
for transforming momenta.  The notation is that of the final amplitude,
i.e. k1,...,kn for the momenta, e1,...,en for the polarization vectors."

PreFunction::usage =
"PreFunction is an option of CalcFeynAmp.  It specifies a function to be
applied to each amplitude before any simplification is made.  This
option is typically used to apply a function to all amplitudes in a
calculation, even in indirect calls to CalcFeynAmp, such as through
CalcRenConst."

PostFunction::usage =
"PostFunction is an option of CalcFeynAmp.  It specifies a function to
be applied to each amplitude after all simplifications have been made. 
This option is typically used to apply a function to all amplitudes in a
calculation, even in indirect calls to CalcFeynAmp, such as through
CalcRenConst."

FileTag::usage =
"FileTag is an option of CalcFeynAmp.  It specifies the middle part of
the name of the FORM file, as in fc-TAG-1.frm."

EditCode::usage =
"EditCode is a debugging option of CalcFeynAmp, HelicityME, and
PolarizationSum.  It determines editing of the intermediate FORM code. 
True invokes the $Editor command, which is supposed to detach from
FormCalc, i.e. the FORM process continues with the unedited code. 
Modal invokes the $EditorModal command, which is supposed to be
modal (non-detached), i.e. continues only after the editor is closed,
thus continuing with possibly modified FORM code."

RetainFile::usage =
"RetainFile is a debugging option of CalcFeynAmp, HelicityME, and
PolarizationSum.  When set to True, the temporary file containing the
FORM code is not removed after running FORM."


(* abbreviationing-related functions *)

Abbr::usage =
"Abbr[] returns a list of all abbreviations introduced so far.
Abbr[patt] returns a list of all abbreviations including the pattern
patt.  Patterns prefixed by ! (Not) are excluded."

GenericList::usage =
"GenericList[] returns a list of the substitutions made for the
computation of generic amplitudes."

ClearProcess::usage =
"ClearProcess[] is necessary to clear internal definitions before
calculating a process with a different kinematical set-up."

ZapFunction::usage =
"ZapFunction is an option of ClearProcess and ClearSubexpr and
determines the function used to clear the definitions of symbols
introduced for abbreviations."

RegisterAbbr::usage =
"RegisterAbbr[abbr] registers a list of abbreviations such that
future invocations of CalcFeynAmp will make use of them.  Note that
abbreviations introduced for different processes are in general not
compatible."

Abbreviate::usage =
"Abbreviate[expr, x] introduces abbreviations for subexpressions
in expr.  If x is an integer, abbreviations are introduced for sums
from level x downward.  Otherwise, x is taken as a function name and
abbreviations are introduced for all subexpressions for which
x[subexpr] returns True."

AbbrevSet::usage =
"AbbrevSet[expr] sets up the AbbrevDo function that introduces
abbreviations for subexpressions.  The expression given here is not
itself abbreviated but used for determining the summation indices."

AbbrevDo::usage =
"AbbrevDo[expr, x] invokes the abbreviationing function defined with
AbbrevSet on expr.  If x is an integer, abbreviations are introduced
from level x downward.  Otherwise, x is taken as a function name and
abbreviations are introduced for all subexpressions for which
x[subexpr] returns True."

Deny::usage =
"Deny is an option of Abbreviate.  It specifies items which must not be
included in abbreviations."

Fuse::usage =
"Fuse is an option of Abbreviate.  It specifies whether adjacent items
for which the selection function is True should be fused into one
abbreviation."

Preprocess::usage =
"Preprocess is an option of Abbreviate.  It specifies a function to be
applied to all subexpressions before introducing abbreviations for
them."

$AbbPrefix::usage =
"$AbbPrefix specifies the prefix for subexpressions introduced by
Abbreviate, i.e. the Sub in Sub123."

Subexpr::usage =
"Subexpr[] returns a list of all subexpressions introduced by
Abbreviate.
Subexpr[args] executes Abbreviate[args] locally, i.e. without
registering the subexpressions permanently and returns a list of the
form {Subexpr[], Abbreviate[args]}."

ClearSubexpr::usage =
"ClearSubexpr[] clears the internal definitions of the subexpressions
introduced by Abbreviate."

RegisterSubexpr::usage =
"RegisterSubexpr[subexpr] registers a list of subexpressions such that
future invocations of Abbreviate will make use of them."

OptimizeAbbr::usage =
"OptimizeAbbr[abbr] optimizes the abbreviations in abbr by eliminating
common subexpressions."

$OptPrefix::usage =
"$OptPrefix specifies the prefix for additional abbreviations introduced
by OptimizeAbbr, i.e. the Opt in Opt123."

SubstSimpleAbbr::usage =
"SubstSimpleAbbr[{expr, abbr}] removes `simple' abbreviations from abbr
and substitutes them back into expr. 
SubstSimpleAbbr[CodeExpr[var, tmpvar, expr]] removes `simple'
abbreviations from expr and deletes them from the variable lists var
and tmpvar."

ExtractInt::usage =
"ExtractInt[expr] extracts the loop integrals from expr for evaluation
with LoopTools/OPP.  It output is a list {lint, abbrexpr}, where lint
is the list of loop integrals and abbrexpr is expr with the loop
integrals substituted by the identifiers in integrals."

Pair::usage =
"Pair[a, b] represents the contraction of the two four-vectors or
Lorentz indices a and b."

Eps::usage =
"Eps[a, b, c, d] represents -I times the antisymmetric Levi-Civita
tensor epsilon_{abcd}.  The sign convention is epsilon^{0123} = +1."

ToSymbol::usage =
"ToSymbol[s...] concatenates its arguments into a new symbol."

NewSymbol::usage =
"NewSymbol[stub, 0] creates a new symbol of the form stubN, where
N is the integer SymbolNumber[stub] + 1.
NewSymbol[stub] furthermore guarantees that stubN is presently not
used elsewhere by increasing N as necessary."

SymbolNumber::usage =
"SymbolNumber[stub] gives the running number last used for creating
a new symbol with prefix stub."

ToArray::usage =
"ToArray[s] turns the symbol s into an array reference by taking it
apart into letters and digits, e.g. Var1234 becomes Var[1234]. 
ToArray[expr, s1, s2, ...] turns all occurrences of the symbols s1NNN,
s2NNN, etc. in expr into s1[NNN], s2[NNN], etc."

Renumber::usage =
"Renumber[expr, var1, var2, ...] renumbers all var1[n], var2[n], ... in
expr."

MaxDims::usage =
"MaxDims[args] returns a list of all distinct functions in args with the
highest indices that appear, e.g. MaxDims[foo[1, 2], foo[2, 1]] returns
{foo[2, 2]}."

Keep::usage =
"Keep[expr, name, path] loads path/name.m if it exists, otherwise it
evaluates expr and stores it (together with the Abbr[] and Subexpr[])
in that file.  path is optional and defaults to $KeepDir. 
Keep[lhs = rhs] is short for lhs = Keep[rhs, \"lhs\"]."

$KeepDir::usage =
"$KeepDir specifies the default directory for storing intermediate
expressions with Keep."


(* miscellaneous functions *)

FormPre::usage =
"FormPre is a function executed immediately before invoking FORM.
It receives the raw amplitudes as argument and is used to initialize
the simplification functions applied to the amplitude when coming back
from FORM."

FormSub::usage =
"FormSub is a function applied to subexpressions extracted by FORM
for simplification."

FormDot::usage =
"FormDot is a function applied to combinations of dot products
extracted by FORM before abbreviationing."

FormMat::usage =
"FormMat is a function applied to the coefficients of matrix elements
(Mat[...]) in the FORM output."

FormNum::usage =
"FormNum is a function applied to numerator functions in the FORM
output (OPP only)."

FormQC::usage =
"FormQC is a function applied to the loop-momentum-independent parts
of the OPP numerator in the FORM output."

RCSub::usage =
"RCSub is a simplification function, it substitutes FormSub during the
calculation of the renormalization constants."

RCInt::usage =
"RCInt is a simplification function applied to the coefficients of loop
integrals during the calculation of the renormalization constants."

DotSimplify::usage =
"DotSimplify[f1, f2][expr] simplifies expr using f2 if expr contains any
of the objects listed in the NoBracket option (i.e. during the execution
of CalcFeynAmp or PolarizationSum), and using f1 otherwise."

TermCollect::usage =
"TermCollect[expr] pulls common factors out of sums."

DenCancel::usage =
"DenCancel[expr] tries to cancel Den[p, m] against factors (p - m)
in the numerator."

DenCollect::usage =
"DenCollect[expr] collects terms in expr whose denominators are
identical up to a numerical constant.  DenCollect[expr, wrap] applies
wrap to the collected numerators."

Pool::usage =
"Pool[expr] combines terms with common factors.  Unlike Factor, it looks
at the terms pairwise and can thus do a b + a c + d -> a (b + c) + d
fast.  Unlike Simplify, it does not modify b and c.  Pool[expr, wrap]
applies wrap to the (b + c) part."

ApplyUnitarity::usage =
"ApplyUnitarity[expr, mat, d, simp] simplifies expr by exploiting the
unitarity of the d-dimensional matrix mat.  The optional argument simp
specifies the simplification function to use internally and defaults to
FullSimplify."

OnSize::usage =
"OnSize[n1, f1, n2, f2, ..., fdef][expr] returns f1[expr] if
LeafCount[expr] < n1, f2[expr] if LeafCount[expr] < n2, etc., and
fdef[expr] if the expression is still larger.  fdef can take the
special value Map which means that fdef[expr] recursively applies
the entire OnSize function to the parts of expr.  If omitted, fdef
defaults to Identity."

OffShell::usage =
"OffShell[amps, i -> mi, ...] returns the FeynAmpList amps with the
mass of the ith external particle set to mi.  This will in general take
particle i off its mass shell since now ki^2 = mi^2 is fulfilled with
the new value of mi."

Combine::usage =
"Combine[amp1, amp2, ...] combines the amplitudes amp1, amp2, ... which
can be either FeynAmpList or Amp objects, i.e. Combine works before and
after CalcFeynAmp."

ExpandSums::usage =
"ExpandSums[expr] turns all pieces of expr multiplied with SumOver
into an actual Sum.  ExpandSums[expr, h] uses h instead of Sum."

MultiplyDiagrams::usage =
"MultiplyDiagrams[func][amp] multiplies the diagrams in amp with the
factor returned by the function func.  The latter is invoked for each
diagram either as func[amplitude] (for a fully inserted diagram), or as
func[generic amplitude, insertion]."

TagDiagrams::usage =
"TagDiagrams[amp] tags each diagram in amp with an identifier of the
form Diagram[number], where number runs sequentially through the
diagrams at all levels.  This makes it possible to locate the
contribution of individual diagrams in the final CalcFeynAmp output. 
TagDiagrams[amp, tag] uses tag rather than Diagram."

TagCollect::usage =
"TagCollect[expr, tag, f] collects expr with respect to powers of tag
and applies f to the term linear in tag.  This function is typically
used to apply f to a tagged part of an expression."

Diagram::usage =
"Diagram[number] is the identifier used to tag a single diagram by
TagDiagrams."

DiagramType::usage =
"DiagramType[diag] returns the number of denominators not containing
the integration momentum."

FermionicQ::usage =
"FermionicQ[diag] gives True for a diagram containing fermions and
False otherwise."

IndexIf::usage =
"IndexIf[cond, a, b] is identical to If[cond, a, b], except that a
and b are not held unevaluated.  IndexIf[cond, a] is equivalent to
IndexIf[cond, a, 0] and IndexIf[cond1, a1, cond2, a2, ...] is
equivalent to IndexIf[cond1, a1, IndexIf[cond2, a2, ...]]."

MapIf::usage =
"MapIf[foo, i] maps foo over the expression parts only (i.e. not over
the conditions) if i is an IndexIf expression.  For all other types of
expressions, MapIf is equivalent to Map."

IndexDiff::usage =
"IndexDiff[i, j] is the same as 1 - IndexDelta[i, j]."

ToIndexIf::usage =
"ToIndexIf[expr] converts all IndexDeltas and IndexDiffs in expr to
IndexIf, which will be written out as if-statements in the generated
code.  ToIndexIf[expr, patt] operates only on indices matching patt.
If patt is a string, e.g. \"Glu*\", it is first expanded to all
matching symbols."

Neglect::usage =
"Neglect[sym] = 0 makes FORM replace sym = 0 except when it appears in
negative powers or in loop integrals."

Square::usage =
"Square[m] = m2 makes FORM replace all m^2 by m2."

NClear::usage =
"NClear[patt] clears the NValues of all symbols matching patt. 
NClear[] is equivalent to NClear[\"Global`*\"]."

ColourSimplify::usage =
"ColourSimplify[expr] simplifies the colour objects in expr.
ColourSimplify[plain, conj] simplifies the colour objects in
(plain conj^*)."

ColourGrouping::usage =
"ColourGrouping[tops] returns a list of parts of the inserted topologies
tops, grouped according to their colour structures."


(* FeynCalc compatibility functions *)

FeynCalcGet::usage =
"FeynCalcGet[mask] reads files produced with FeynCalc.  mask is taken
as input to the Mathematica function FileNames, so it might be
FeynCalcGet[\"file.m\"] or FeynCalcGet[\"*.m\"] or
FeynCalcGet[\"*.m\", \"~/feyncalcfiles\"]."

FeynCalcPut::usage =
"FeynCalcPut[expr, file] writes expr to file in FeynCalc format."


(* finiteness checks *)

UVDivergentPart::usage =
"UVDivergentPart[expr] returns expr with all loop integrals replaced
by their UV-divergent part.  The divergence itself is denoted by
Divergence."

UVSeries::usage =
"UVSeries[expr] expands expr into a Laurent series in Dminus4.
UVSeries[expr, pow] returns the coefficient of Dminus4^pow."

Divergence::usage =
"Divergence represents the dimensionally regularized divergence
2/(4 - D) of loop integrals.  It is used by the function
UVDivergentPart."

Finite::usage =
"Finite is a symbol with which the local terms (resulting from D times
divergent integral) are multiplied with.  This is to be able to remove
these terms when evaluating the eps^-1 or eps^-2 coefficients."


(* matrix elements *)

HelicityME::usage =
"HelicityME[plain, conj] calculates the helicity matrix elements for all
combinations of spinor chains that appear in the expression
(plain conj^*).  Terms of this kind arise in the calculation of the
squared matrix element, where typically plain is the one-loop result
and conj the tree-level expression.  The arguments do not necessarily
have to be amplitudes since they are only used to determine which spinor
chains to select from the abbreviations.  The symbol All can be used to
select all spinor chains currently defined in the abbreviations."

ColourME::usage =
"ColourME[plain, conj] calculates the colour matrix elements.  ColourME
is very similar to HelicityME, except that it computes the matrix
elements for SU(N) objects, not for spinor chains."

All::usage =
"All as an argument of HelicityME and ColourME indicates that all spinor
chains or SUNT objects currently defined in the abbreviations should be
used instead of just those appearing in the argument."

Hel::usage =
"Hel[i] is the helicity of the ith external particle.  It can take the
values +1, 0, -1, where 0 stands for an unpolarized particle."

s::usage =
"s[i] is the ith helicity reference vector."

Mat::usage =
"Mat[Fi SUNi] is a matrix element in an amplitude, i.e. an amplitude is
a linear combination of Mat objects.\n Mat[Fi, Fj] appears in the
squared matrix element and stands for the product of the two arguments,
Fi Fj^*.  Such expressions are calculated by HelicityME and ColourME."

Lor::usage =
"Lor[i] is a contracted Lorentz index in a product of Dirac chains."

SquaredME::usage =
"SquaredME[plain, conj] returns the matrix element
(plain Conjugate[conj]).  This performs a nontrivial task only for
fermionic amplitudes: the product of two fermionic amplitudes\n
    M1 = a1 F1 + a2 F2 + ... and\n
    M2 = b1 F1 + b2 F2 + ... is returned as\n
    M1 M2^* = a1 b1^* Mat[F1, F1] + a2 b1^* Mat[F2, F1] + ...\n
The special case of plain === conj can be written as SquaredME[plain]
which is of course equivalent to SquaredME[plain, plain]."

RealQ::usage =
"RealQ[sym] is True if sym represents a real quantity which means in
particular that Conjugate[sym] = sym."

PolarizationSum::usage =
"PolarizationSum[expr] sums expr over the polarizations of external
gauge bosons.  It is assumed that expr is the squared amplitude into
which the helicity matrix elements have already been inserted. 
Alternatively, expr may also be given as an amplitude directly, in which
case PolarizationSum will first invoke SquaredME and HelicityME (with
Hel[_] = 0) to obtain the squared amplitude."

GaugeTerms::usage =
"GaugeTerms is an option of PolarizationSum.  With GaugeTerms -> False,
the gauge-dependent terms in the polarization sum, which should
eventually cancel in gauge-invariant subsets of diagrams, are omitted
from the beginning."

eta::usage =
"eta[i] is a vector that defines a particular gauge via Pair[eta[i],
e[i]] = 0.  It is introduced by PolarizationSum for massless particles,
where the sum over e[i][mu] ec[i][nu] is gauge dependent.  Only for
gauge-invariant subsets of diagrams should the dependence on eta[i]
cancel.  eta obeys Pair[eta[i], k[i]] != 0 and Pair[eta[i], e[i]] = 0."


(* writing out code *)

SetupCodeDir::usage =
"SetupCodeDir[dir] installs the driver programs necessary to compile the
code generated by WriteSquaredME and WriteRenConst in the directory dir.
Customized versions of the drivers are taken from the directory pointed
to by the Drivers option and take precedence over the default versions
from $DriversDir.  Drivers already in dir are not overwritten."

Drivers::usage =
"Drivers is an option of SetupCodeDir.  Drivers points to a directory
containing customized versions of the driver programs necessary for
compiling the generated code.  This directory need not contain all
driver programs: files not contained therein are taken from the default
directory $DriversDir."

WriteSquaredME::usage =
"WriteSquaredME[tree, loop, me, abbr, ..., dir] writes out code to
compute the squared matrix element for a process whose tree-level and
one-loop contributions are given in the first and second argument,
respectively.  All further arguments except the last specify the
necessary matrix elements and abbreviations.  The last argument dir
finally gives the path to write the generated code to."

ExtraRules::usage =
"ExtraRules is an option of WriteSquaredME.  Rules given here will be
applied before the loop integrals are abbreviated."

TreeSquare::usage =
"TreeSquare is an option of WriteSquaredME.  It specifies whether to
add the square of the tree-level amplitude, |M_0|^2, to the result 
in the SquaredME subroutine."

LoopSquare::usage =
"LoopSquare is an option of WriteSquaredME.  It specifies whether to
add the square of the 1-loop amplitude, |M_1|^2, to the result in the
SquaredME subroutine.  This term is of order alpha^2 with respect to the
tree-level contribution, |M_0|^2.  Usually one takes into account only
the interference term, 2 Re M_0^* M_1, which is of order alpha."

Folder::usage =
"Folder is an option of WriteSquaredME and WriteRenConst.  It specifies
the folder into which the generated files are written."

FilePrefix::usage =
"FilePrefix is an option of WriteSquaredME and WriteRenConst.  It
specifies a string to be prepended to the filenames of the generated
code."

SymbolPrefix::usage =
"SymbolPrefix is an option of WriteSquaredME and WriteRenConst.  It
specifies a string which is prepended to externally visible symbols in
the generated code to prevent collision of names when several processes
are linked together."

FileIncludes::usage =
"FileIncludes is an option of WriteSquaredME and WriteRenConst. 
It specifies per-file #include statements (or other declarations). 
Like SubroutineIncludes, it accepts also a list of two elements
with which one can separate declarations from (inline) function
definitions."

SubroutineIncludes::usage =
"SubroutineIncludes is an option of WriteSquaredME and WriteRenConst. 
It specifies per-subroutine #include statements (or other declarations).  
If a list of two elements is given, they are included before and after
local variable declarations so that the second may include code such
as inline function definitions."

FileHeader::usage =
"FileHeader is an option of WriteSquaredME and WriteRenConst and
specifies a file header.  This string may contain %f, %d, and %t, which
are substituted at file creation by file name, description, and time
stamp, respectively."


(* renormalization constants *)

FindRenConst::usage =
"FindRenConst[expr] returns a list of all renormalization constants
found in expr including those needed to compute the former."

CalcRenConst::usage =
"CalcRenConst[expr] calculates the renormalization constants appearing
in expr."

WriteRenConst::usage =
"WriteRenConst[expr, dir] calculates the renormalization constants
appearing in expr and generates code from the results.  The resulting
files (the Fortran program itself and the corresponding declarations)
are written to the directory dir.  The names of the files are
determined by the RenConstFile option."

InsertFieldsHook::usage =
"InsertFieldsHook[tops, proc] is the function called by SelfEnergy and
DSelfEnergy to insert fields into the topologies tops for the process
proc.  It is normally equivalent to InsertFields, but may be redefined
to change the diagram content of certain self-energies."

CreateFeynAmpHook::usage =
"CreateFeynAmpHook[diags, opt] is the function called by SelfEnergy and
DSelfEnergy to create the amplitudes for diagrams diags.  It is normally
equivalent to CreateFeynAmp, but may be redefined to modify the
amplitudes."

ClearSE::usage =
"ClearSE[] clears the internal definitions of already calculated
self-energies."

OptPaint::usage =
"OptPaint[ins, True] invokes Paint[ins]. 
OptPaint[ins, pre] invokes Paint[ins], saving the graphics in files
with names constructed from ProcessName[ins], prefixed with the
string pre (which may include a directory). 
OptPaint[ins, pre, suf] further appends suf to the file name. 
OptPaint[ins] executes OptPaint[ins, $PaintSE]."

$PaintSE::usage =
"$PaintSE determines whether SelfEnergy paints the diagrams it generates
to compute the self-energies.  $PaintSE can be True, False, or a string
which indicates that the output should be saved in a PostScript file
instead of being displayed on screen, and is prepended to the filename."

$LongitudinalSE::usage =
"$LongitudinalSE specifies that the longitudinal rather than the
transverse part of the vector-boson self-energies is taken in
SelfEnergy and DSelfEnergy."


(* low-level code output functions *)

ToList::usage =
"ToList[expr] returns a list of summands of expr."

MkDir::usage =
"MkDir[\"dir1\", \"dir2\", ...] makes sure the directory dir1/dir2/...
exists, creating the individual subdirectories dir1, dir2, ... as
necessary."

ToForm::usage =
"ToForm[expr] returns the FORM form of expr as a string."

OpenForm::usage =
"OpenForm[stub] opens a temporary FORM file with a unique name
starting with stub (\"fc\" by default) for writing."

FormFile::usage =
"FormFile[n] returns the name of the temporary FORM file name with
index n."

ToCode::usage =
"ToCode[expr] returns the Fortran or C form of expr as a string."

ToFortran::usage =
"ToFortran has been superseded by ToCode."

SetLanguage::usage =
"SetLanguage[lang] sets the language for source code, currently
\"Fortran\" or \"C\"."

$Code::usage =
"$Code is the current language for writing out source code."

$CodeExt::usage =
"$CodeExt is the current filename extension for source code."

OpenCode::usage =
"OpenCode[file] opens file for writing out source code in the current
output language set by SetLanguage."

OpenFortran::usage =
"OpenFortran[file] opens file for writing in Fortran format."

OpenC::usage =
"OpenC[file] opens file for writing in C99 format."

TimeStamp::usage =
"TimeStamp[] returns a string with the current date and time."

BlockSplit::usage =
"BlockSplit[var -> expr] tries to split the calculation of expr into
subexpressions each of which has a leaf count less than $BlockSize."

FileSplit::usage =
"FileSplit[exprlist, mod, writemod, writeall] splits exprlist into
batches with leaf count less than $FileSize.  If there is only one
batch, writemod[batch, mod] is invoked to write it to file.  Otherwise,
writemod[batch, modN] is invoked on each batch, where modN is mod
suffixed by a running number, and in the end writeall[mod, res] is
called, where res is the list of writemod return values.  The optional
writeall function can be used e.g. to write out a master subroutine
which invokes the individual modules."

FortranNames::usage =
"FortranNames[base, ind] constructs two names out of base and ind,
where the latter are typically indices.  The first is the direct
concatenation of the elements, separated by underscores, and is for
use in file names and similar uncritical places.  The second is
limited to a maximum of $MaxFortranName characters and should be
used for Fortran symbols to comply with compiler limits.  Truncation
is done by leaving out delimiting underscores and if that is not
enough, contracting index names down to 2 characters."

$MaxFortranName::usage =
"$MaxFortranName specifies the maximum length of a symbol name in
Fortran."

RuleAdd::usage =
"RuleAdd[var, expr] is written out by WriteExpr as var -> var + expr."

PrepareExpr::usage =
"PrepareExpr[{var1 -> expr1, var2 -> expr2, ...}] prepares a list of
variable assignments for write-out to a Fortran file. Expressions with a
leaf count larger than $BlockSize are split into several pieces, as in\n
\tvar = part1\n\
\tvar = var + part2\n\
\t...\n
thereby possibly introducing temporary variables for subexpressions. 
The output is a CodeExpr[vars, tmpvars, exprlist] object, where vars
are the original and tmpvars the temporary variables introduced by
PrepareExpr."

WriteExpr::usage =
"WriteExpr[file, exprlist] writes a list of variable assignments in
Fortran format to file.  The exprlist can either be a CodeExpr object
or a list of expressions of the form {var1 -> expr1, var2 -> expr2,
...}, which is first converted to a CodeExpr object using
PrepareExpr.  WriteExpr returns a list of the subexpressions that were
actually written."

HornerStyle::usage =
"HornerStyle is an option of WriteExpr.  It specifies whether
expressions are arranged in Horner form before writing them out as
Fortran code."

FinalCollect::usage =
"FinalCollect is an option of WriteExpr.  It specifies whether common
factors are collected in the final expression, just before write-out to
Fortran or C."

FinalFunction::usage =
"FinalFunction is an option of WriteExpr.  It specifies a function to be
applied to the final expressions, just before write-out to Fortran or C."

Type::usage =
"Type is an option of WriteExpr.  If a string is given, e.g. Type ->
\"double precision\", WriteExpr writes out declarations of that type for
the given expressions.  Otherwise no declarations are produced."

TmpType::usage =
"TmpType is an option of WriteExpr.  It is the counterpart of Type for
the temporary variables.  TmpType -> Type uses the settings of the Type
option."

IndexType::usage =
"IndexType is an option of WriteExpr.  It is the counterpart of Type
for do-loop indices.  IndexType -> Type uses the settings of the Type
option."

RealArgs::usage =
"RealArgs is an option of WriteExpr.  It specifies a list of functions
whose numerical arguments must be of a guaranteed type (usually real). 
For example, if foo expects a real argument, it must be invoked as
foo(0D0), not foo(0) in Fortran.
RealArgs[foo] := ... defines the actual conversion for foo.  The
default is to map NArgs over all arguments, which turns the integers
into reals."

NArgs::usage =
"NArgs[args] returns args with integers turned into reals.  Note that
NArgs is not quite the same as N: while it changes 1 to 1., it leaves
m[1] intact so that array indices remain integers."

Newline::usage =
"Newline is an option of WriteExpr.  It specifies a string to be printed
after each Fortran statement."

Optimize::usage =
"Optimize is an option of PrepareExpr.  With Optimize -> True, variables
are introduced for subexpressions which are used more than once."

Expensive::usage =
"Expensive is an option of PrepareExpr.  It specifies patterns of objects
whose evaluation is expensive in terms of CPU time and which should be
hoisted from inner do-loops if possible."

MinLeafCount::usage =
"MinLeafCount is an option of PrepareExpr and Abbreviate.  It specifies
the minimum LeafCount a common subexpression must have in order that a
variable is introduced for it."

DebugLines::usage =
"DebugLines is an option of PrepareExpr.  It specifies whether debugging
statements are written out for each variable.  If instead of True a
string is given, the debugging messages are prefixed by this string. 
The actual debugging statements are constructed from the items in
$DebugCmd."

FinalTouch::usage =
"FinalTouch is an option of PrepareExpr.  It specifies a function which
is applied to each final subexpression, such as will then be written out
to the Fortran file."

ResetNumbering::usage =
"ResetNumbering is an option of PrepareExpr.  It restarts numbering
of variable names for temporary and duplicate expressions at zero."

CodeExpr::usage =
"CodeExpr[vars, tmpvars, exprlist] is the output of PrepareExpr and
contains a list of expressions ready to be written to a Fortran file,
where vars are the original variables and tmpvars are temporary
variables introduced in order to shrink individual expressions to a size
small enough for Fortran."

DebugLine::usage =
"DebugLine[var] emits a debugging statement (print-out of variable var)
when written to a Fortran file with WriteExpr.  DebugLine[var, tag]
prefixes the debugging message with the string \"tag\"."

$DebugCmd::usage =
"$DebugCmd is a list of four strings used to construct debugging
statements in Fortran.  The first two strings provide commands
issued before and after the debugging statement, the second two
strings provide the opening and closing of the debugging statement.
The debugging statement has two arguments, a string (the variable
name) and a number (the variable value)."

SplitSums::usage =
"SplitSums[expr] splits expr into a list of expressions such that index
sums (marked by SumOver) always apply to the whole of each part. 
SplitSums[expr, wrap] applies wrap to the coefficients of the SumOver."

ToDoLoops::usage =
"ToDoLoops[list, ifunc] splits list into patches which must be summed
over the same set of indices.  ifunc is an optional argument, where
ifunc[expr] must return the indices occurring in expr."

DoLoop::usage =
"DoLoop[expr, ind] is a symbol introduced by ToDoLoops indicating that
expr is to be summed over the indices ind."

Dim::usage =
"Dim[i] returns the highest value the index i takes on, as determined
from the amplitude currently being processed."

NDim::usage =
"var[NDim[n]], used on the l.h.s. of a variable definition, indicates
that var is an array of dimension n."

MoveDepsRight::usage =
"MoveDepsRight[r1, ..., rn] shuffles variable definitions (var -> value)
among the lists of rules ri such that the definitions in each list do
not depend on definitions in ri further to the left.  For example,
MoveDepsRight[{a -> b}, {b -> 5}] produces {{}, {b -> 5, a -> b}}, i.e.
it moves a -> b to the right list because that depends on b."

MoveDepsLeft::usage =
"MoveDepsLeft[r1, ..., rn] shuffles variable definitions (var -> value)
among the lists of rules ri such that the definitions in each list do
not depend on definitions in ri further to the right.  For example,
MoveDepsLeft[{a -> b}, {b -> 5}] produces {{b -> 5, a -> b}, {}}, i.e.
it moves b -> 5 to the left list because that depends on b."

OnePassOrder::usage =
"OnePassOrder[r] orders a list of interdependent rules such that the
definition of each item (item -> ...) comes before its use in the
right-hand sides of other rules."

$OnePassDebug::usage =
"When OnePassOrder detects a recursion among the definitions of a list,
it deposits the offending rules in an internal format in $OnePassDebug
as debugging hints."

Tag::usage =
"Tag[t, expr] tags expr with t (possibly empty).  This tag is 
transparent to the functions MoveDepsLeft, MoveDepsRight, OnePassOrder."

SubroutineDecl::usage =
"SubroutineDecl[name] returns a string with the declaration of the
Fortran or C subroutine name.  SubroutineDecl[name[args], decl]
declares the subroutine name with arguments args, where decl is a
string with the declaration of args."

SubroutineEnd::usage =
"SubroutineEnd[] returns a string with the proper closing of a
subroutine in Fortran or C."

VarDecl::usage =
"VarDecl[v, t] returns a string with the declaration of v as variables
of type t in Fortran or C.  VarDecl[v1, t1, v2, t2, ...] does the same
for several variable types.  VarDecl[Common[b][v, t...]] puts the
variables inside common block b.  VarDecl[NameList[b][v, t]] generates
an array definition plus a name-index map.  Any other strings are
written out verbatim."

Common::usage =
"Common[b][v, t] is used inside VarDecl to indicate that the variables
v of type t are to be put inside the common block b."

NameMap::usage =
"NameMap[b][v, t] is used inside VarDecl to indicate that the variables
v of type t are to be put inside an array together with preprocessor
definitions to map variable names onto array indices."

NotEmpty::usage =
"NotEmpty[vars] is used inside VarDecl and outputs its contents only
if at least one variable list is not empty."

Extern::usage =
"Extern is a variable `type' recognized by VarDecl which is translated
into a declaration for an external function."

DoDecl::usage =
"DoDecl[v, m] returns two strings with the do/enddo declarations of a
Fortran loop over v from 1 to m.  DoDecl[v, {a, b}] returns the same for
a loop from a to b.  DoDecl[v] invokes Dim[v] to determine the upper
bound on v."

CallDecl::usage =
"CallDecl[names] returns a string with the invocations of the
subroutines names in Fortran, taking into account possible loops
indicated by DoLoop."

$SymbolPrefix::usage =
"$SymbolPrefix is a string prepended to all externally visible symbols
in the generated Fortran code to avoid symbol collisions."


(* symbols used in the Fortran code *)

Ctree::usage =
"Ctree[Fi] is the ith form factor (the coefficient of Fi) of the
tree-level amplitude."

Cloop::usage =
"Cloop[Fi] is the ith form factor (the coefficient of Fi) of the
one-loop amplitude."

SInvariant::usage =
"SInvariant[ki, kj] represents the s-type (invariant-mass type)
invariant formed from the momenta ki and kj, i.e. s_{ij} = (ki + kj)^2."

TInvariant::usage =
"TInvariant[ki, kj] represents the t-type (momentum-transfer type)
invariant formed from the momenta ki and kj, i.e. t_{ij} = (ki - kj)^2."

exp::usage =
"exp[x] is the exponential function in Fortran."

cI::usage =
"cI represents the imaginary unit in Fortran."

SplitChain::usage =
"SplitChain[w] prepares WeylChain[w] for numerical evaluation."

ChainHead::usage =
"ChainHead[om, n] constructs the header of a Weyl chain with n
vectors preceded by chirality projector om, i.e. of the form
<s1|om|v1...vn|s2>."

ChainV0::usage = ChainB0::usage =
ChainV1::usage = ChainB1::usage =
ChainV2::usage = ChainB2::usage =
ChainV3::usage = ChainB3::usage =
ChainV4::usage = ChainB4::usage =
ChainV5::usage = ChainB5::usage =
ChainV6::usage = ChainB6::usage =
"ChainXn[sL,epsL, vec..., epsR,sR] is the Fortran representation of
the spinor chain <sL|epsL|vec...|epsR|sR>.  The spinors sL and
sR are multiplied by the spinor metric if epsL,R = 1, respectively. 
ChainVn starts with a sigma, ChainBn with a sigma-bar."

MomEncoding::usage =
"MomEncoding[f, i] is the encoded version of momentum i with (integer)
prefactor f."


(* system variables *)

$Editor::usage =
"$Editor specifies the command line used for editing FORM code in a
detached window."

$EditorModal::usage =
"$Editor specifies the command line used for editing FORM code in a
non-detached (modal) window."

$FormCalc::usage =
"$FormCalc contains the version number of FormCalc."

$FormCalcDir::usage =
"$FormCalcDir is the directory from which FormCalc was loaded."

$FormCalcSrc::usage =
"$FormCalcSrc is the directory containing the FormCalc source files."

$FormCalcBin::usage =
"$FormCalcBin is the directory containing the FormCalc binary files."

$ReadForm::usage =
"$ReadForm contains the location of the ReadForm executable."

$FormCmd::usage =
"$FormCmd specifies the invocation of the FORM executable.  Arguments
are separated by a vertical bar (|)."

$DriversDir::usage =
"$DriversDir is the path where the driver programs for the generated
Fortran code are located."

$BlockSize::usage =
"$BlockSize is the maximum LeafCount a single Fortran statement written
out by WriteExpr may have.  Any expression with LeafCount > $BlockSize
will be chopped up before being written to the Fortran file."

$FileSize::usage =
"$FileSize gives the maximum LeafCount the expressions in a single
Fortran file may have.  If the expressions grow larger than $FileSize,
the file is split into several pieces."


Begin["`Private`"]

$FormCalc = 8.1

$FormCalcDir = DirectoryName[ File /.
  FileInformation[System`Private`FindFile[$Input]] ]

$FormCalcSrc = ToFileName[{$FormCalcDir, "FormCalc"}]

$FormCalcBin = ToFileName[{$FormCalcDir, $SystemID}]

$DriversDir = ToFileName[{$FormCalcDir, "drivers"}]

$ReadForm = ToFileName[$FormCalcBin, "ReadForm"];

Check[
  Install[$ReadForm],
  ReadForm::notcompiled = "The ReadForm executable `` could not be \
installed.  Did you run the compile script first?";
  Message[ReadForm::notcompiled, $ReadForm];
  Abort[] ]

$FormCmd = ToFileName[$FormCalcBin, "tform"] <>
  "|-w" <> ToString[Max[1, $ProcessorCount]]

FormCode[file_] := FormCode[file] = "##\n\n" <> ReadList[
  ToFileName[$FormCalcSrc, file],
  Record, RecordSeparators -> {} ] <> "\n##\n"

(*
FormCode[file_] := "#include " <> file <> "\n";
$FormCmd = $FormCmd <> "|-p|" <> $FormCalcSrc
*)

If[ StringMatchQ[$SystemID, "Windows*"],
  $EditorModal = "notepad ``";
  $Editor := "start /b " <> $EditorModal;
  Escape[s_] := StringReplace[s, " " -> "^ "],
(* else *)
  $EditorModal = "${VISUAL:-xterm -e nano} ``";
  $Editor := $EditorModal <> " &";
  Escape[s_] := "\"" <> s <> "\"" ]


$NumberMarks = False

Off[General::spell1, General::spell]


If[ $VersionNumber < 6,
  Needs["Algebra`Horner`"];

  (* actually load the Horner package so that the Off works: *)
  Algebra`Horner[1];
  Off[Algebra`Horner::fail];

  System`HornerForm = Algebra`Horner`Horner;

  Unprotect[StringMatchQ];
  StringMatchQ[s_, l_List] := !VectorQ[l, !StringMatchQ[s, #]&];
  Protect[StringMatchQ]
]


SetOptions[ToString, CharacterEncoding -> "ASCII"]

SetOptions[OpenWrite, CharacterEncoding -> "ASCII"]

FilterOptions = System`Utilities`FilterOptions


(* generic functions *)

ParseOpt[func_, opt___] :=
Block[ {names = First/@ Options[func]},
  Message[func::optx, #, func]&/@
    Complement[First/@ {opt}, names];
  names //. {opt} //. Options[func]
]


ToSymbol[x__] := ToExpression[ StringJoin[ToString/@ Flatten[{x}]] ]


_SymbolNumber = 0

next[stub_] := stub <> ToString[++SymbolNumber[stub]]

NewSymbol[stub_, 0] := ToExpression[next[stub]]

NewSymbol[stub_] :=
  Block[{sym = next[stub]}, ToExpression[sym] /; !NameQ[sym]]

NewSymbol[stub_] := NewSymbol[stub]


Attributes[ToArray] = {Listable}

ToArray[sym_Symbol] :=
Block[ {t = ToString[sym], p},
  p = Position[DigitQ/@ Characters[t], False][[-1, 1]];
  ToExpression[StringTake[t, p]] @ ToExpression[StringDrop[t, p]]
]

ToArray[other_] = other

ToArray[expr_, vars__] :=
Block[ {v = ToString/@ Flatten[{vars}], rules, t, p, h},
  rules = (
    t = ToString[#];
    p = Position[DigitQ/@ Characters[t], False][[-1, 1]];
    h = StringTake[t, p];
    If[ MemberQ[v, h],
      # -> ToExpression[h] @ ToExpression[StringDrop[t, p]],
      {} ]
  )&/@ Symbols[expr];
  expr /. Dispatch[Flatten[rules]]
]


Renumber[expr_, vars__] :=
Block[ {v = Flatten[{vars}], old, new},
  old = Union[Cases[expr, #[__], Infinity]]&/@ v;
  new = MapThread[Array, {v, Length/@ old}];
  expr /. Dispatch[Thread[Flatten[old] -> Flatten[new]]]
]


ExpandSums[x_, ___] := x /; FreeQ[x, SumOver]

ExpandSums[a_. IndexDelta[i_, j_] SumOver[i_, _], h___] :=
  ExpandSums[a /. i -> j, h]

ExpandSums[a_. s__SumOver, h_:Sum] :=
  h[a, Evaluate[Sequence@@ List@@@ {s}]]

ExpandSums[other_, h___] := ExpandSums[#, h]&/@ other


Kind[Tag[___, x_]] := Kind[x]

Kind[h_ -> _] := kind[h]

_Kind = Sequence[]

kind[h_[___], ___] := h

kind[h_, ___] := h


Alt[l_List] := Alt@@ Flatten[l]

Alt[s_] := s

Alt[s__] := Alternatives[s]


MaxDims[args__] :=
  listmax/@ Split[Union[Flatten[{args}]], Head[#1] === Head[#2]&]

listmax[{s__Symbol}] := s

listmax[{s__String}] := s

listmax[l_] := l[[1, 0]]@@ MapThread[Max, List@@@ l]


Attributes[Keep] = {HoldFirst}

Keep[lhs_ = rhs_] := lhs = Keep[rhs, Block[{lhs}, ToString[lhs]]]

Keep[lhs_ = rhs_, other__] := lhs = Keep[rhs, other]

Keep[cmd_, name_String, prefix_String:$KeepDir] :=
Block[ {file = ChkExist[prefix, name <> ".m"]},
  If[ FileType[file] === File,
    Print["loading ", file];
    (RegisterAbbr[#1]; RegisterSubexpr[#2]; #3)&@@ Get[file],
  (* else *)
    (Put[{Abbr[], Subexpr[], #}, file]; #)& @ cmd ]
]

$KeepDir = "keep/"


Attributes[ToForm] = {Listable}

ToForm[x_String] := x

ToForm[x_] := ToString[x, InputForm]


ToSeq[li_List, opt___] :=
  StringTake[ToString[li, InputForm, opt], {2, -2}]

ToSeq[x_, ___] := ToForm[x]


ToBool[True] = "1"

ToBool[___] = "0"


ToCode[x_String] := x

ToCode[x_List] := StringTake[ToString[x, CodeForm], {6, -2}]

ToCode[x_] := ToString[x, CodeForm]

ToFortran[x_] := ToCode[x]	(* legacy *)


ToCat[n_, {}] := Table[{}, {n}]

ToCat[_, li_] := Flatten/@ Transpose[li]


Symbols[expr_] := Union[Cases[expr, _Symbol, {-1}]]


FromPlus[h_, p_Plus] := h@@ p

FromPlus[_, other_] = other


nobrk = Alternatives[]

DotSimplify[f_, _][expr_] := f[expr] /; FreeQ[expr, nobrk]

DotSimplify[_, f_][expr_] := f[expr]


(* should be better but locks up sometimes:
TermCollect[x_] := x //. a_ b_. + a_ c_ -> a (b + c) *)

TermCollect[x_] := x //. a_ b_ + a_ c_ :> a (b + c)


DenCancel[x_] := x /. {
  Den[a_, b_]^n_. (a_ - b_) :> Den[a, b]^(n - 1),
  Den[a_, b_]^n_. (b_ - a_) :> -Den[a, b]^(n - 1) }


numadd[term_] := numadd[
  Numerator[term],
  Denominator[term] //Simplify ]

numadd[n_, x_?NumberQ d_] := numer[d] += n/x

numadd[n_, d_] := numer[d] += n

DenCollect[p_Plus, wrap_:Identity] :=
Block[ {numer},
  _numer = 0;
  numadd/@ p;
  _numer =.;
  Plus@@ (wrap[#2]/#1[[1, 1]] &)@@@ DownValues[numer]
]

DenCollect[x_, wrap_:Identity] := wrap[x] /; FreeQ[x, Plus]

DenCollect[x_, wrap___] := DenCollect[#, wrap]&/@ x


Pool[expr_, wrap_:Identity] := expr /.
  p_Plus :> ploos[wrap]@@ Cases[p, _Times] +
    DeleteCases[p, _Times] /; LeafCount[p] > 10

ploos[wrap_][a_, r__] :=
Block[ {pos = 0, lcmin, lcmax, lc, ov, ovmax, ovpos, ploos},
  lcmin = Floor[LeafCount[a]/3];
  lcmax = -Infinity;
  Scan[ (
    ++pos;
    lc = LeafCount[ov = Intersection[a, #]];
    If[ lc > lcmax,
      lcmax = lc; ovpos = pos; ovmax = ov;
      If[ lc > lcmin, Return[] ] ] )&, {r} ];
  If[ lcmax < 5,
    a + ploos[wrap][r],
  (* else *)
    ovmax wrap[a/ovmax + {r}[[ovpos]]/ovmax] +
      Drop[ploos[wrap][r], {ovpos}] ]
]

ploos[_][other___] := Plus[other]


Attributes[usq] = {Orderless}

usq/: usq[a__][ik__] + usq[b__][ik__] := usq[a, b][ik]

usq/: usq[a__][h_, i_, i_] + usq[i_][h_, k_, k_] := usq[a, k][h, i, i]


usqplus[simp_][p__] := simp[Plus[p]] /;
  Length[Position[{p}, usq[__][__], 2, 1]] > 0

usqplus[_][p__] := Plus[p]


ApplyUnitarity[expr_, U_, dim_Integer, simp_:FullSimplify] :=
Block[ {ux, uexpr, pos},
  Evaluate[ux@@ Range[dim]] = KroneckerDelta[##2] &;
  ux[a__] := With[ {cpl = usq@@ Complement[Range[dim], {a}]},
    ux[a] = KroneckerDelta[##2] - cpl[##] &
  ] /; Length[{a}] > dim/2;
  ux[other__] := usq[other];

  uexpr = expr //. {
    (u:U[i_, j_])^n_. (uc:Conjugate[U[k_, j_]])^nc_. :>
      (usq[j][1, i, k]^# u^(n - #) uc^(nc - #) &) @ Min[n, nc],
    (u:U[j_, k_])^n_. (uc:Conjugate[U[j_, i_]])^nc_. :>
      (usq[j][2, i, k]^# u^(n - #) uc^(nc - #) &) @ Min[n, nc]
  } /. usq -> ux;

  uexpr /. Plus -> usqplus[simp] /. usq -> ux /. {
    usq[a__][1, i_, k_] -> Plus@@ (U[i, #] Conjugate[U[k, #]] &)/@ {a},
    usq[a__][2, i_, k_] -> Plus@@ (U[#, k] Conjugate[U[#, i]] &)/@ {a} }
]


Attributes[osfun] = {HoldAll}

osfun[w___][Map] := With[ {simp = Unique[simp]},
  simp = Which[w, True, simp/@ #] & ]

osfun[w___][] := Which[w, True, #] &

osfun[w___][f_] := Which[w, True, f[#]] &

osfun[w___][n_, f_, r___] := osfun[w, LeafCount[#] < n, f[#]][r]

OnSize[args__] := osfun[][args]


DiagramType[a_FeynAmp] := Exponent[
  a[[3]] /. _FeynAmpDenominator -> 1 /. _PropagatorDenominator -> pd,
  pd ]


FermionicQ[a_FeynAmp] := !FreeQ[a[[3]], FermionChain | MatrixTrace]


Attributes[IndexDiff] = Attributes[IndexDelta] = {Orderless}

IndexDelta[i_, i_] = 1

IndexDelta[i_Integer, _Integer] = 0

IndexDiff[i_, i_] = 0

IndexDiff[i_Integer, _Integer] = 1


IndexIf[] = 0

IndexIf[a_] := a

IndexIf[True, a_, ___] := a

IndexIf[False, _, b___] := IndexIf[b]

IndexIf[a__] := IndexIf[a, 0] /; EvenQ[Length[{a}]]

IndexIf[cond_, 0, a_] := IndexIf[!cond, a, 0]

IndexIf[cond1_, IndexIf[cond2_, a_, 0], 0] := IndexIf[cond1 && cond2, a, 0]

IndexIf[a___, i_IndexIf] := IndexIf[a, Sequence@@ i]


MapIf[foo_, i_IndexIf] := MapAt[foo, i,
  List/@ Append[Range[2, # - 1, 2], #]]& @ Length[i]

MapIf[foo_, other_] := foo/@ other


Off[Optional::opdef]

ToIndexIf[expr_, s_String] := ToIndexIf[expr, Alt[Names[s]]]

ToIndexIf[expr_, patt_:_] :=
  Fold[ singleIf, expr, Union @ Cases[expr,
    (IndexDelta | IndexDiff)[i:patt..] -> {i}, Infinity] ]

singleIf[expr_, {i__}] := expr /.
  { IndexDelta[i] -> suck[1, 1],
    IndexDiff[i] -> suck[2, 1] } /.
  a_. suck[1, x_] + a_. suck[2, y_] -> a IndexIf[Equal[i], x, y] /.
  suck[h_, x_] :> IndexIf[{Equal, Unequal}[[h]][i], x]

suck/: r_ suck[h_, x_] := suck[h, r x] /; FreeQ[r, Mat (*| SumOver*)]

suck/: suck[h_, x_] + suck[h_, y_] := suck[h, x + y]

(* suck/: suck[h_, x_] + y_ := suck[3 - h, y] /; x + y == 0 *)


(* preparations for FORM *)

ExtWF = {e, ec, z, zc, eT, eTc}

ConjWF = ({#1 -> #2, #2 -> #1}&)@@@ Partition[ExtWF, 2] //Flatten

KinVecs = {eta, e, ec, z, zc, k}

KinTens = {eT, eTc}

KinObjs = Join[KinVecs, KinTens]
	(* note: KinObjs determines the order in which vectors
	   and tensors appear in FORM and hence ultimately in
	   functions like Pair and Eps *)


Attributes[KinFunc] = {HoldAll}

KinFunc[args__] := Function[Evaluate[KinObjs], args]


FromFormRules = Outer[ ToSymbol["FormCalc`", #2, #1] -> #2[#1] &,
  Range[8], KinObjs ]

FormKins = Apply[#1&, FromFormRules, {2}]

FromFormRules = Flatten[FromFormRules]


MomThread[i_Index, foo_] := foo[i]

MomThread[p_Symbol, foo_] := foo[p]

MomThread[p_, foo_] := Replace[MomReduce[p], k_Symbol :> foo[k], {-1}]


MomReduce[p_Plus] := Fewest[p, p + MomSum, p - MomSum]

MomReduce[p_] := p


Fewest[a_, b_, r___] := Fewest[a, r] /; Length[a] <= Length[b]

Fewest[_, b__] := Fewest[b]

Fewest[a_] := a


fvec[p_] := (vecs = {vecs, Symbols[p]}; MomReduce[p])

fvec[p_, mu_] := (vecs = {vecs, Symbols[p]}; MomThread[p, #[mu]&])


iname[type_, n_] := iname[type, n] =
  ToSymbol[StringTake[ToString[type], 3], n]


Attributes[idelta] = {Orderless}

idelta[c1:Index[Colour, _], c2_] := SUNT[c1, c2]

idelta[g1:Index[Gluon, _], g2_] := 2 SUNT[g1, g2, 0, 0]

idelta[x__] := IndexDelta[x]


ieps[c__] := Block[{eps = SUNEps[c]}, eps /; !FreeQ[eps, Index[Colour, _]]]

ieps[x__] := IndexEps[x]


psum/: psum[n_][x__]^k_ := psum[n k][x]

_isum[expr_, {i_, f_:1, t_}] := (t - f + 1) expr /; FreeQ[expr, i]

isum[opt___][expr_, {i_, r__}] :=
Block[ {dummy = NewSymbol[ToString[i]]},
  indices = {indices, dummy};
  (expr /. i -> dummy) SumOver[dummy, r, opt]
]

_isum[expr_, i_] :=
Block[ {dummy = NewSymbol[ToString[i]]},
  indices = {indices, dummy};
  ranges = {dummy -> dummy, ranges};
  (expr /. i -> dummy)
]


sumover[i_, 1, r_] := sumover[i, r]

sumover[i:Index[Colour | Gluon, _], r__] := SUNSum[i, r]

sumover[other__] := SumOver[other]

SUNSum[_, _, External] = 1


KinFunc[
  pol[_, k, mu:Index[EpsilonScalar, _]] = z[mu];
  polc[_, k, mu:Index[EpsilonScalar, _]] = zc[mu];
  pol[_, k, mu_] = e[mu];
  polc[_, k, mu_] = ec[mu];
  pol[_, k, mu__] = eT[mu];
  polc[_, k, mu__] = eTc[mu] ]@@@ FormKins

Conjugate[pol] ^= polc


Attributes[scalar] = {Orderless}

scalar[0, _] = 0

scalar[a_Symbol, b_Symbol] := a . b

scalar[a_, p:_[__]] := MomThread[p, scalar[a, #]&]


prop[0, m_Symbol, d___] := -m^(-2 d)

prop[p_, m__] := prop[-p, m] /; !FreeQ[p, -q1]

prop[p_, m_, d___] := Den[p, m^2, d]


loop[a___, Den[p_, m_], b___, Den[p_, x_ m_], c___] :=
  loop[a, Den[p, m] - Den[p, x m], b, c]/(m (1 - x))

(*
loop[a___, Den[p_, m_, d1___], b___, Den[p_, m_, d2___], c___] :=
  loop[a, Den[p, m, 1 d1 + 1 d2], b, c]
*)

loop[a___, Den[p_, m_, d_], b___] := loop[a, ##, b]&@@ Table[Den[p, m], {d}]

loop[a___, d1_Den - d2_Den, b___] := loop[a, d1, b] - loop[a, d2, b]

loop[d__] := I Pi^2 intM[d]


noncomm[p_Plus] := noncomm/@ p

noncomm[g_] := g

noncomm[g__] := NonCommutativeMultiply[g]


Neglect[m_] := m

FormPatt[_[_, m (* verbatim m, matches rhs in above Neglect def *)]] = {}

FormPatt[_?NumberQ, _] = {}

FormPatt[_[_[_[lhs_]], rhs_]] := FormPatt[lhs, rhs]

FormPatt[lhs_Alternatives, rhs_] := FormPatt[#, rhs]&/@ List@@ lhs

FormPatt[lhs_, rhs_] :=
Block[ {c = 96, newlhs, newrhs = rhs},
  newlhs = lhs /.
    {Blank -> FormArg, BlankSequence | BlankNullSequence -> FormArgs} /.
    Pattern -> ((newrhs = newrhs /. #1 -> #2; #2)&);
  (newlhs /. patt[x_] :> x <> "?") -> (newrhs /. patt[x_] -> x)
]

FormArg[h_:Identity] := h[patt["ARG" <> FromCharacterCode[++c]]]

FormArgs[h_:Identity] := h["?" <> FromCharacterCode[++c]]


OrdSq[r:_[_, rhs_]] := {{}, r} /; VectorQ[lhs, FreeQ[rhs, #]&]

OrdSq[r_] := {r, {}}

SortSq[dv_] :=
Block[ {lhs = #[[1, 1, 1]]&/@ dv},
  Flatten[Transpose[OrdSq/@ dv]]
]


Attributes[Inv] = {Listable}

Inv[i_, j_] :=
Block[ {s = signs[[i]] signs[[j]], ki = Moms[[i]], kj = Moms[[j]], inv},
  inv = If[ Legs === 3, dot[#, #]& @ Moms[[3]], Invariant[s, i, j] ];
  dot[ki, kj] = s/2 inv - s/2 dot[ki, ki] - s/2 dot[kj, kj];
  inv
]


Invariant[1, 1, 2] = S

Invariant[-1, 1, 3] = T

Invariant[-1, 2, 3] = U

Invariant[s_, i_, j_] := Invariant[s, i, j] =
  ToSymbol["FormCalc`", FromCharacterCode[(167 - s)/2], i, j]


OtherProd[{k1___, k_}, {s1___, s_}, zero_] :=
  MapThread[
    (dot[#1, k] = -s (Plus@@ ({s1} dot[#1, {k1}]) - #2 zero))&,
    {{k1}, {s1}} ]


NLegs[n_Integer] := n

NLegs[f_Integer -> t_Integer] := f + t

NLegs[f_ -> t_] := Length[f] + Length[t]


Signs[_Integer] = 0

Signs[f_Integer -> t_Integer] := Join[Table[1, {f}], Table[-1, {t}]]

Signs[f_ -> t_] := Join[1&/@ f, -1&/@ t]


Masses[{f___} -> {t___}] := #3&@@@ {f, t}

_Masses = {}


FormMom[{f___} -> {t___}] :=
  Cases[Thread[(#2&@@@ {f, t}) -> Moms], _[_FourMomentum, _]]

_FormMom = {}


Attributes[GetProc] = {Listable}

GetProc[Amp[proc_][___]] := proc

GetProc[FeynAmpList[info___][___]] := Process /. {info}

GetProc[proc_Rule] := proc


General::incomp = "Warning: incompatible processes ``."

Attributes[ChkProc] = {HoldRest}

ChkProc[amp_, h_, fail___] := (
  If[ !SameQ@@ Apply[#3&, #, {3}], Message[h::incomp, Union[#]]; fail ];
  #[[1]] )& @ GetProc[amp]


(* global variables set here:
   CurrentProc, CurrentOptions,
   FormProcs, FormSymbols,
   Legs, Moms, MomSubst, MomSum, InvSum, PairRules, LastAmps *)

CurrentProc = CurrentOptions = {}

Options[DeclareProcess] = {
  OnShell -> True,
  Invariants -> True,
  Transverse -> True,
  Normalized -> True,
  InvSimplify -> True,
  Antisymmetrize -> True }

Attributes[DeclareProcess] = {Listable}

DeclareProcess::incomp =
"Calculation of incompatible process(es) attempted.  If you want to \
calculate a new process, run ClearProcess[] first."

DeclareProcess::syntax =
"Wrong syntax: DeclareProcess expects FeynAmpList objects as arguments."

DeclareProcess[fal:FeynAmpList[__][___].., opt___Rule] := (
  DeclP[#, ParseOpt[DeclareProcess, opt], CurrentOptions];
  Level[LastAmps, {2}, FormAmp[#]]
)& @ ChkProc[LastAmps = {fal}, DeclareProcess, Abort[]]

DeclareProcess[l__List, opt___Rule] := DeclareProcess@@ Flatten[{l, opt}]

_DeclareProcess := (Message[DeclareProcess::syntax]; Abort[])


DeclP[_, opt_, opt_] = Null

DeclP[proc_, opt:{onshell_, inv_, transv_, norm_, invsimp_, antisymm_}, ___] :=
Block[ {signs, masses, kins, dot, n, neglect, square, id,
invproc = {}, invs = {}, kikj = {}, eiki = {}, eiei = {}},
  CurrentProc = proc;
  CurrentOptions = opt;
  MomSum = InvSum = 0;
  FormProcs = {};
  FormSymbols = {
    neglect = FormPatt/@ DownValues[Neglect],
    square = FormPatt/@ SortSq[DownValues[Square]] };

  signs = Signs[CurrentProc];
  masses = Masses[CurrentProc]^2 /. Index -> iname;
  Legs = NLegs[CurrentProc];
  kins = Take[FormKins, Legs];

  Moms = KinFunc[k]@@@ kins;
  MomSubst = {};
  FormVectors =
    Level[{q1, MapThread[List, KinFunc[Evaluate[KinVecs]]@@@ kins]}, {-1}];
  FormTensors =
    Level[MapThread[List, KinFunc[Evaluate[KinTens]]@@@ kins], {-1}];

  Attributes[dot] = {Orderless, Listable};

  kikj = dot[Moms, Moms];
  If[ onshell && Length[kins] === Length[masses],
    kikj = MapThread[Set, {kikj, masses}] ];

  If[ transv, eiki = KinFunc[{e.k -> 0, ec.k -> 0}]@@@ kins ];

  If[ norm, eiei = KinFunc[e.ec -> -1]@@@ kins ];

  Switch[ Length[masses],
    0 | 1,
      Null,
    2,
      (*dot[ Moms[[2]], Moms[[2]] ] =.;*)
      MomSubst = {s_Spinor -> s, Moms[[2]] -> Moms[[1]]};
      eiki = eiki /. MomSubst;
      FormSymbols = {FormSymbols, Spinor[]},
    _,
      If[ inv,
        invs = Flatten[Array[Inv[Range[# - 1], #]&, Legs - 2, 2]];
        Scan[(RealQ[#] = True)&, invs];
        InvSum = kikj[[-1]] + Plus@@ ((Legs - 3) Drop[kikj, -1]);

	(* The number of invariants is ninv = (Legs - 1)(Legs - 2)/2 in
	   total and ndot = (Legs - 2) in a dot product (of the terms in
	   pi.pLegs = Sum[pi.pj, {j, Legs - 1}], pi.pi is a mass^2).
	   Thus, the number of invariants in pi.pLegs can be reduced by
	   using the Mandelstam relation only if ndot > ninv/2, that is
	   Legs < 5.  The case Legs = 3 is handled specially already by
	   Inv, so only Legs = 4 remains. *)
        OtherProd[ Moms, signs,
          If[Legs === 4, Distribute[(Plus@@ invs - InvSum)/2], 0] ];

        If[ invsimp && onshell && Legs =!= 3,
          invproc = ToForm[MapIndexed[
            { "id `foo'(ARG?) = `foo'(ARG, ARG*replace_(",
              #1, ", ", InvSum - Plus@@ Drop[invs, #2],
              "));\n#call Fewest(`foo')\n" }&, invs ]] ];

	  (* not used anywhere in FormCalc, but useful for debugging: *)
        InvSum = Plus@@ invs - InvSum
      ];

      n = signs Moms;
      MomSum = Plus@@ n;
      FormProcs = {FormProcs,
        "#define MomSum \"", MomSum, "\"\n",
        Array[{"#define ", Moms[[#]], " \"",
          -signs[[#]] Plus@@ Drop[n, {#}], "\"\n"}&, Legs]};
  ];

  kikj = ReleaseHold[DownValues[dot] /. dot -> Dot];
  FormSymbols = Flatten[{FormSymbols, Last/@ kikj}];

  FormProcs = "\n\
#define Legs \"" <> ToString[Legs] <> "\"\n\
#define Invariants \"" <> ToSeq[Cases[invs, _Symbol]] <> "\"\n\
#define OnShell \"" <> ToBool[onshell] <> "\"\n\
#define Antisymmetrize \"" <> ToBool[antisymm] <> "\"\n\n" <>
    ToForm[FormProcs] <> "\n\
#procedure Neglect\n" <>
    FormId[neglect] <> "\
#endprocedure\n\n\
#procedure Square\n" <>
    FormSq[square] <> "\
#endprocedure\n\n\
#procedure eiei\n" <>
    FormId[eiei] <> "\
#endprocedure\n\n\
#procedure eiki\n" <>
    FormId[eiki] <>
    FormId[Cases[kikj, _[_, 0]]] <> "\
#endprocedure\n\n\
#procedure kikj\n" <>
    FormId[kikj] <> "\
#endprocedure\n\n\
#procedure InvSimplify(foo)\n" <>
    invproc <> "\
#endprocedure\n\n";

	(* not used anywhere in FormCalc, but useful for debugging: *)
  PairRules = Flatten[{eiei, eiki, kikj}] /. Dot -> Pair /. FromFormRules;
]


Options[ToFeynAmp] = {
  Process -> Automatic
}

ToFeynAmp[amps___, opt___Rule] :=
Block[ {proc},
  {proc} = ParseOpt[ToFeynAmp, opt] /. Automatic :>
    Level[ Select[FromFormRules, !FreeQ[ {amps}, #[[1]] ]&],
      {3}, Max];
  MapIndexed[ToAmp, FeynAmpList[Process -> proc][amps] /.
    Reverse/@ FromFormRules /.
    PolarizationVector | PolarizationTensor -> pol /.
    Den -> PropagatorDenominator /. {
    g_MetricTensor :> ToIndex[Lorentz]/@ g,
    v:(Alt[First/@ FromFormRules][__]) :> ToIndex[Lorentz]/@ v,
    FourVector[v_, mu_] :> FourVector[v, ToIndex[Lorentz][mu]],
    e_LeviCivita :> ToIndex[Lorentz]/@ e,
    SUNT[g__, c1_, c2_] :> SUNT[ToIndex[Gluon][g], ToIndex[Colour][c1, c2]],
    SUNF[g__] :> SUNF[ToIndex[Gluon][g]]
  }]
]


ToAmp[FeynAmp[g_, q_, amp_, r___], _] := FeynAmp[g, q, ToDen[amp], r]

ToAmp[amp_, {n_}] := FeynAmp[GraphID[Generic == n], Integral[], ToDen[amp]]


ToDen[amp_] := amp /.
  d:PropagatorDenominator[p_, __] :> FeynAmpDenominator[d] /; !FreeQ[p, q1] /.
  HoldPattern[Times[d__FeynAmpDenominator]] :> Join[d]


_ToIndex[i_Index] := i

_ToIndex[v_FourVector] := v

(_ToIndex[#] := #)& @@@ FromFormRules

ToIndex[t_][n_] := Index[t, n]

t_ToIndex[i__] := Sequence@@ t/@ {i}


LevelSelect[_][id_, _, amp_] :=
  {# -> TrivialSums[ampden[GenNames[amp]]], #}& @
    AmpName[Select[id, FreeQ[#, Number]&]]

LevelSelect[Automatic][r__, gm_ -> ins_] := LevelSelect[
  Which[
    !FreeQ[ins, Particles], Particles,
    !FreeQ[ins, Classes], Classes,
    True, Generic ]
][r, gm -> ins]

LevelSelect[Generic][id_, _, gen_, ___] :=
  {# -> ampden[GenNames[gen]], #}& @
    AmpName[Select[id, FreeQ[#, Classes | Particles | Number]&]]

LevelSelect[lev_][id_, _, gen_, gm_ -> ins_] :=
Block[ {old, new, amp, name = AmpName[id], pc},
  _pc = 0;
  old = TrivialSums/@ Cases[{ins}, Insertions[lev][r__] :> r, Infinity];
  old = ampden[GetDen[gen], gm]/@ old;
  new = Thread[Flatten[ReduceIns[gm, Transpose[old]]], Rule];
  amp = gen /. new[[1]] /.
    {d_Den :> (d /. small[m_] -> m), _small -> 0};
  { name -> amp,
    If[ Length[ new[[2]] ] === 0,
      Length[old] name,
    (* else *)
      "i" <> name -> name Plus@@ (Level[#, {2}, "replace_"]&)/@
        Transpose[ Thread/@ new[[2]] ] ] }
]


AmpDen[amp_] := amp GetDen[amp]

AmpDen[1, _] = Identity

AmpDen[den_, gm_] = Block[ {ins = #1},
  ins[[-1]] *= den /. Thread[gm -> #1];
  ins ]&


GetDen[amp_] := Times@@ ExtDen/@
  Union[Cases[amp, d_intM /; Length[d] >= opp, Infinity]]


DenList = {}

Attributes[DenMatch] = {Orderless}

ppatt = {v1_, v2_, v3_, v4_, v5_, v6_, v7_, v8_}

DenExtend[fad_] := DenExtend[fad, DenCases@@ Transpose @
  MapThread[{ReplacePart[##, 1], #1[[1]], #2[[1]]}&,
    {List@@ fad, Take[ppatt, Length[fad]]}]]

DenExtend[fad_, {}] := (AppendTo[DenList, DenMatch@@ fad]; 1)

DenExtend[fad_, {extM[___, q1]}] = 1

DenExtend[fad_, {e_}] := Prepend[e, fad]


DenCases[{d__}, p_, v_] := Cases[DenList,
  x:DenMatch[d, r___] :> Block[ {q},
    extM[ intM@@ x, r, q[[1]] ] /;
      SameQ@@ (q = MomReduce/@ (q1 - p + v)) ||
      SameQ@@ (q = MomReduce/@ (q1 - p - v)) ],
  {1}, 1]


ReduceIns[{g_, rg___},
  {{ins_ /; (* LeafCount[ins] < 10 && *)
		(* Commuting functions in FORM are only commuting in
		   nonnegative powers.  The following FreeQ makes sure
		   that all expressions containing such are always
		   inserted *after* the kinematical simplification so
		   that FORM won't screw up the spinor chains. *)
            FreeQ[ins, _[__]^_?Negative], r___} /; ins === r, rins___}] :=
Block[ {smallg},
  smallg = If[ Neglect[ins] === 0, small[ins], ins ];
  { (g -> smallg) -> Sequence[], ReduceIns[{rg}, {rins}] }
]

ReduceIns[{g_, rg___}, {ins_List, rins___}] :=
Block[ {instype = InsType[g], newg, smallg},
  newg = ToSymbol[
    If[FreeQ[ins, DenyHide], "Form`p", "FormCalc`c"],
    instype, ++pc[instype] ];
  smallg = If[ Union[Neglect/@ ins] === {0}, small[newg], newg ];
  { (g -> smallg) -> (newg -> inssym/@ ins),
    ReduceIns[{rg}, Replace[{rins}, {ins -> smallg, -ins -> -smallg}, 1]] }
]

ReduceIns[{g_, rg___}, {newg_, rins___}] :=
  { (g -> newg) -> Sequence[], ReduceIns[{rg}, {rins}] }

ReduceIns[{}, {}] = {}


InsType[_Mass] = "M"

InsType[RelativeCF] = "R"

InsType[_GaugeXi] = "X"

InsType[_] = "G"


AmpName[g_] :=
Block[ {t, c},
  t = StringJoin@@ ({StringTake[ToString[#1], 1], ToString[#2]}&)@@@ g;
  If[ (c = ++uniq[t]) === 1, t, t <> "v" <> ToString[c] ]
]


TrivialSums[ins_ -> _] := TrivialSums[ins]

TrivialSums[ins_] := ins /; FreeQ[ins, SumOver]

TrivialSums[ins_] :=
Block[ {test = ins /. _SumOver -> 1},
  ins /. SumOver -> CarryOut
]

CarryOut[i_, n_, r___] := (
  ranges = {i -> i == n, ranges};
  sumover[i, n, r]
) /; !FreeQ[test, i]

CarryOut[_, v_, External] := Sqrt[v]

CarryOut[_, v_, ___] := v


Attributes[OffShell] = {Listable}

OffShell[fal:FeynAmpList[__][___].., extmass__Rule] :=
Block[ {subst},
  setsubst@@@ {extmass};
  _subst = Identity;
  Sequence@@ ({fal} /. (Process -> p_) :>
    Block[{c = 0}, Process -> Map[MapAt[subst[++c], #, 3]&, p, {2}]])
]

setsubst[n_, f_Function] := subst[n] = f

setsubst[n_, m_] := subst[n] = m &


Combine[fal:FeynAmpList[__][___]..] := (
  ChkProc[#, Combine];
  Level[#, {1}, #[[1, 0]]] )& @ {fal}

Combine[amp:Amp[_][___]..] :=
Block[ {comp},
  ChkProc[#, Combine];
  _comp = 0;
  Map[Component, #, {2}];
  _comp =.;
  #[[1, 0]]@@ Cases[DownValues[comp], _[_[_[x___]], r_] -> r x]
]& @ {amp}

Component[r_ s__SumOver] := comp[s] += r

Component[other_] := comp[] += other


m_MultiplyDiagrams[fal:FeynAmpList[__][___]] := m/@ fal

m_MultiplyDiagrams[fal:FormAmp[_][___]] := m/@ fal

MultiplyDiagrams[f_][FeynAmp[id_, q_, gen_, gm_ -> ins_]] :=
  FeynAmp[id, q, gen, gm -> MultiplyIns[f, gen, gm]/@ ins]

MultiplyDiagrams[f_][FeynAmp[id_, q_, amp_]] :=
  FeynAmp[id, q, f[amp] amp]

i_MultiplyIns[ins_ -> more_] := i[ins] -> i/@ more

MultiplyIns[f_, gen_, gm_][{r__, fac_}] :=
  {r, f[gen, Thread[gm -> {r, fac}]] fac}


TagDiagrams[diags_, tag_:Diagram] :=
Block[ {diag = 0},
  MultiplyDiagrams[tag[++diag]&][diags]
]


TagCollect[amp:Amp[_][___], tag__] := TagCollect[#, tag]&/@ amp

TagCollect[other_, tag_, foo_:Identity] :=
  Collect[other, tag, foobar] /.
    tag foobar[x_] :> foo[tag x] /.
    foobar -> Identity


IndexDim[Lorentz, i_] := i -> i

IndexDim[Lorentz4, i_] := i -> i == 4

IndexDim[EpsilonScalar, i_] := i -> i == Dminus4

_IndexDim = {}


helM[_, _Hel] = 1


(* the main function CalcFeynAmp *)

SUNObjs = SUNSum | SUNT | SUNTSum | SUNF | SUNEps

DenyNoExp = Level[{ga, Spinor, Den, A0i, B0i, C0i, D0i, E0i, F0i,
  SumOver, PowerOf, SUNObjs, "replace_"}, {-1}]

DenyHide = Level[{SumOver, PowerOf, IndexDelta, IndexEps, SUNObjs},
  {-1}, Alternatives]

FinalFormRules = {
  (f:_[__])^(n_?Negative) -> powM[f, n],
  x_^n_ :> powM[x, n] /; !IntegerQ[n],
  Complex[a_, b_] -> a + "i_" b }


Attributes[CalcFeynAmp] = {Listable}

Options[CalcFeynAmp] = {
  CalcLevel -> Automatic,
  Dimension -> D,
  MomElim -> Automatic,
  DotExpand -> False,
  NoCostly -> False,
  FermionChains -> Weyl,
  FermionOrder -> Automatic,
  Evanescent -> False,
  InsertionPolicy -> Default,
  SortDen -> True,
  CombineDen -> Automatic,
  PaVeReduce -> False,
  SimplifyQ2 -> True,
  OPP -> 100,
  OPPQSlash -> False,
  Gamma5Test -> False,
  Gamma5ToEps -> False,
  NoExpand -> {},
  NoBracket -> {},
  MomRules -> {},
  PreFunction -> Identity,
  PostFunction -> Identity,
  FileTag -> "amp",
  EditCode -> False,
  RetainFile -> False }

CalcFeynAmp::ddim = "Warning: `` \
implies the use of Fierz identities which are in general not \
applicable in D dimensions.  You know what you are doing."

CalcFeynAmp[FormAmp[proc_][amp___], opt___Rule] :=
Block[ {lev, dim, momelim, dotexp, nocost, fchain, forder, evanes,
inspol, sortden, combden, pavered, simpq2, opp, oppqsl, g5test,
noexp, nobrk, momrulz, pre, post, edit, retain,
uniq, vecs, kc = 0, hd, ic, inssym, mmains,
indices = {}, ranges = {}, indsym, haveferm,
intmax, extmax = 0, ampden, vars, hh, amps, res, traces = 0},

  {lev, dim, momelim, dotexp, nocost, fchain, forder, evanes,
    inspol, sortden, combden, pavered, simpq2, opp, oppqsl,
    g5test, g5eps, noexp, nobrk, momrulz,
    pre, post, tag, edit, retain} = ParseOpt[CalcFeynAmp, opt];

  opp = opp /. {True -> 1, Rational -> -1, False -> 100};

  If[ dim === 0,
    If[ Length[#] > 0, Message[CalcFeynAmp::ddim, #] ]&[{
      If[ forder =!= None, FermionOrder -> forder, {} ],
      If[ fchain === Weyl, FermionChains -> fchain, {} ]
    } //Flatten] ];

(*
  NoExp[p_] := NoExp[p] = NewSymbol["noexp"];
  NoExpandRule = (p_Plus :> NoExp[p] /; !FreeQ[p, #1] && FreeQ[p, #2])&[
    Alt[noexp],
    Alt[{FormVectors, FormTensors, DenyNoExp}] ];
*)

  NoExpandRule = (p_Plus :> (p /. Plus -> addM) /;
    !FreeQ[p, #1] && FreeQ[p, #2])&[
      Alt[noexp],
      Alt[{FormVectors, FormTensors, DenyNoExp}] ];

  If[ NumberQ[inspol],
    _ic = 0;
    inssym[ins_] := ins /; LeafCount[ins] < inspol || !FreeQ[ins, DenyHide];
    inssym[ins_] := inssym[ins] =
      ToSymbol["FormCalc`i", instype, ++ic[instype]],
  (* else *)
    inssym = Identity ];

  _uniq = 0;
  vecs = {};

  amps = pre/@ {amp} /.
    { g:G[_][_][__][_] -> g,
      IndexDelta -> idelta,
      IndexEps -> ieps,
      IndexSum -> psum[1] } //.
    psum[n_] -> (Product[Fold[isum[Renumber], #1, {##2}], {n}]&) /.
    FormMom[proc] /.
    FourMomentum[Internal, 1] -> q1 /.
    PolarizationVector | PolarizationTensor -> pol /.
    Prepend[momrulz,
     (DiracSpinor | MajoranaSpinor)[s_. p_Symbol, m_] :>
        Spinor[p, Neglect[m], s]];
  amps = amps /. {
    MetricTensor[mu_, _]^2 :> "d_"[mu, mu],
    MetricTensor -> "d_",
    LeviCivita -> Eps,
    ScalarProduct -> scalar,
    PropagatorDenominator -> prop,
    FeynAmpDenominator -> loop,
    ChiralityProjector[c_] :> ga[(13 - c)/2],
    (DiracSlash | DiracMatrix)[p_] :> MomThread[p, ga],
    NonCommutative -> noncomm,
    FourVector -> fvec };

  ampden = If[ combden === False, Identity,
    DenList = Sort[DenList, Length[#1] >= Length[#2] &];
    AmpDen ];
  amps = LevelSelect[lev]@@@ amps;
  Apply[ (amps[[##, 0]] = DenExtend)&,
    Sort[Thread[-Length@@@ Extract[amps, #] -> #]& @
      Position[amps, _ExtDen]], {2} ];

  intmax = Max[0, Cases[amps, i_intM :> Length[i], Infinity]];

  FormVectors = Join[FormVectors,
    Complement[Flatten[vecs], FormVectors]];

  indsym[type_, n_] := (
    indices = {indices, #};
    ranges = {IndexDim[type, #], ranges};
    indsym[type, n] = #
  )& @ iname[type, n];

  hd = Amp[Apply[{#1, k[++kc], ##3}&, proc /. Index -> indsym, {2}]];

  amps = DeleteCases[amps, {___, _ -> 0, __}];
  If[ Length[amps] === 0, Return[hd[0]] ];

  mmains = Cases[DownValues[inssym], _[_[_[ins_]], s_Symbol] :> s == ins];

  amps = ToCat[2, amps] /. NoExpandRule /. FinalFormRules;
  If[ FreeQ[amps[[2]], Rule], inspol = Begin ];

  haveferm = !FreeQ[amps, FermionChain];

  amps = FLines[amps /. Index -> indsym];
  mmains = mmains /. Index -> indsym;
  ranges = Append[Union[Flatten[ranges]] /. Index -> indsym, i_ -> i == 0];

  indices = Union[Flatten[ {indices, sM,
    Cases[DownValues[indsym], _[_, s_Symbol] -> s],
	(* possibly some colour or gluon indices have been fixed by
	   the user in FeynArts; these would by default be declared
	   as symbols and not be seen by FORM: *)
    Cases[amps, SUNObjs[o__] :> Cases[{o}, _Symbol], Infinity]} ]];

  If[ Head[forder] === Colour,
    indices = Join[#, Complement[indices, #]]&[
      iname[Colour, #]&/@ List@@ forder ];
    forder = Head[forder] ];

  If[ TrueQ[g5eps],
    amps = amps /. {"g5_" -> g5M, "g6_" -> g6M, "g7_" -> g7M} ];

  vars = FormVars[ If[dim === 4, 4, D],
    {amps, SUNN, Hel[], dirM[]}, indices, ranges ];

  hh = OpenForm["fc-" <> tag <> "-"];
  WriteString[hh, "\
#define Dim \"", ToString[dim], "\"\n\
#define InsertionPolicy \"" <> ToForm[inspol] <> "\"\n\
#define IntMax \"" <> ToForm[intmax] <> "\"\n\
#define ExtMax \"" <> ToForm[Max[intmax, extmax]] <> "\"\n\
#define MomElim \"" <> ToString[momelim && Length[MomSubst] === 0] <> "\"\n\
#define DotExpand \"" <> ToBool[dotexp] <> "\"\n\
#define NoCostly \"" <> ToBool[nocost] <> "\"\n\
#define HaveFermions \"" <> ToBool[haveferm] <> "\"\n\
#define FermionOrder \"" <> ToSeq[forder] <> "\"\n\
#define FermionChains \"" <> ToForm[fchain] <> "\"\n\
#define Gamma5Test \"" <> ToBool[g5test] <> "\"\n\
#define Evanescent \"" <> ToBool[evanes] <> "\"\n\
#define HaveSUN \"" <> ToBool[!FreeQ[amps, SUNObjs]] <> "\"\n\
#define SUNN \"" <> ToForm[SUNN] <> "\"\n\
#define SortDen \"" <> ToBool[sortden] <> "\"\n\
#define CombineDen \"" <> ToForm[combden] <> "\"\n\
#define PaVeReduce \"" <> ToForm[pavered] <> "\"\n\
#define SimplifyQ2 \"" <> ToBool[simpq2] <> "\"\n\
#define OPP \"" <> ToForm[Max[2, Abs[opp]]] <> "\"\n\
#define OPPRat \"" <> ToBool[opp < 0] <> "\"\n\
#define OPPQSlash \"" <> ToBool[oppqsl] <> "\"\n\n" <>
    vars[[1]] <> "\n\
table HEL(1:" <> ToString[Legs] <> ");\n" <>
    Array[{"fill HEL(", ToString[#], ") = ", ToForm[Hel[#]], ";\n"}&,
      Legs] <> "\n\n"];

  FormWrite[hh, amps[[1]] /. MomSubst];

  WriteString[hh,
    FormProcs <>
    FormConst[vars, nobrk] <>
    "#procedure Insertions\n"];
  res = Plus@@ FormWrite[ hh, amps[[2]] ];
  WriteString[hh, ".sort\ndrop;\n\n"];
  FormWrite[hh, FormCalc`Result -> res];
  WriteString[hh, "#endprocedure\n\n" <>
    FormCode["CalcFeynAmp.frm"]];
  Close[hh];

  nobrk = Alt[nobrk];	(* for DotSimplify *)
  FormPre[amps];

  hd@@ post/@ FormOutput[mmains][edit, retain][[1]]
]

CalcFeynAmp[amps__, opt___Rule] := CalcFeynAmp[
  DeclareProcess[amps, FilterOptions[DeclareProcess, opt]],
  FilterOptions[CalcFeynAmp, opt] ]


(* FORM interfacing *)

trace[g__] := (traces = Max[traces, ++fline];
  NonCommutativeMultiply[ga[], g] /. {
    a_. ga[li__] ** ga[6] + a_. ga[li__] ** ga[7] -> a ga[li],
    a_. ga[6] + a_. ga[7] -> a
  } /. ga -> om)

chain[g1_, g___] := (fline = Max[fline + 1, 100];
  NonCommutativeMultiply[g1, ga[], g] /. ga -> om)

dobj[g1_, g___][di__] :=
Block[ {fline = sM},
  dirM[NonCommutativeMultiply[g1, ga[], g] /. ga -> om, di]
]


om[1] := "g_"[fline]

om[5] := "g5_"[fline]

om[6] := "g6_"[fline]/2

om[7] := "g7_"[fline]/2

om[li___] := "g_"[fline, li]


Attributes[FLines] = {Listable}

FLines[lhs_ -> rhs_] :=
Block[ {fline = 0},
  lhs -> (rhs /.
    MatrixTrace -> trace /.
    FermionChain -> chain /.
    DiracObject -> dobj /.
    NonCommutativeMultiply[a_] -> a)
]

FLines[other_] := other


Attributes[FormWrite] = {Listable}

FormWrite[hh_, lhs_ -> rhs_] :=
Block[ {fline = 0},
  Write[hh, "L ", lhs == rhs, ";"];
  WriteString[hh, "\n"];
  lhs
]

FormWrite[_, lhs_] := lhs


FormDecl[_, _[]] = {}

FormDecl[type_, _[f_, v___]] :=
Block[ {l, ll, t},
  ll = StringLength[t = type <> ToForm[f]];
  { Fold[
      ( l = StringLength[t = ToForm[#2]] + 2;
        {#1, If[(ll += l) > 70, ll = l; ",\n  ", ", "], t} )&,
      t, {v} ],
    ";\n" }
]


Attributes[FormId] = Attributes[FormSq] = {Listable}

FormId[_[lhs_, rhs_]] :=
  {"id ", ToForm[lhs], " = ", ToForm[rhs], ";\n"}

FormSq[_[lhs_, h_[x___]]] :=
  {"id ", #1, "^2 = ", #2, ";\n",
   "id ", #1, "^-2 = 1/(", #2, ");\n"}&[ ToForm[lhs], ToForm[h[x]] ] /;
  Context[h] === "System`"

FormSq[_[lhs_, rhs_]] :=
  {"id ", #1, "^2 = ", #2, ";\n",
   "id ", #1, "^-2 = ", #2, "^-1;\n"}&[ ToForm[lhs], ToForm[rhs] ]


NCFuncs = Spinor | g5M | g6M | g7M

FormVars[dim_, expr_, inds_:{}, ranges_:{}, cvars_:{},
  vecs_:FormVectors, tens_:FormTensors] :=
Block[ {theexpr, vars, func},
  theexpr = {expr, dim, FormSymbols};

  vars = Complement[Symbols[theexpr],
    cvars, inds, FormVectors, FormTensors];

  func = Complement[
    Cases[Head/@ Level[theexpr, {1, -2}, Heads -> True], _Symbol],
    cvars, FormVectors, FormTensors, ExtWF,
    {ga, MatrixTrace, FermionChain, NonCommutativeMultiply,
     Rule, Equal, Plus, Times, Power, Dot, Rational} ];

  { { FormDecl["s ", vars],
      FormDecl["cf ", DeleteCases[func, NCFuncs]],
      FormDecl["f ", Cases[func, NCFuncs]],
      If[ dim === 4, "", "d D;\n"],
      FormDecl["i ", Replace[inds, ranges, 1]],
      "\n#define Vectors \"", ToSeq[vecs], "\"\nv `Vectors';\n",
      "\n#define Tensors \"", ToSeq[tens], "\"\nt `Tensors';\n" },
    inds, vars, func }
]


FormConst[vars_, nobrk_] := {
  "#procedure ConstBracket\n",
  FormDecl["ab ", DeleteCases[
    Complement[Flatten @ vars[[{3, 4}]], Flatten[{D, Dminus4, nobrk}]],
      (x_ /; MemberQ[{"FeynArts`", "FormCalc`", "LoopTools`", "Form`"},
        Context[x]]) | SUNObjs | Conjugate ]],
  "#endprocedure\n\n" }


tempnum = 1

FormFile[stub_, n_] :=
  ToFileName[Directory[], stub <> ToString[n] <> ".frm"]

OpenForm[stub_:"fc"] :=
Block[ {hh},
  While[ FileType[tempfile = FormFile[stub, tempnum]] =!= None,
    ++tempnum ];
  Print[""];
  Print["preparing FORM code in ", tempfile];
  hh = OpenWrite[toform <> Escape[tempfile],
    FormatType -> InputForm, PageWidth -> 73];
  WriteString[hh, FormSetup];
  hh
]

toform = "!" <> Escape[ToFileName[$FormCalcBin, "ToForm"]] <> " > "


Attributes[FormExpr] = {HoldAll}

FormExpr[x__] := Block[{addM = Plus, mulM = Times, powM = Power, subM}, {x}]


FormPre = AbbrevSet[#, Preprocess -> FormMat]&

FormSub = AbbrevDo[#, 0]&

FormDot = DotSimplify[Simplify, TermCollect]

FormMat = TermCollect

FormNum = TermCollect

FormQC = DenCancel


Attributes[FormExec] = {HoldAll}

FormExec[s__Set] := Block[{s}, ReadForm[$FormCmd, tempfile]]


FormOutput[r___][edit_, retain_] :=
Block[ {res},
  Switch[ edit,
    True, Pause[1]; Run[StringForm[$Editor, tempfile]]; Pause[3],
    Modal, Pause[1]; Run[StringForm[$EditorModal, tempfile]] ];
  WriteString["stdout", "running FORM... "];
  If[ Check[
        res = Set@@@ Level[
          {r, FromFormRules, {Dot -> p$$}},
          {2}, FormExec ];
        True,
        False] &&
    !retain, DeleteFile[tempfile] ];
  Print["ok"];
  res
]


(* things to do when the amplitude comes back from FORM *)

_dummy$$ = 1

d$$ = MetricTensor

i$$ = I

e$$ = Eps

	(* different operator priority in FORM: *)
p$$[a_, b_^n_] := Pair[a, b]^n

p$$[x__] := Pair[x]


paveM[n_, {i__}, args__] := PaVeIntegral[[n]][paveid[n, i], args]

paveid[n_, i__] := paveid[n, i] =
  ToSymbol["LoopTools`", #, #, i]& @ FromCharacterCode[n + 96]

cutM[n_, {num__}, hel_, {args__}, pm__] :=
  CutIntegral[[n]][args][hel, num, pm]


A0[0] = 0

B0i[id:bb0 | dbb0, p_, m1_, m2_] :=
  B0i[id, p, m2, m1] /; !OrderedQ[{m1, m2}]

(Derivative[0, 1, 0, 0][B0i][#1, args__] := B0i[#2, args])&@@@
  {{bb0, dbb0}, {bb1, dbb1}, {bb00, dbb00}, {bb11, dbb11}}


(* abbreviationing business *)

Mat[0] = 0

fermM[x_] := fermM[x] = NewSymbol["F"]


pair[x_] := pair[x] = NewSymbol["Pair"]

eps[x_] := eps[x] = NewSymbol["Eps"]


abb[x_] := x /; Depth[x] < 4

abb[n_?NumberQ x_] := n abb[x]

(*abb[x_?AtomQ] = x*)

abb[x_^n_Integer] := abb[x]^n

abb[x_] := abb[x] = NewSymbol["Abb"]


abbM[x_] := x /; FreeQ[x, Pair | Eps]

abbM[x_] := abb[x /. {p_Pair :> pair[p], e_Eps :> eps[e]}]

dotM[x_] := abbM[FormDot[x]]


sunM[x_] := sunM[x] = NewSymbol["SUN"]


numM[x_, i___] := num2[i]@@ CoefficientList[
  x /. Den[a_, b_]^n_. (a_ - b_) :> Den[a, b]^(n - 1),
  dm4M ]

num2[i___][x_, y_:0, ___] := { qcount[x],
  num[Num[FormNum[x]], i],
  num[Num[FormNum[y]], i] }

num[x_] := num[x] = NewSymbol["Num"]

num[x_, i_List] := Level[{{x}, Select[i, !FreeQ[x, #]&]}, {2}, num]

num[x_, i__] := num[x, i] = NewSymbol["Num"][i]


qcount[p_Plus] := Max[qcount/@ List@@ p]

qcount[t_Times] := Plus@@ qcount/@ List@@ t

qcount[x_] := Count[x, q1, {-1}]


qfM[x_] := x /. {p_Pair :> pair[p] /; FreeQ[p, q1],
                 e_Eps :> eps[e] /; FreeQ[e, q1]}


qcM[x_] := qc[FormQC[x]]

qc[x_] := x /; LeafCount[x] < 10

qc[n_?NumberQ x_, i___] := n qc[x, i]

qc[p_Plus, i___] := -qc[-p, i] /; MatchQ[p[[1]], _?Negative _.]

qc[x_] := qc[x] = NewSymbol["QC"]

qc[x_, i_List] := Level[{{x}, Select[i, !FreeQ[x, #]&]}, {2}, qc]

qc[x_, i__] := qc[x, i] = NewSymbol["QC"][i]


Attributes[dv] = {HoldAll}

dv[x_, _] := {} /; !FreeQ[x, Pattern]

dv[_[_[x_, ___]], s_] := s -> x


Abbr[] := Flatten @ Apply[dv, DownValues/@
  {fermM, pair, eps, abb, sunM, num, qc, xi}, {2}]

Abbr[patt__] :=
Block[ {all = Abbr[], need = Flatten[{patt}], omit},
  omit = Cases[need, !x_ -> x];
  need = DeleteCases[need, !_];
  FixedPoint[ abbsel[First/@ #]&, abbsel[need, omit] ]
]

abbsel[{need__}, {omit__}] := Select[all,
  !FreeQ[#, Alternatives[need]] && FreeQ[#, Alternatives[omit]]&]

abbsel[{need__}, ___] := Select[all, !FreeQ[#, Alternatives[need]]&]


IndexHeader[h_, {}] := h

IndexHeader[h_, {i__}] := h[i]

IndexHeader[h_[i___], expr__] :=
  IndexHeader[h, Flatten[{i,
    Select[ Union @ Symbols[Level[{expr}, {-2}]],
      Head[Dim[#]] =!= Dim & ]}]]


MomEncode[other_] := other /; FreeQ[other, k]

MomEncode[p_Plus] := MomEncode/@ p

MomEncode[f_. k[i_]] := MomEncoding[f, i]


(* AbbrevInt introduces abbreviations for the loop integrals.
   They fall into two categories:
   1. A0, A00,
   2. B0i, C0i, D0i, E0i, F0i.
   For the latter the LoopTools functions [BCDEF]put is used to
   compute all tensor coefficients at once. *)

Function[ {X0i, Xput, Nxx},
  AbbrevInt[X0i[i_, args__]] :=
  Block[ {abb},
    abbint[X0i[id_, args]] = abb[id_] =
      IndexHeader[ToSymbol["pave", ++pavec][id], args];
    lint = {lint, abb[NDim[Nxx]] -> (Xput[#, args]&)};
    abb[i]
  ] ]@@@
{ {B0i, Bput, Nbb},
  {C0i, Cput, Ncc},
  {D0i, Dput, Ndd},
  {E0i, Eput, Nee},
  {F0i, Fput, Nff} }

Function[ {Xcut, Xmas, Mxx},
  AbbrevInt[Xcut[margs__][nargs__]] :=
    abbint[Xcut[abbint[Xmas[margs]], nargs]];
  AbbrevInt[Xcut[args__]] :=
  Block[ {abb},
    abbint[Xcut[args]] = abb =
      IndexHeader[ToSymbol["cut", ++cutc][], args];
    lint = {lint, abb -> MomEncode/@ Xcut[args]};
    abb
  ];
  AbbrevInt[Xmas[args__]] :=
  Block[ {abb},
    abb[id_] = IndexHeader[ToSymbol["mas", ++masc][id], args];
    lint = {lint, abb[NDim[Mxx]] -> (Xmas[#, args]&)};
    abbint[Xmas[args]] = abb[1]
  ] ]@@@
{ {Bcut, Bmas, Mbb},
  {Ccut, Cmas, Mcc},
  {Dcut, Dmas, Mdd},
  {Ecut, Emas, Mee},
  {Fcut, Fmas, Mff} }

AbbrevInt[func_] :=
Block[ {abb = IndexHeader[ToSymbol["pave", ++pavec][], func]},
  lint = {lint, abb -> func};
  abbint[func] = abb
]


ExtractInt[expr_] :=
Block[ {abbint, new, lint = {}, lc = 0},
  abbint[f_] := AbbrevInt[f];
  new = expr /. int:LoopIntegral[__] :> abbint[int];
  {Flatten[lint], new}
]


GenNames[amp_] := amp /. {
  g:G[_][_][__][__] :> coup[g],
  m_Mass :> mass[m],
  x_GaugeXi :> xi[x],
  v:VertexFunction[_][__] :> vf[v] }

coup[g_] := coup[g] = NewSymbol["Coupling"]

mass[m_] := mass[m] = NewSymbol["Mass"]

xi[x_] := xi[x] = NewSymbol["GaugeXi"]

vf[v_] := vf[v] = NewSymbol["VertexFunction"]


GenericList[] := Flatten @
  Apply[dv, DownValues/@ {coup, mass, xi, vf}, {2}]


Options[ClearProcess] = {
  ZapFunction -> Remove
}

ClearProcess[opt___Rule] :=
Block[ {zapf},
  {zapf} = ParseOpt[ClearProcess, opt];
  Apply[zap, DownValues/@ {fermM, pair, eps, abb, sunM,
    num, qc, coup, mass, xi, vertex}, {2}];
  CurrentProc = CurrentOptions = DenList = {};
  Clear[SymbolNumber];
  _SymbolNumber = 0;
]


Attributes[zap] = {HoldAll}

zap[p_, s_Symbol] := (Unset[p]; zapf[s]) /; FreeQ[p, Pattern]

zap[p_, s_Symbol[__]] := (Unset[p]; zapf[s]) /; FreeQ[p, Pattern]


subdef[minleaf_, deny_, fuse_, pre_, ind_] := (
  DownValues[subadd] = {
    HoldPattern[subadd[x_^2 - y_^2]] :> subadd[x - y] subadd[x + y],
    HoldPattern[subadd[x_]] :> (abbsub@@ Select[ind, !FreeQ[x, #]&])[x] /;
      LeafCount[x] > minleaf && FreeQ[x, deny],
    HoldPattern[subadd[x_]] :> x
  } /. {
    HoldPattern[abbsub@@ Select[{}, _]] -> abbsub[],
    HoldPattern[x_ && FreeQ[_, _[]]] -> x };

  DownValues[subfun] = {
    HoldPattern[subfun[x_]] :>
      (abbsub@@ Select[ind, !FreeQ[x, #]&])[pre[x]] /;
      subF[x] && LeafCount[x] > minleaf && FreeQ[x, deny],
    HoldPattern[subfun[s_SumOver x_]] :> s subfun[x],
    HoldPattern[subfun[h_[x_]]] :> h[subfun[x]],
    HoldPattern[subfun[h_[x__]]] :> h[
      subfun[Select[h[x], !subF[#]&], h],
      subfuse[Select[h[x], subF], h]
    ] /; Length[Intersection[Attributes[h], {Flat, Orderless}]] === 2,
    HoldPattern[subfun[x_]] :> subfun/@ x,
    HoldPattern[subfun[h_[x___], h_]] :> subfun/@ h[x],
    HoldPattern[subfun[x_, _]] :> subfun[x]
  } /. {
    HoldPattern[abbsub@@ Select[{}, _]] -> abbsub[],
    HoldPattern[x__ && FreeQ[_, _[]]] :> And[x],
    HoldPattern[subfuse[x_, _]] :> subadd[x] /; fuse };

  DownValues[AbbrevDo] = {
    HoldPattern[AbbrevDo[expr_, lev_Integer]] :>
      Replace[expr, Plus -> sublev, {lev, Infinity}, Heads -> True],
    HoldPattern[AbbrevDo[expr_, f_]] :>
      Block[{subF = f}, subfun[expr]]
  } /. sublev :> (pre[Plus[##]] /. Plus -> sublev &) /; pre =!= Identity
)


AbbrevSet[expr_, opt___Rule] :=
Block[ {minleaf, deny, fuse, pre},
  {minleaf, deny, fuse, pre} = ParseOpt[Abbreviate, opt];
  subdef[minleaf, Alt[deny], fuse, pre,
    Cases[expr, SumOver[i_, ___] -> i, Infinity] //Union];
]


AbbrevDo::unset = "Run AbbrevSet before AbbrevDo."

_AbbrevDo := (Message[AbbrevDo::unset]; Abort[])


Options[Abbreviate] = {
  MinLeafCount -> 10,
  Deny -> {k, q1},
  Fuse -> True,
  Preprocess -> Identity }

Abbreviate[expr_, x_:2, opt___Rule] :=
Block[ {subadd, subfun, AbbrevDo},
  AbbrevSet[expr, opt];
  AbbrevDo[expr, x]
]


sublev[x:(_Integer | _Rational) _., r__] :=
 (#1 subadd[Block[{plus = Plus}, #2]]&)@@
    FactorTermsList[Plus@@ ({x, r} /. Plus -> plus)]

sublev[p___] := subadd[Plus[p]]


subfuse[h_[x___], h_] := subadd/@ h[x]

subfuse[x_, _] := subadd[x]


$AbbPrefix = "Sub"

abbsub[][x_] := abbsub[][x] = NewSymbol[$AbbPrefix]

abbsub[i__][x_] := abbsub[i][x] = NewSymbol[$AbbPrefix][i]


Subexpr[args__] :=
Block[ {abbsub, res},
  abbsub[][x_] := abbsub[][x] = NewSymbol[$AbbPrefix];
  abbsub[i__][x_] := abbsub[i][x] = NewSymbol[$AbbPrefix][i];
  res = Abbreviate[args];
  {Subexpr[], res}
]

Subexpr[] := Cases[ SubValues[abbsub],
  _[_[_[x_]], s_] :> (s -> x) /; FreeQ[x, Pattern] ]


Options[ClearSubexpr] = {
  ZapFunction -> Remove
}

ClearSubexpr[opt___Rule] :=
Block[ {zapf},
  {zapf} = ParseOpt[ClearProcess, opt];
  zap@@@ SubValues[abbsub];
]


likeQ[patt_][s_] := StringMatchQ[ToString[s], patt]


General::unknown = "Don't know how to register ``."

General::conflict = "Definition of `` conflicts with ``."

Attributes[RegisterAbbr] = {Listable}

RegisterAbbr[s_ -> x_Num] := regset[num, x, s]

RegisterAbbr[s_ -> other_] := regabb[s -> 4711 other]

RegisterAbbr[other_] := (Message[RegisterAbbr::unknown, other]; other)

regabb[s_ -> x_. p___Pair e___Eps d___DiracChain w___WeylChain t___eT tc___eTc] :=
  regabb[s, x/4711, Times[p], Times[e], d w, t tc]

regabb[s_, x_, 1, 1, 1, t_] := regabb[s, x t]

  regabb[s_?(likeQ["SUN*"]), x_] := regset[sunM, x, s]

  regabb[s_?(likeQ["Abb*"]), x_] := regset[abb, x, s]

  regabb[s_?(likeQ["QC*"]), x_] := regset[qc, x, s]

regabb[s_, x_, p_, 1, 1, 1] := regabb[s, x, p, pair]

regabb[s_, x_, 1, e_, 1, 1] := regabb[s, x, e, eps]

regabb[s_, x_, p_, e_, f_, t_] := regabb[s, x, p e f t, fermM]

  regabb[s_, 1, t_, h_] := regset[h, t, 1, s]

  regabb[s_, x_, t_, h_] := regset[h, t, x^(-2 Count[4711 t, _k, {2}]), s]

regabb[s_, x__] :=
  (Message[RegisterAbbr::unknown, #]; #)&[ s -> Times[x] ]


regset[h_, args__, s_] :=
  newset[h[args], Cases[DownValues[h], _[_[_[args]], x_] -> x], s]

Attributes[newset] = {HoldFirst}

newset[x_, {}, s_] := x = s

newset[x_, {s_}, s_] := s

newset[x_, {s1_}, s2_] :=
  (Message[RegisterAbbr::conflict, s1, s2]; s1)


Attributes[RegisterSubexpr] = {Listable}

RegisterSubexpr[s_Symbol[i__] -> x_] :=
  regset[abbsub[i], x, s]

RegisterSubexpr[s_Symbol -> x_] :=
  regset[abbsub[], x, s]

RegisterSubexpr[other_] :=
  (Message[RegisterSubexpr::unknown, other]; other)


$OptPrefix = "Opt"

newvar[] := NewSymbol[$OptPrefix]

newvar[i__] := NewSymbol[$OptPrefix][i]


Attributes[set] = {Flat, Orderless}


ptdef[0][x_, t_] := lotdef[x, lot@@ t]

ptdef[1][x_, p_] := (lotdef[-x, FromPlus[lot, -p]]; lotdef[x, lot@@ p])

lotdef[x_, l_lot] := (l = x; {})

lotdef[x_, _Integer x_] := x -> 0

lotdef[x_, other_] := x -> other


lotget[_, _Times] := {{}, {}}

lotget[_[_[terms__]], x_] := {x, set[terms]}


Attributes[ptcom] = {Listable}

ptcom[a_set, b_set, 0] := Intersection[a, b]

ptcom[a_set, b_set, 1] :=
Block[ {cp, cm},
  cp = Intersection[a, b];
  If[ Length[a] > 2 Length[cp],
    cm = Intersection[a, Thread[-b, set]];
    If[ Length[cm] > Length[cp], cp = cm ] ];
  cp
]

_ptcom = set[]


ptrul[a_, b_, 0] := a -> b

ptrul[a_, b_, 1] := {a -> b, Thread[-a, set] -> -b}


Overlap[] = set[]

Overlap[x__] := Intersection[x]


CSE[_, {}] = {}

CSE[pt_, abbr_] :=
Block[ {lot, alias, var, def, com, tmp, new = {}},
  Attributes[lot] = {Flat, Orderless};
  alias = Flatten[ ptdef[pt]@@@
    Sort[abbr, Length[ #1[[2]] ] <= Length[ #2[[2]] ] &] ];
  {var, def} = ToCat[2, Apply[lotget, DownValues[lot], 1]];
  def = def /. alias;

  Do[
    While[ Length[com = ptcom[def[[i]], def[[i + 1]], pt]] > 3,
      tmp = Ceiling[Length[com]/2];
      tmp = Overlap@@ Select[
        ptcom[com, Drop[def, {i, i + 1}], pt],
        Length[#] > tmp & ];
      If[ Length[tmp] > 3, com = tmp ];
      tmp = Position[def, com, 1, 1, Heads -> False];
      If[ Length[tmp] > 0,
        tmp = tmp[[1, 1]];
        def[[tmp, 0]] = List;
        def = def /. ptrul[com, var[[tmp]], pt];
        def[[tmp, 0]] = set,
      (* else *)
        tmp = newvar@@ Select[
          Union[Level[{var[[i]], var[[i + 1]]}, {2}]],
          !FreeQ[com, #]& ];
        new = {new, tmp -> com};
        def = def /. ptrul[com, tmp, pt] ]
    ],
  {i, Length[def] - 1}];

  Flatten[{new, Thread[var -> def], alias}]
]


AbbrCat[rul:_[_, _Plus]] := {{}, {}, rul}

AbbrCat[rul:_[_, t_Times]] := {{}, rul, {}} /; FreeQ[t, DiracChain]

AbbrCat[rul_] := {rul, {}, {}}


OptimizeAbbr[{}] = {}

OptimizeAbbr[rul:{__Rule}] := Flatten[
  {#1, CSE[0, #2] /. set -> Times, CSE[1, #3] /. set -> Plus}&@@
    ToCat[3, AbbrCat/@ rul] ]


Attributes[simple] = {Listable}

simple[r:(_ -> _?NumberQ)] := subst[r] = Sequence[]

simple[r:(_ -> n_. _Symbol)] := (subst[r] = Sequence[]) /; NumberQ[n]

simple[CodeExpr[vars__, expr_]] := CodeExpr[vars, simple[expr]]

simple[other_] := other

SubstSimpleAbbr[arg_] :=
Block[ {subst, new = arg, rul},
  While[ ( new = simple[new];
           Length[rul = (#[[1, 1, 1]]&)/@ DownValues[subst]] > 0 ),
    new = new /. CodeExpr[var_, tmpvar_, expr_] :>
      (CodeExpr[DeleteCases[var, #], DeleteCases[tmpvar, #], expr]&[
        First/@ rul ]) //. rul;
    Clear[subst]
  ];
  new
]


(* UV and IR finiteness checks *)

ToNewBRules = {
  B0[args__] -> B0i[bb0, args],
  B1[args__] -> B0i[bb1, args],
  B00[args__] -> B0i[bb00, args],
  B11[args__] -> B0i[bb11, args],
  B001[args__] -> B0i[bb001, args],
  B111[args__] -> B0i[bb111, args],
  DB0[args__] -> B0i[dbb0, args],
  DB1[args__] -> B0i[dbb1, args],
  DB00[args__] -> B0i[dbb00, args],
  DB11[args__] -> B0i[dbb11, args],
  C0[args__] -> C0i[cc0, args],
  D0[args__] -> D0i[dd0, args],
  E0[args__] -> E0i[ee0, args],
  F0[args__] -> F0i[ff0, args] }

ToOldBRules = {
  B0i[bb0, args__] -> B0[args],
  B0i[bb1, args__] -> B1[args],
  B0i[bb00, args__] -> B00[args],
  B0i[bb11, args__] -> B11[args],
  B0i[bb001, args__] -> B001[args],
  B0i[bb111, args__] -> B111[args],
  B0i[dbb0, args__] -> DB0[args],
  B0i[dbb1, args__] -> DB1[args],
  B0i[dbb00, args__] -> DB00[args],
  B0i[dbb11, args__] -> DB11[args] }


Attributes[FindDiv] = {Listable}

FindDiv[rc_ -> x_] := rc -> FindDiv[x]

FindDiv[i_IndexIf] := MapIf[FindDiv, i]

FindDiv[amp:Amp[_][___]] := FindDiv/@ amp

FindDiv[other_] := foo[other]


SubDiv[p_Plus] := SubDiv/@ p

SubDiv[x_ r_] := x SubDiv[r] /; FreeQ[x, div]

SubDiv[x_] := x - (x /. div -> 0)


SerDiv[x_] := Series[x /. Divergence -> -2/Dminus4, {Dminus4, 0, 0}]


MapDiv[f_, expr_, fin_] :=
Block[ {div = RCPattern[Divergence], foo = f, Finite = fin},
  FindDiv[expr /. ToNewBRules /.
    {int:PaVeIntegral[__] :> UVDiv[int] + fin int, D -> Dminus4 + 4} /.
    Dminus4^(n_?Negative) -> (-2/Divergence)^n]
]


UVDivergentPart[expr_] := MapDiv[SubDiv, expr, 0]


UVSeries[expr_, pow_] := Coefficient[UVSeries[expr], Dminus4, pow]

UVSeries[expr_] := MapDiv[SerDiv, expr, 1]


UVDiv[A0[m_]] = m Divergence

UVDiv[A00[m_]] = m^2/4 Divergence

UVDiv[B0i[bb0, __]] = Divergence

UVDiv[B0i[bb1, __]] = -1/2 Divergence

UVDiv[B0i[bb00, p_, m1_, m2_]] := -(p - 3 m1 - 3 m2)/12 Divergence

UVDiv[B0i[bb11, __]] = 1/3 Divergence

UVDiv[B0i[bb001, p_, m1_, m2_]] := (p - 2 m1 - 4 m2)/24 Divergence

UVDiv[B0i[bb111, __]] = -1/4 Divergence

UVDiv[B0i[dbb00, __]] = -1/12 Divergence

UVDiv[C0i[cc00, __]] = 1/4 Divergence

UVDiv[C0i[cc001 | cc002, __]] = -1/12 Divergence

UVDiv[C0i[cc0000, p1_, p2_, p1p2_, m1_, m2_, m3_]] :=
  -(p1 + p2 + p1p2 - 4 (m1 + m2 + m3))/96 Divergence

UVDiv[C0i[cc0011 | cc0022, __]] = 1/24 Divergence

UVDiv[C0i[cc0012, __]] = 1/48 Divergence

UVDiv[D0i[dd0000, __]] = 1/24 Divergence

UVDiv[D0i[dd00001 | dd00002 | dd00003, __]] = -1/96 Divergence

_UVDiv = 0


(* FeynCalc compatibility functions *)

PaVe[i__Integer, {p__}, {m__}] :=
  ToExpression[#1 <> "0i"][ ToSymbol[#2, #2, Sort[{i}]], p, m ]&[
    FromCharacterCode[Length[{m}] + 64],
    FromCharacterCode[Length[{m}] + 96] ]


FeynCalcGet[mask___] :=
Block[ {Global`GraphName, Global`Momentum = Identity},
  _Global`GraphName = 0;
  Plus@@ Get/@ FileNames[mask] /. ep_Eps :> I ep
]


FeynCalcPut[expr_, file_] :=
Block[ {C0i, D0i, E0i, F0i, PaVe},
  C0i[cc0, args___] := C0[args];
  D0i[dd0, args___] := D0[args];
  E0i[ee0, args___] := E0[args];
  F0i[ff0, args___] := F0[args];
  C0i[i_, p__, m1_, m2_, m3_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3} ];
  D0i[i_, p__, m1_, m2_, m3_, m4_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3, m4} ];
  E0i[i_, p__, m1_, m2_, m3_, m4_, m5_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3, m4, m5} ];
  F0i[i_, p__, m1_, m2_, m3_, m4_, m5_, m6_] := PaVe[
    Sequence@@ (ToExpression/@ Drop[Characters[ToString[i]], 2]),
    {p}, {m1, m2, m3, m4, m5, m6} ];
  Put[expr /. ToOldBRules, file]
]


(* helicity matrix elements *)

ToTrace[fi_ -> plain_, fj_ -> conj_] := Mat[fi, fj] ->
  plain (conj /. DiracChain -> ConjChain /. ep_Eps -> -ep /. ConjWF)

ConjChain[s1_Spinor, om_Integer, g___, s2_Spinor] :=
  (#1 Reverse[DiracChain[s1, g, #2, s2]])&@@
     ConjOm[EvenQ[Length[{g}]], om]

ConjOm[False, om_] := {1, om};
ConjOm[_, 6] = {1, 7};
ConjOm[_, -6] = {1, -7};
ConjOm[_, 7] = {1, 6};
ConjOm[_, -7] = {1, -6};
ConjOm[_, 5] = {-1, 5};
ConjOm[_, -5] = {-1, -5};
ConjOm[_, 1] = {1, 1}
ConjOm[_, -1] = {1, -1}


HelTab[Hel[i_]] := helM[i, e]

HelTab[h_] := h + e


SelectAbbr[abbr_, All] = abbr

SelectAbbr[abbr_, expr_] := Select[abbr, !FreeQ[ expr, #[[1]] ]&]

SelectAbbr[abbr_, plain_, _[0]] :=
  {#, #}& @ SelectAbbr[abbr, plain]

SelectAbbr[abbr_, plain_, conj_] :=
  {Union[SelectAbbr[abbr, plain], #], #}& @ SelectAbbr[abbr, conj]

SelectAbbr[abbr_, patt_, plain_, conj_] :=
  SelectAbbr[Select[abbr, !FreeQ[#, patt]&], plain, conj]


Rep[Spinor[k[i_], _, s_, 2, _]] :=
  {Overscript[FromCharacterCode[(235 - s)/2], "."], i}

Rep[Spinor[k[i_], _, s_, ___]] :=
  {FromCharacterCode[(235 - s)/2], i}

Format[DiracChain[s1_Spinor, om_, g___, s2_Spinor]] :=
  SequenceForm@@ Flatten[{"<", Rep[s1], "|", om,
    {",", #}&/@ {g}, "|", Rep[s2], ">"}]

Format[WeylChain[s1_Spinor, om_, g___, s2_Spinor]] :=
  SequenceForm@@ Flatten[{"(", Rep[s1], "|", om,
    {",", #}&/@ {g}, "|", Rep[s2], ")"}]


Options[HelicityME] = {
  Dimension -> 4,
  MomElim -> MomElim,
  DotExpand -> DotExpand,
  EditCode -> False,
  RetainFile -> False }

HelicityME::noprocess = "No process defined so far.  \
HelicityME works only after DeclareProcess or CalcFeynAmp."

HelicityME::weyl = "Warning: HelicityME does not work on WeylChains.  \
CalcFeynAmp uses DiracChains with the option FermionChains -> Chiral or VA."

General::nomat = "Warning: No matrix elements to compute."

HelicityME[plain_, opt___Rule] := HelicityME[plain, plain, opt]

HelicityME[plain_, conj_, opt___Rule] :=
Block[ {dim, momelim, dotexp, edit, retain,
abbr, part, ind, ic = 0, vars, hels, helM, hh, e, mat},

  If[ CurrentProc === {},
    Message[HelicityME::noprocess];
    Abort[] ];

  {dim, momelim, dotexp, edit, retain} =
    ParseOpt[HelicityME, opt] /. Options[CalcFeynAmp];

  abbr = Abbr[];
  If[ !FreeQ[abbr, WeylChain], Message[HelicityME::weyl] ];

  abbr = SelectAbbr[abbr, DiracChain[_Spinor, ___], plain, conj];
  If[ Times@@ Length/@ abbr === 0,
    Message[HelicityME::nomat];
    Return[{}] ];

  part = Cases[abbr, Spinor[k[i_], __] -> i, Infinity] //Union;

  ind = Map[# -> "N" <> ToString[++ic] <> "_?" &,
    Union[Cases[#, _Lor, Infinity]]&/@ abbr, {2}];

  mat = Flatten[Outer[ ToTrace,
      abbr[[1]] /. ind[[1]],
      abbr[[2]] /. ind[[2]] ]] /.
    Reverse/@ FromFormRules /. Eps -> "e_" /.
    FinalFormRules;

  Print["> ", Length[mat], " helicity matrix elements"];

  hels = HelTab/@ Hel/@ part;
  dim = If[dim === 0, D, 4];
  vars = FormVars[dim, {Last/@ mat, hels}];

  hh = OpenForm["fc-hel-"];
  WriteString[hh, "\
#define MomElim \"" <> ToString[momelim && Length[MomSubst] === 0] <> "\"\n\
#define DotExpand \"" <> ToBool[dotexp] <> "\"\n\n" <>
    vars[[1]] <> "\n\
table HEL(" <> ToString[Min[part]] <> ":" <>
               ToString[Max[part]] <> ", e?);\n" <>
    MapThread[{"fill HEL(", ToForm[#1], ") = ", ToForm[#2], ";\n"}&,
      {part, hels}] <> "\n" <>
    FormProcs <>
    FormCode["HelicityME.frm"]];

  FormWrite[hh, "Mat" <> ToString/@ List@@ #1 -> #2]&@@@ mat;

  WriteString[hh, "#call Emit\n"];
  Close[hh];

  (e[#] = s[#])&/@ part;
  helM[i_, s_] := Hel[i] + s;

  FormPre[mat];

  Thread[First/@ mat -> Plus@@@ FormOutput[][edit, retain]]
]


(* colour matrix elements *)

sunT[i_Symbol, i_] := SUNN

sunT[_, i_Symbol, i_] = 0

sunT[a_Integer, a_, i_Symbol, i_] = 1/2

sunT[_Integer, _Integer, i_Symbol, i_] = 0

sunT[t1___, a_Symbol, t2___, a_, t3___, i_, j_] :=
  (sunT[t1, t3, i, j] sunTr[t2] -
    sunT[t1, t2, t3, i, j]/SUNN)/2

sunT[t1___, a_Symbol, t2___, i_, j_] sunT[t3___, a_, t4___, k_, l_] ^:=
  (sunT[t1, t4, i, l] sunT[t3, t2, k, j] -
    sunT[t1, t2, i, j] sunT[t3, t4, k, l]/SUNN)/2

sunT[a___, i_, j_Symbol] sunT[b___, j_, k_] ^:=
  sunT[a, b, i, k]

sunT[a___, i_, j_Symbol] sunT[b___, k_, j_] ^:=
  Level[{{a}, Reverse[{b}], {i, k}}, {2}, sunT]

sunT[a___, j_Symbol, i_] sunT[b___, j_, k_] ^:=
  Level[{Reverse[{a}], {b}, {i, k}}, {2}, sunT]

sunT/: sunT[a___, i_, j_Symbol]^2 :=
  Level[{{a}, Reverse[{a}], {i, i}}, {2}, sunT]

sunT/: sunT[a___, i_Symbol, j_]^2 :=
  Level[{Reverse[{a}], {a, j, j}}, {2}, sunT]


sunTr[] := SUNN

sunTr[_] = 0

sunTr[a__] := sunT[a, #, #]& @ NewSymbol["col"]


(* we assume that structures of the form delta[a, a] indicate
   summations over external colour/gluon indices *)

sunText[i_Symbol, i_] := Sqrt[SUNN]

sunText[a_Symbol, a_, 0, 0] := Sqrt[SUNN^2 - 1]/2

sunText[a___, 0, 0] := sunTrace[a]

sunText[other__] := sunT[other]


sunF[a_, b_, c_] := 2 I (sunTrace[c, b, a] - sunTrace[a, b, c])

sunF[a__, b_, c_] := (sunF[a, #] sunF[#, b, c])& @ NewSymbol["glu"]


sunEps/: sunEps[a___, _Symbol, b___]^2 :=
  (sunT[#1, #1] sunT[#2, #2] - sunT[#1, #2]^2)&[a, b]

sunEps[a___, i_Symbol, b___] sunEps[c___, i_, d___] ^:=
  (-1)^Length[{a, c}] *
  ((sunT[#1, #3] sunT[#2, #4] -
    sunT[#1, #4] sunT[#2, #3])&[a, b, c, d])


Conjugate[t_SUNT] ^:= RotateLeft[Reverse[t], 2]

Conjugate[f_SUNF] ^= f


ColourSimplify[plain_, conj_:1] :=
Block[ {res, ind},
  {res, ind} = Reap[ plain Conjugate[conj] /.
    {IndexDelta -> idelta, SumOver -> sumover} /.
    SUNSum[i_, _] :> (Sow[i]; 1) ];
  Simplify[ Expand[ res /.
    Cases[ind, i_Index :> i -> iname@@ i, {2}] /.
    {SUNT -> sunText, SUNF -> sunF, SUNEps -> sunEps} /.
    sunTrace[a__]^n_. :> Times@@ Table[sunTr[a], {n}]
  ] /. {sunT[i_, j_] :> Sort[SUNT[i, j]],
        sunT[g__, i_, i_] :> SUNT[g, 0, 0],
        sunT -> SUNT, sunEps -> SUNEps} ]
]


ColourGrouping[tops_] := DiagramGrouping[ tops,
  Replace[
    ColourSimplify[Times@@
      FeynAmpCases[_[Index[Colour | Gluon, _], ___]][##]],
    _?NumberQ r_ -> r]& ]


ColourFactor[fi_ -> plain_, fj_ -> conj_] :=
  Mat[fi, fj] -> ColourSimplify[plain, conj]


ColourME[plain_] := ColourME[plain, plain]

ColourME[plain_, conj_] :=
Block[ {abbr},
  abbr = SelectAbbr[Abbr[], SUNObjs, plain, conj];
  If[ Times@@ Length/@ abbr === 0,
    Message[ColourME::nomat];
    Return[{}] ];

  Outer[ ColourFactor, abbr[[1]], abbr[[2]] ] //Flatten
]


(* squaring the matrix element *)

UniqIndices[conj_, plain_] :=
Block[ {ind},
  ind = Intersection@@
    (Union[Cases[#, SumOver[x_, ___] -> x, Infinity]]&)/@ {conj, plain};
  conj /. Thread[ind -> (ToSymbol[#, "c"]&)/@ ind]
]


SquaredME[amp_] := SquaredME[amp, amp]

SquaredME[plain:Amp[_][___], conj:Amp[_][___]] := (
  ChkProc[{plain, conj}, SquaredME];
  squaredME[Plus@@ plain, Plus@@ conj]
)

squaredME[0, _] = squaredME[_, 0] = 0

squaredME[plain_, conj_] :=
  plain Conjugate[UniqIndices[conj, plain]] /;
  FreeQ[{plain, conj}, Mat]

squaredME[plain_, conj_] :=
  Apply[Plus,
    Outer[ ToSquared,
      MatList[plain], MatList[UniqIndices[conj, plain]] ],
    {0, 1}]


MatList[expr_] := ToList[Collect[expr, _Mat, Hold]]

ToList[p_Plus] := List@@ p

ToList[other_] := {other}


ToSquared[Mat[m1_] x1_., Mat[m2_] x2_.] :=
  ReleaseHold[x1 Conjugate[x2]] ToMat[m1, m2]

ToMat[m1_Symbol, m2_Symbol] := Mat[m1, m2]

ToMat[m1_, m2_] := Inner[Mat, m1, m2, Times]


Unprotect[Conjugate]

Format[ Conjugate[x_] ] := Superscript[x, "*"]

(*
Format[ Conjugate[t_Times] ] :=
  Superscript[SequenceForm["(", t, ")"], "*"]
*)

Conjugate[D] = D

Conjugate[Dminus4] = Dminus4

Re[Dminus4] ^= Dminus4

Re[Finite] ^= Finite

Re[Divergence] ^= Divergence

Conjugate[Finite] ^= Finite

Conjugate[Divergence] ^= Divergence

Conjugate[o_SumOver] := o

Conjugate[d_IndexDelta] := d

Conjugate[p_Plus] := Conjugate/@ p

Conjugate[t_Times] := Conjugate/@ t

Conjugate[d_Den] := Conjugate/@ d

Conjugate[p_Pair] := Conjugate/@ p

Conjugate[ep_Eps] := -Conjugate/@ ep

Conjugate[e[n_]] := ec[n]

Conjugate[ec[n_]] := e[n]

Conjugate[z[n_]] := zc[n]

Conjugate[zc[n_]] := z[n]

Conjugate[eT[n_][i__]] := eTc[n][i]

Conjugate[eTc[n_][i__]] := eT[n][i]

Conjugate[k[n_]] := k[n]

Conjugate[s[n_]] := s[n]

Conjugate[Lor[n_]] := Lor[n]

Conjugate[x_?RealQ] := x

Conjugate[(x_?RealQ)^n_] := x^n

Protect[Conjugate]


Unprotect[Re]

Re[D] = D

Re[Dminus4] ^= Dminus4

Re[Divergence] ^= Divergence

Re[o_SumOver] := o

Re[d_IndexDelta] := d

Re[p_Plus] := Re/@ p

Re[d_Den] := Re/@ d

Re[x_?RealQ y_.] := x Re[y]

Re[(x_?RealQ)^n_ y_.] := x^n Re[y]

Protect[Re]


(* performing the polarization sum analytically *)

Options[PolarizationSum] = {
  Dimension -> 4,
  GaugeTerms -> True,
  MomElim -> MomElim,
  DotExpand -> DotExpand,
  NoBracket -> NoBracket,
  EditCode -> False,
  RetainFile -> False }

PolarizationSum::noprocess = "No process defined so far.  \
PolarizationSum works only after DeclareProcess or CalcFeynAmp."

PolarizationSum::incomp = "PolarizationSum used on an amplitude \
other than the last one set up by DeclareProcess or CalcFeynAmp."

PolarizationSum[amp:Amp[_][___].., opt___?OptionQ] :=
Block[ {Hel},
  ChkProc[{amp, CurrentProc}, PolarizationSum, Abort[]];
  _Hel = 0;
  PolarizationSum[
    SquaredME[amp] /.
      HelicityME[amp, FilterOptions[HelicityME, opt]] /.
      ColourME[amp, FilterOptions[ColourME, opt]],
    opt ]
]

PolarizationSum[expr_, opt___?OptionQ] :=
Block[ {abbr, dim, gauge, momelim, dotexp, nobrk, edit, retain,
fullexpr, lor, indices, legs, masses, vars, hh},

  If[ CurrentProc === {},
    Message[PolarizationSum::noprocess];
    Abort[] ];

  {dim, gauge, momelim, dotexp, nobrk, edit, retain} =
    ParseOpt[PolarizationSum, opt] /. Options[CalcFeynAmp];

  fullexpr = expr //. Dispatch[Subexpr[]] //. Dispatch[Abbr[]] /.
    FinalFormRules;
  lor = Cases[fullexpr, _Lor, Infinity] //Union;
  indices = FormIndices[[ Level[lor, {2}] ]];
  fullexpr = fullexpr /. Thread[lor -> indices];

  legs = Cases[fullexpr, Alt[ExtWF][i_] -> i,
    Infinity, Heads -> True] //Union;
  masses = Masses[CurrentProc][[legs]];

  fullexpr = fullexpr /. s -> e /. Reverse/@ FromFormRules /.
    {Eps -> "e_", MetricTensor -> "d_", Pair -> Dot} /.
    NoExpandRule /.
    FinalFormRules;

  dim = If[dim === 0, D, 4];
  vars = FormVars[dim, {fullexpr, masses}, indices];

  hh = OpenForm["fc-pol-"];
  WriteString[hh, "\
#define Dim \"", ToString[dim], "\"\n\
#define GaugeTerms \"" <> ToBool[gauge] <> "\"\n\n\
#define MomElim \"" <> ToString[momelim && Length[MomSubst] === 0] <> "\"\n\
#define DotExpand \"" <> ToBool[dotexp] <> "\"\n\n" <>
    vars[[1]] <>
    FormProcs <>
    FormConst[vars, nobrk] <>
    FormCode["PolarizationSum.frm"]];

  Write[hh, "L SquaredME = ", fullexpr, ";"];

  WriteString[hh,
    "\n#call Prepare\n" <>
    MapThread[{"\n#call PolSum(", ToString[#1], ", ", ToForm[#2], ", ",
        ToString[If[FreeQ[fullexpr, (z | zc)[#1]], dim, Dminus4]], ")"}&,
      {legs, masses}] <>
    "\n\n#call Emit\n"];
  Close[hh];

  nobrk = Alt[nobrk];	(* for DotSimplify *)
  FormPre[fullexpr];

  Plus@@ FormOutput[][edit, retain][[1]]
]


FormIndices = Array[ToSymbol["FormCalc`Lor", #]&, 10]


(* set up a directory for the Fortran code *)

ChkExist[file__] := (MkDir[DirectoryName[#]]; #)& @ StringJoin[file]

(*ChkExist[dir_, file_] := (MkDir[dir]; ToFileName[dir, file])*)


MkDir[""] = "./"

MkDir[dir_String] := ToFileName[{dir}] /; FileType[dir] === Directory

MkDir[dir_String] := Check[ToFileName[{CreateDirectory[dir]}], Abort[]]

MkDir[dirs_List] := Fold[MkDir[ToFileName[##]]&, {}, dirs]

MkDir[dirs__String] := MkDir[{dirs}]


Off[CopyFile::filex]

Options[SetupCodeDir] = {Drivers -> "drivers"}

SetupCodeDir[dir_, opt___Rule] :=
Block[ {drivers, path, files = {}},
  {drivers} = ParseOpt[SetupCodeDir, opt];
  path = SetDirectory[MkDir[dir]];
  ResetDirectory[];

  If[ FileType[drivers] === Directory,
    SetDirectory[drivers];
    files = FileNames["*", "", Infinity];
    CopyFile[#, ToFileName[path, #]]&/@
      (files = FileNames["*", "", Infinity]);
    ResetDirectory[]
  ];

  If[ FileType[$DriversDir] === Directory,
    SetDirectory[$DriversDir];
    CopyFile[#, ToFileName[path, #]]&/@
      Complement[FileNames["*", "", Infinity], files];
    ResetDirectory[]
  ];

  CopyFile[
    ToFileName[$FormCalcBin, "util.a"],
    ToFileName[path, "util.a"] ];

  path
]


(* Fortran code generation *)

Dim[i_Integer] := i


FortranNames[base_, ind___] := {
  #1 <> abridge[0, #2, Infinity],
  #1 <> abridge[0, #2, $MaxFortranName -
          StringLength[$SymbolPrefix <> #1 <> #2]]
}&[ ToString[base], ToString/@ Flatten[{ind}] ]

abridge[0, ind_, n_] :=
  Insert[ind, "_", Array[List, Min[Length[ind], n]]] /; n > 0

abridge[0, ind_, n_] := {"_", abridge[0, ind, {}, n - 1]}

abridge[_, {i___}, {j___}, 0] := {i, j}

abridge[f_, {i___, s_}, {j___}, n_] :=
  If[ StringLength[s] > 2 && !DigitQ[StringTake[s, {2}]],
    abridge[1, {i}, {StringDrop[s, {2}], j}, n + 1],
    abridge[f, {i}, {s, j}, n] ]

abridge[1, {}, ind_, n_] := abridge[0, ind, {}, n]

abridge[0, {}, ind_, _] = ind


Attributes[WriteFF] = {HoldFirst, Listable}

WriteFF[x_Symbol, array_] :=
  FFPut[x, array, Block[{x}, ToString[x]]]

WriteFF[amp_, array_] :=
  FFPut[amp, array, ToString[array] <> ToString[++modnum]]


InvDef[h_[i_], h_[j_]] := Invariant[1, i, j] -> SInvariant[k[i], k[j]]

InvDef[_[i_], _[j_]] := Invariant[-1, i, j] -> TInvariant[k[i], k[j]]


InvList[n_, r__] := {InvDef[n, #]&/@ {r}, InvList[r]}

InvList[_] = {}


ProcCheck[p_] := (
  proc = p;
  name = ToString[Map[First, p, {2}], InputForm];
  legs = Plus@@ Length/@ p;
  invs = If[ legs < 4, {},
    InvList@@ MapIndexed[Apply, Drop[Signs[p], -1]] //Flatten ];
)

ProcCheck[p_, p_] = 0

_ProcCheck := Message[WriteSquaredME::incomp]


DefModName[dir_] := (
  ModName[mod_, ext_:$CodeExt] := ModName[mod, ext] =
    ToFileName[dir, file = prefix <> mod <> ext];
  header = StringReplace[header, "%t" -> TimeStamp[]];
  Hdr[desc_] := StringReplace[header, {"%f" -> file, "%d" -> desc}];
)


FFPut[Amp[p_][amp__], array_, file_] := (
  ProcCheck[p, proc];
  FFWrite[#, array, file]&/@ {amp}
)

FFPut[_[], __] = {}

FFPut[other_, __] := (Message[WriteSquaredME::noamp, other]; Abort[])


FFWrite[0, __] = {}

FFWrite[amp_, array_, file_] :=
Block[ {ind, ff, mods},
  ind = Cases[amp, SumOver[i_, r_] :> (Dim[i] = r; i), Infinity];
  ff = FFList[amp /. unused[array] -> 0 /. fcs /. xrules /. {
    _SumOver -> 1,
    int:LoopIntegral[__] :> abbint[int] }, array];
  mods = FileSplit[ff, FortranNames[file, ind], Delayed[FFMod]];
  (Indices[#] = ind)&/@ (#[[2, 2]]&)/@ mods;
  mods
]

FFMod[ff_, {fmod_, mod_}] :=
Block[ {hh},
  hh = OpenCode[ModName[fmod]];
  WriteString[hh,
    Hdr["form factors for " <> name] <> "\
#include \"" <> prefix <> "vars.h\"\n" <>
    fincl[[2]] <> "\n\n" <>
    SubroutineDecl[mod] <> "\
#include \"" <> prefix <> "vars.h\"\n"];
  WriteExpr[hh, {sincl[[2]], ff},
    Optimize -> True, DebugLines -> fmod];
  WriteString[hh, SubroutineEnd[]];
  Close[hh];
  mod
]


FFList[0, _] = {}

FFList[amp_, array_] :=
  ( maxmat[array] = {Mat[1]};
    RuleAdd[array[1], amp] ) /; FreeQ[amp, Mat]

FFList[amp_Plus, array_] := FFList[#, array]&/@ List@@ amp

FFList[Mat[m_] x_., array_] :=
  ( maxmat[array] = MaxDims[maxmat[array], Level[m, {-2}]];
    RuleAdd[Level[m, {-1}, array], x] )


(* Calculating the abbreviations in a clever way is key to a decent
   performance of the generated code.  Therefore, the abbreviations are
   split into three categories:
   1. objects that depend only on model constants and S
      -> subroutine abbr*s,
   2. objects that depend on other phase-space variables (angles etc.)
      -> subroutine abbr*angle,
   3. objects that depend on the helicities
      -> subroutine abbr*hel.
   The master subroutine SquaredME takes care to invoke these abbr_nnn
   subroutines only when necessary. *)

Category[rul_] := {{}, {}, rul} /;
  !FreeQ[rul, Hel | e | ec | Spinor | q1 | MuTildeSq]

Category[rul_] := {{}, rul, {}} /; !FreeQ[rul, angledep]

Category[rul_] := {rul, {}, {}}


	(* caveat: a single pattern "d_" is fast but does not cover
	   multiple instances of NDim *)
VarPatt[x_] := Block[{NDim}, _NDim = d_; x]


setdep[Tag[___, r_]] := setdep[r]

setdep[v_ -> _] := v = TAG[++c]

markdep[Tag[___, r_]] := markdep[r]

markdep[v_ -> x_] := (v = TAG[++c]; {}) /; !FreeQ[x, TAG]

markdep[r_] := {} /; !FreeQ[r, TAG]

MoveDepsRight[li_List] := {li}

MoveDepsRight[f__List, li_List] :=
Block[ {pos = {f}, c = 0, cc = 0},
  Block[ #,
    setdep/@ VarPatt[li];
    pos = Map[markdep, pos, {2}];
    While[c != cc, cc = c; pos = Evaluate/@ pos]
  ]&[ Union@@ Map[Kind, {f, li}, {2}] ];
  pos = Position[pos, {}, {2}, Heads -> False];
  Append[
    MoveDepsRight@@ Delete[{f}, pos],
    Flatten[{li, Extract[{f}, pos]}] ]
]


depcat[expr_][r:Tag[___, v_ -> _]] := (cpos = 1; {r, {}}) /;
  FreeQ[expr, VarPatt[v]]

depcat[expr_][r:(v_ -> _)] := (cpos = 1; {r, {}}) /;
  FreeQ[expr, VarPatt[v]]

depcat[_][r_] := (cnew = 1; {{}, r})

MoveDepsLeft[li_List] := {li}

MoveDepsLeft[li_List, f__List] :=
Block[ {pos = {f}, old = {}, new = li, cpos = 1, cnew = 1},
  While[ cpos + cnew > 1,
    old = {new, old};
    cpos = cnew = 0;
    {pos, new} = Transpose[ ToCat[2, depcat[new]/@ #1]&/@ pos ];
    new = Flatten[new] ];
  Prepend[MoveDepsLeft@@ pos, Flatten[{new, old}]]
]


toTicket[Tag[___, x_]] := toTicket[x]

toTicket[v_ -> x_] := (
  v = Dep[v] /. HoldPattern[Pattern][a_, _] -> a;
  Ticket[Dep[v], x] )

toTicket[x_] := Ticket[x]

fromTicket[v_, x_] := v -> x

fromTicket[x_] := x

OnePassOrder::recurs = "Recursive definition among `` \
(full list in $OnePassDebug).  Returning list unordered."

OnePassOrder[li_List] :=
Block[ {c = 0, l = Length[li], Dep, Ticket, posmap, prev},
  Attributes[Dep] = {HoldFirst};
  Block[#,
    posmap = Hold@@ Evaluate[toTicket/@ VarPatt[li]]
  ]& @ Union[Kind/@ li];
  Ticket[v_, x_] := (v = Random[]; ++c) /; FreeQ[x, Dep];
  Ticket[x_] := ++c /; FreeQ[x, Dep];
  While[ c < l,
    prev = c;
    posmap = Evaluate/@ posmap;
    If[ c === prev,
      $OnePassDebug = Cases[posmap, t_Ticket :> fromTicket@@ t] /.
        Dep -> Identity;
      Message[OnePassOrder::recurs, $OnePassDebug];
      Return[li] ]
  ];
  li[[ Level[Sort[MapIndexed[List, posmap]], {3}] ]]
]

$OnePassDebug = {}


AbbrMod[{}, _] := Sequence[]

AbbrMod[abbr_, mod_] :=
Block[ {hh},
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["abbreviations for " <> name] <> "\
#include \"" <> prefix <> "vars.h\"\n" <>
    fincl[[2]] <> "\n\n" <>
    SubroutineDecl[mod] <> "\
#include \"" <> prefix <> "vars.h\"\n"];
  WriteExpr[hh, {sincl[[2]], abbr}];
  WriteString[hh, SubroutineEnd[]];
  Close[hh];
  mod
]


NumMod[{}, _] := Sequence[]

NumMod[expr_, mod_] :=
Block[ {hh},
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["numerators for " <> name] <> "\
#include \"" <> prefix <> "vars.h\"\n\
#include \"num.h\"\n"];
  WriteNum[hh, #]&/@ expr;
  Close[hh];
  mod
]


NumName[h_[___], ___] := _h -> $SymbolPrefix <> ToCode[h]

NumName[h_, ___] := h -> $SymbolPrefix <> ToCode[h]


WriteNum[hh_, var_ -> Num[expr_]] :=
Block[ {ex, name = NumName[var][[2]], $AbbPrefix = "sp"},
  ex = Subexpr[expr, MatchQ[#, _WeylChain]&,
    Deny -> {}, Fuse -> False];
  WriteString[hh, "\n\n\
\tNumeratorFunction(" <> name <> ")\n\
\timplicit none\n"];
  WriteExpr[hh, {"\
#include \"" <> prefix <> "vars.h\"\n\
#include \"num.h\"\n",
      ex[[1]],
      "Result(" <> name <> ")" -> ex[[2]]},
    Type -> "ComplexType",
    FinalTouch -> Simplify,
    FinalCollect -> True];
  WriteString[hh, "\tend\n"];
]


$SymbolPrefix = ""


VarDecl[vars__] := StringJoin @ (varArgs@@ ({vars} /. (x_ -> _) -> x))[]

varArgs[{}, _, r___] := varArgs[r]

varArgs[{l__}, r___] := varArgs[l, r] /; !FreeQ[{l}, List, {2, Infinity}]

varArgs[s_String, r___][d___] := varArgs[r][d, s]

varArgs[vars_List, type_, r___][d___] := varArgs[r][d, {
  DeleteCases[
    Replace[vars /. NDim -> Identity,
      i_Symbol :> Dim[i], {2, Infinity}] /. Dim -> Identity,
    _[0] | _[] ],
  type }]

varArgs[foo_[vars__], r___][d___] := {
  varDecl@@@ {d},
  varArgs[vars, 0, foo, r][],
  varArgs[r][] }

varArgs[0, Common[com_String], r___][d__] := {
  varDecl[Level[{d}, {3}], varDecl@@@ {d}, com],
  If[ MatchQ[{r}, {_String, ___}], {}, "\n" ] }

varArgs[0, com_NameMap, r___][d__] := {
  varDecl[{d}, com],
  If[ MatchQ[{r}, {_String, ___}], {}, "\n" ] }

varArgs[0, NotEmpty, ___][d__] := varDecl@@@ {d} /;
  !FreeQ[{d}, {__}, {2, Infinity}]

varArgs[0, ___][___] = {}

varArgs[][d___] := varDecl@@@ {d}


varDeclF[vars_, Extern] := varDeclF[vars, "external"]

varDeclF[vars_, NameMap[com_]] :=
Block[ {arr = {}},
  { varMapF[com]@@@ vars,
    varDeclF[Flatten[arr], "common /" <> $SymbolPrefix <> com <> "/"] }
]

varDeclF[vars_, type_String] :=
Block[ {Continuation},
  Format[_Continuation] = "\t";
  { "\t", type, " ", StringReplace[
    StringTake[ToString[ToCode/@ vars, OutputForm,
      PageWidth -> 63 - StringLength[type]], {2, -2}],
    ", \n \n\t" -> "\n\t" <> type ], "\n" }
]

varDeclF[vars_, decl_, com_String] := {
  decl,
  varDeclF[kind/@ vars, "common /" <> $SymbolPrefix <> com <> "/"] }


varDeclC[vars_, Extern] :=
  varDeclC[ToFunction/@ vars, "extern void"]

varDeclC[vars_, NameMap[com_]] :=
  {"struct v", com, " {\n", varMapC[com]@@@ vars, "} v", com, ";\n"}

varDeclC[vars_, type_String] :=
  { "  ", type, " ", StringReplace[
    ToSeq[Replace[vars, h_[i_, j__] :>
        Fold[#1[#2]&, h, Reverse[{i, j}]], 1],
      PageWidth -> 63 - StringLength[type]],
    ", \n" -> ";\n  " <> type ], ";\n" }

varDeclC[vars_, decl_, com_String] := {
  "struct ", com, " {\n", decl, "} ", $SymbolPrefix, com, ";\n",
  {"\n#define ", #1, " ", $SymbolPrefix, com, ".", #2}&@@@ varSubstC/@ vars,
  "\n"
}


varMapF[com_][vars_, type_] :=
Block[ {off = 1, v = com <> StringTake[type, 1]},
  arr = {arr, v};
  { varMap[v <> "(", ")\n"]/@ vars,
    "#define ", n = v <> "Len", " ", Strip[ToCode[off - 1]], "\n",
    varDeclF[{v[n]}, type] }
]

varMapC[com_][vars_, type_] :=
Block[ {off = 0, v = StringTake[type, 1], n},
  { varMap[com <> "." <> v <> "[", "]\n"]/@ vars,
    "#define ", n = com <> v <> "Len", " ", Strip[ToCode[off]], "\n",
    varDeclC[{v[n]}, type] }
]

varMap[arrL_, arrR_][var_[i___]] :=
Block[ {ind = off, stride = 1, c = 104, lhs},
  lhs = Reap[Scan[
    ( ind += stride Sow[FromCharacterCode[++c]] - stride;
      stride *= # )&, {i} ]];
  off += stride;
  { "#define ", Strip[ToCode[Level[lhs, {3}, var]]], " ",
     arrL, Strip[ToCode[ind]], arrR }
]

varMap[arrL_, arrR_][var_] :=
  {"#define ", ToCode[var], " ", arrL, ToCode[off++], arrR}


Strip[s_String] := StringReplace[s, {" " -> "", "\"" -> ""}]


varSubstC[var_[i___]] :=
  {varSeq[#1, "(", #2, ",", ")"], varCarr["-1"][##]}&[
    ToCode[var],
    Array[FromCharacterCode, Length[{i}], 105], {i} ]

varSubstC[var_] := {#, #}& @ ToCode[var]


varSeq[s1___, {i___, j_}, c_, s2___] := s1 <> ({#, c}&/@ {i}) <> j <> s2

varCarr[x___][h_, i__] := h <> Reverse[MapThread[varCind[x], {i}]]

varCind["-1"][i_, n_Symbol] := varCind[][i, n] /;
  Context[n] === "LoopTools`"

varCind[x___][i_, ___] := {"[", i, x, "]"}


ToFunction[h:_[___]] := h

ToFunction[h_] := h[]


subroutineDeclF[name_, decl___String] :=
  "\tsubroutine " <> $SymbolPrefix <> ToCode[name] <>
    "\n\timplicit none\n" <> decl <> "\n"

subroutineDeclC[name_, decl___String] :=
  "void " <> $SymbolPrefix <> ToCode[ToFunction[name]] <> " {\n" <>
    decl <> "\n"


CallDecl[li_List] := StringJoin[CallDecl/@ li]

CallDecl[DoLoop[name_, ind__]] :=
  ("\n" <> #1 <> CallDecl[name] <> #2 &)@@ DoDecl[ind]

CallDecl[name_] := callDecl[name]

callDeclC[name_] :=  "  " <> $SymbolPrefix <> name <> "();\n"

callDeclF[name_] := "\tcall " <> $SymbolPrefix <> name <> "\n"


DoDecl[{var_}] := DoDecl[{var, Dim[var]}]

DoDecl[{_, _Dim}] := {{}, {}}

DoDecl[{var_, from_:1, to_, step_:1}] := {
  "!LOOP(" <> ToCode[var] <> ", " <> ToCode[{from, to, step}] <> ")\n",
  "!ENDLOOP(" <> ToCode[var] <> ")\n" }

DoDecl[var_] := DoDecl[{var, Dim[var]}]

DoDecl[vars__] := {StringJoin[#1], StringJoin[Reverse[#2]]}&@@
  Transpose[DoDecl/@ ReverseDo[{vars}]]

ReverseDo = (*Identity*) Reverse


	(* LoopElemF gives back e.g.
1. {F[jFtree], SUN[jSUNtree]}
2. "tree(jFtree,jSUNtree)"
3. "tree(nFtree,nSUNtree)"
4. "nFtree,nSUNtree,"
5. "\n\tinteger jFtree, nFtree"
6. "\n\tinteger jFtree, nFtree
    \n\tparameter (nFtree = 5)\n"
7. "\n\tLOOP(jFtree, 1,nFtree,1)
    \n\tLOOP(jSUNtree, 1,nSUNtree,1)"
8. "\n\tENDLOOP(jSUNtree)
    \n\tENDLOOP(jFtree)" *)

LoopVarF[_][h_[-1]] = {h[1], "1", "1", "", "", "", "", ""}

LoopVarF[type_][h_[n_]] :=
Block[ {jv = "j" <> ToString[h] <> type, nv = "n" <> ToString[h] <> type, s},
  { h[jv], jv, nv,
    nv <> ",",
    s = "\n\tinteger " <> jv <> ", " <> nv,
    s <> "\n\tparameter (" <> nv <> " = " <> ToString[n] <> ")\n",
    "\n\tLOOP(" <> jv <> ", 1," <> nv <> ",1)",
    "\n\tENDLOOP(" <> jv <> ")" }
]


LoopElemF[_, {}] = ""

LoopElemF[type_, maxmat_] :=
  { #1, varSeq[type, "(", #2, ",", ")"],
    varSeq[type, "(", #3, ",", ")"], ##4 }&@@
  Transpose[LoopVarF[type]/@ maxmat]


	(* LoopElemC gives back e.g.
1. {F[jFtree], SUN[jSUNtree]}
2. "tree(jFtree,jSUNtree)"
3. "tree[nSUNtree][nFtree]"
4. "tree[jSUNtree-1][jFtree-1]"
5. "nFtree,nSUNtree,"
6. "int nFtree, int nSUNtree"
7. "\n  int jFtree;"
8. "\n  int jFtree;
    \n  enum { nFtree = 5 };\n"
9. "\n  LOOP(jFtree, 1,nFtree,1)
    \n  LOOP(jSUNtree, 1,nSUNtree,1)"
10. "\n  ENDLOOP(jSUNtree)
     \n  ENDLOOP(jFtree)" *)

LoopVarC[_][h_[-1]] = {h[1], "1", "1", "", "", "", "", "", ""}

LoopVarC[type_][h_[n_]] :=
Block[ {jv = "j" <> ToString[h] <> type, nv = "n" <> ToString[h] <> type, s},
  { h[jv], jv, nv,
    nv <> ",",
    "int " <> nv <> ", ",
    s = "\n  int " <> jv <> ";",
    s <> "\n  enum { " <> nv <> " = " <> ToString[n] <> " };\n",
    "\n  LOOP(" <> jv <> ", 1," <> nv <> ",1)",
    "\n  ENDLOOP(" <> jv <> ")" }
]


LoopElemC[_, {}] = ""

LoopElemC[type_, maxmat_] :=
  { #1, varSeq[type, "(", #2, ",", ")"], varCarr[][type, #3], 
    varCarr["-1"][type, #2], ##4 }&@@
  Transpose[LoopVarC[type]/@ maxmat]


Invoke[mod0_, {}] := CallDecl[mod0]

Invoke[mod0_, mod1_] := {
  CallDecl[mod0],
  "!TEST(flags, BIT_LOOP)\n",
  CallDecl[mod1],
  "!ENDTEST(flags, BIT_LOOP)\n" }


LoopReduce[m_] := Transpose[MapThread[LoopNeed, m]] /; SameQ@@ Length/@ m

LoopReduce[m_] := m /. 1 -> -1

LoopNeed[h_[1], h_[1]] := {h[-1], h[-1]}

LoopNeed[other__] := {other}


MatType[_Mat, _] = 1

MatType[h_[i_], h_[j_]] := MatType[h][i, j]

MatType[h_] := MatType[h] = ToSymbol["Mat", h]


Assort[m_Mat -> x_] := {m -> x, {}, {}, {}}

Assort[f_ -> x_] := {{}, {}, {}, f -> x} /;
  !FreeQ[x, WeylChain]

Assort[f_ -> x_] := {{}, f -> ToArray[f], {}, {}} /;
  !FreeQ[x, DiracChain | SUNT]

Assort[x_ -> n_Num] := {{}, {}, x -> n, {}}

Assort[other_] := {{}, {}, {}, other}


ChainHead[om_, n_] := ChainHead[om, n] =
  ToSymbol["Chain", {"V", "B"}[[om - 5]], n]

SplitChain[Spinor[_[i_], _, s1_, d1_, e1_], om_, g___,
    Spinor[_[j_], _, s2_, d2_, e2_]] :=
  afsign[i, s1, OddQ[om - 1 - e1]] *
  afsign[j, s2, OddQ[om + Length[{g}] + e2]] *
  ChainHead[om, Length[{g}]][Spinor[i, s1, d1], e1, g,
    e2, Spinor[j, s2, d2]]


afsign[i_, -1, False] = -Hel[i]

afsign[i_, -1, True] = Hel[i]

_afsign = 1


DefFilter[h_][m_[i_, j_]] :=
  h[_[m[x_, y_], _]] := x <= i && y <= j

DefCat[h_][m_[i_, j_]] :=
  h[r:_[m[x_, y_], _]] := If[x <= i && y <= j, {r, {}}, {{}, r}]


Attributes[Cond] = {HoldRest}

Cond[True, args__] := {args}

Cond[False, __] = {}


abbrtype["0h"] = abbrtype["1h"] = "HelType"

_abbrtype = "ComplexType"


Attributes[WriteSquaredME] = {HoldAll}

Options[WriteSquaredME] = {
  TreeSquare -> True,
  LoopSquare -> False,
  Folder -> "squaredme",
  ExtraRules -> {},
  FilePrefix -> "",
  SymbolPrefix -> "",
  FileHeader -> "#if 0\n* %f\n* %d\n* generated by FormCalc " <>
    ToString[NumberForm[$FormCalc, {Infinity, 1}]] <> " on %t\n#endif\n\n",
  FileIncludes -> {"#include \"decl.h\"\n", "#include \"inline.h\"\n"},
  SubroutineIncludes -> FileIncludes
}

WriteSquaredME::noamp = "`` is not an amplitude."

WriteSquaredME::empty = "Warning: no amplitudes were specified."

WriteSquaredME::badmat = "Incompatible matrix elements `` and ``."

WriteSquaredME[tree_, loop_, dir_, opt___Rule] :=
  WriteSquaredME[tree, loop, Abbr[], dir, opt]

WriteSquaredME[tree_, loop_, abbr__, dir_, opt___Rule] :=
Block[ {treesq, loopsq, folder, xrules, prefix, $SymbolPrefix,
header, fincl, sincl,
ModName, Hdr, proc = Sequence[], name, legs, invs,
mat, nums, fcs, abrs, matsel, treecat, angledep, abbrsel,
Dim, abbint, lint = {}, pavec = 0, cutc = 0, masc = 0, defs,
Indices, Hel, hels, pos, file, files, hh,
unused, maxmat, ntree, nloop, mats, com, loops,
ffmods, nummods, abbrmods, abbrmap},

  {treesq, loopsq, folder, xrules, prefix, $SymbolPrefix,
    header, fincl, sincl} =
    ParseOpt[WriteSquaredME, opt] //. Options[WriteRenConst];

  fincl = Flatten[{fincl, "", ""}];
  sincl = Flatten[{sincl, "", ""}];

  DefModName[MkDir[dir, folder]];

  abbint[f_] := AbbrevInt[f];

  {mat, fcs, nums, abrs} = ToCat[4, Assort/@ Flatten[{abbr}]];

  mats = First/@ DeleteCases[mat, _ -> 0];
  unused[Ctree] = Alt@@
    Select[Union[#[[1, 2]]&/@ mat], FreeQ[mats, Mat[_, #]]&];
  unused[Cloop] = Alt@@
    Select[Union[#[[1, 1]]&/@ mat], FreeQ[mats, Mat[#, _]]&];

(* Part 1: the form factors *)

  _maxmat = {};
  ffmods = Flatten/@ {
    Block[{modnum = 0}, WriteFF[tree, Ctree]],
    Block[{modnum = 0}, WriteFF[loop, Cloop]] };
  If[ Plus@@ Length/@ ffmods === 0,
    Message[WriteSquaredME::empty]; Return[{}] ];

  abrs = abrs /. xrules /. int:LoopIntegral[__] :> abbint[int];
  lint = Flatten[lint];

  ntree = maxmat[Ctree];
  nloop = maxmat[Cloop];
  If[ Length[ntree] === 0, ntree = nloop,
    If[ Length[nloop] === 0, nloop = ntree,
      If[ Head/@ ntree =!= Head/@ nloop,
        Message[WriteSquaredME::badmat, ntree, nloop];
        Abort[] ];
      nloop = MaxDims[nloop, ntree];
      If[ loopsq, ntree = nloop ];
  ] ];
  mats = Select[MapThread[MatType, {nloop, ntree}], Length[#] > 0 &];

(* Part 2: the numerators and abbreviations *)

  Scan[DefFilter[matsel], mats];
  mat = Select[mat /. fcs /. Mat -> MatType, matsel];

  defs = Flatten[{abrs, mat, nums, lint}];

  (* split into tree/loop *)
  Scan[DefCat[treecat][MatType[#, #]]&, maxmat[Ctree]];
  treecat[r:_[v_, _]] := {r, {}} /; !FreeQ[tree, v];
  treecat[r_] := {{}, r};
  defs = Join[#1, Tag/@ #2]&@@
    MoveDepsLeft@@ ToCat[2, treecat/@ defs];

  (* split into s/angle/hel *)
  angledep = Alt[(Range[#2] + #1)&@@ Length/@ proc];
  angledep = Alt[{
    Cases[invs, _[x_, r_] :> x /; MemberQ[r, angledep, {2}]],
    (k | s)[angledep],
    MomEncoding[_, angledep] }];
  defs = OnePassOrder/@ MoveDepsRight@@ ToCat[3, Category/@ defs];
  pos = Take[#, 2]&/@ Position[defs, _Num];
  nums = Extract[defs, pos] /. Tag -> Identity;
  defs = ToCat[2, #]&/@ Replace[ Delete[defs, pos],
    {Tag[r_] -> {{}, r}, r_ -> {r, {}}}, {2} ];

  nummods = FileSplit[nums, "num", Delayed[NumMod]];
  nums = NumName@@@ nums;

  abbrmap[foo_] := MapThread[foo,
    { defs, {{"0s", "1s"},
             {"0a", "1a"},
             {"0h", "1h"}} }, 2];

  abbrmods = abbrmap[FileSplit[
    ToDoLoops[#1],
    "abbr" <> #2,
    Delayed[AbbrMod] ]&] /. nums;

  abrs = Join[abrs, lint];
  abbrsel[sel_] := Select[abrs, MemberQ[sel, #]&];

  com = VarDecl[
    abbrmap[Common["var" <> #2][abbrsel[#1], abbrtype[#2]]&],
    Common["kinvars"][invs, "RealType"],
    Common["helvars"][{"seq"[2], Hel[legs]}, "integer"],
    Common["indices"][#[[1, 1, 1]]&/@ DownValues[Dim], "integer"],
    Common["formfactors"][
      Flatten[{
        Level[maxmat[Ctree], {2}, Ctree],
        Level[maxmat[Cloop], {2}, Cloop], mats }], "ComplexType" ],
    NotEmpty["#ifndef NUM_H\n",
      Last/@ nums, Extern,
      "#endif\n\n"] ];

  {ffmods, abbrmods, nummods} = {ffmods, abbrmods, nummods} /.
    Delayed -> Identity;

(* Part 3: the variable declarations *)

  hh = OpenCode[ModName["vars", ".h"]];

  WriteString[hh, Hdr["variable declarations"] <> "\
#ifndef VARS_H\n\
#define VARS_H\n\n\
#define LEGS " <> ToString[legs] <> "\n\n" <>
    fincl[[1]] <> "\n" <> If[$Code === "C", com, {"\
#ifdef PARALLEL\n\
#define ", $SymbolPrefix, "var1s ", $SymbolPrefix, "var0s\n\
#define ", $SymbolPrefix, "var1a ", $SymbolPrefix, "var0a\n\
#define ", $SymbolPrefix, "kinvars ", $SymbolPrefix, "var0a\n\
#endif\n\n"}] <>
    "#else\n\n" <>
    sincl[[1]] <> "\n" <> If[$Code === "C", "", {"\
#ifdef PARALLEL\n\
\tmarker b0s(0), b0a(0), bhel(0)\n\
\tcommon /", $SymbolPrefix, "var0s/ b0s\n\
\tcommon /", $SymbolPrefix, "var0a/ b0a\n\
\tcommon /", $SymbolPrefix, "helvars/ bhel\n\
#endif\n\n\
", com, "\
#ifdef PARALLEL\n\
\tmarker e0s(0), e0a(0), ehel(0)\n\
\tcommon /", $SymbolPrefix, "var0s/ e0s\n\
\tcommon /", $SymbolPrefix, "var0a/ e0a\n\
\tcommon /", $SymbolPrefix, "helvars/ ehel\n\
#endif\n\n"}] <>
    "#endif\n"];

  Close[hh];

(* Part 4: the makefile *)

  hh = OpenWrite[ModName["makefile", ""]];

  WriteString[hh, "\
NUMS :=" <> ({" \\\n  $(DIR)/", #, ".o"}&)/@
    nummods <> "\n\n\
OBJS := $(NUMS)" <> ({" \\\n  $(DIR)/", #, ".o"}&)/@
    Flatten[{abbrmods, ffmods, prefix <> "SquaredME"}] <> "\n\n\
$(LIB): $(LIB)($(OBJS))\n\n\
$(LIB)($(OBJS)): $(DIR)/" <> prefix <> "vars.h " <> $MakeDeps[[1]] <> "\n\n\
$(LIB)($(NUMS)): " <> $MakeDeps[[2]] <> "\n\n\
LIBS += $(LIB)\n\n"];

  Close[hh];

(* Part 5: the master subroutine SquaredME *)

  hels = Array["Hel" <> ToString[#] -> ToString[#]&, legs];
  {maxmat[Ctree], maxmat[Cloop]} = LoopReduce[{maxmat[Ctree], maxmat[Cloop]}];

  ntree = stree = loopElem["tree", maxmat[Ctree]];
  nloop = sloop = loopElem["loop", maxmat[Cloop]];
  loops = DeleteCases[{ntree, nloop}, ""];
  If[ stree === "", stree = loopElem["tree", maxmat[Cloop]] ];
  If[ sloop === "", sloop = loopElem["loop", maxmat[Ctree]] ];

  hh = OpenCode[ModName["SquaredME"]];
  writeSquaredME[hh];
  Close[hh];

  Cases[DownValues[ModName], _[_, s_String] -> s]
]


writeSquaredMEF[hh_] := (
  (* a) declarations *)
  WriteString[hh, "\
*#define CHECK\n\n" <>
    Hdr["assembly of squared matrix element"] <> "\
#include \"" <> prefix <> "vars.h\"\n" <>
    fincl[[2]] <> "\n\n\
************************************************************************\n\n\
\tRealType function " <> $SymbolPrefix <> "sumup(" <>
    sloop[[4]] <> "CCloop, " <> stree[[4]] <> "CCtree)\n\
\timplicit none\n\n\
#include \"" <> prefix <> "vars.h\"\n" <>
    sloop[[5]] <>
    stree[[5]] <> "\n\
\tComplexType CC" <> sloop[[3]] <> ", CC" <> stree[[3]] <> "\n\
\tComplexType m\n\n\
\t" <> $SymbolPrefix <> "sumup = 0\n" <>
    stree[[7]] <> "\n\
\tm = 0" <>
    sloop[[7]] <> "\n\
\tm = m + CC" <> sloop[[2]] <> "*" <>
      ToCode[Inner[MatType, sloop[[1]], stree[[1]], Times]] <>
    sloop[[8]] <> "\n\
\t" <> $SymbolPrefix <> "sumup = " <>
      $SymbolPrefix <> "sumup + Re(Conjugate(CC" <> stree[[2]] <> ")*m)" <>
    stree[[8]] <> "\n\
\tend\n
************************************************************************\n\n\
\tsubroutine " <> $SymbolPrefix <> "SquaredMEHel(result, flags)\n\
\timplicit none\n\
\tRealType result(*)\n\
\tinteger flags\n\n\
#include \"" <> prefix <> "vars.h\"\n\n\
\tRealType " <> $SymbolPrefix <> "sumup\n\
\texternal " <> $SymbolPrefix <> "sumup\n" <>
    (#6 &)@@@ loops <> "\n\
* " <> prefix <> "BEGIN ABBR_HEL\n" <>
    Invoke@@ abbrmods[[3]] <> "\
* " <> prefix <> "END ABBR_HEL\n\n\
* " <> prefix <> "BEGIN FF_INI" <>
    ({#7, "\n\tC", #2, " = 0", #8, "\n"}&)@@@ loops <> "\
* " <> prefix <> "END FF_INI\n" <>
  ({"\n\
* ", prefix, "BEGIN FF_TREE\n",
    CallDecl[ToDoLoops[ffmods[[1]], Indices]],
    Cond[ TrueQ[treesq],
      "\n\tresult(1) = result(1) + ", $SymbolPrefix, "sumup(",
        #4, "Ctree, ", #4, "Ctree)" ], "\n\
* ", prefix, "END FF_TREE\n"}&)@@ ntree <>
  ({"\n\
\tTEST(flags, BIT_LOOP)\n\
* ", prefix, "BEGIN FF_LOOP\n",
    CallDecl[ToDoLoops[ffmods[[2]], Indices]],
    Cond[ ntree =!= "",
      "\n\tresult(2) = result(2) + 2*", $SymbolPrefix, "sumup(",
        #4, "Cloop, ", ntree[[4]], "Ctree)" ],
    Cond[ ntree === "" || TrueQ[loopsq],
      "\n\tresult(2) = result(2) + ", $SymbolPrefix, "sumup(",
        #4, "Cloop, ", #4, "Cloop)" ], "\n\
* ",  prefix, "END FF_LOOP\n\
\tENDTEST(flags, BIT_LOOP)\n"}&)@@ nloop <> "\
\tend\n\n\
************************************************************************\n\n\
\tsubroutine " <> $SymbolPrefix <> "SquaredME(result, helicities, flags)\n\
\timplicit none\n\
\tRealType result(*)\n\
\tinteger*8 helicities\n\
\tinteger flags\n\n\
#include \"" <> prefix <> "vars.h\"\n\n\
* " <> prefix <> "BEGIN VARDECL\n\
\texternal " <> $SymbolPrefix <> "SquaredMEHel\n\n" <>
    VarDecl[hels, "integer"] <>
    ({"\tequivalence (Hel(", #2, "), ", #1, ")\n"}&)@@@ hels <> "\n" <>
    ({"\tdata ", ToString[Head[#]],
      " /", ToString[Times@@ #], "*bogus/\n"}&)/@ mats <> "\
* " <> prefix <> "END VARDECL\n\n" <>
    sincl[[2]] <> "\n\
\tPREP(bhel,ehel, vec,vec_end, b0a,e0a, b0s,e0s)\n"];

  (* b) definitions of the invariants *)
  WriteExpr[hh, {"\
* " <> prefix <> "BEGIN INVARIANTS\n",
    invs, "\
* " <> prefix <> "END INVARIANTS\n\n"}, Newline -> ""];

  (* c) calculation of the abbreviations *)
  WriteString[hh, "\
\tTEST(flags, BIT_RESET)\n\
* " <> prefix <> "BEGIN ABBR_S\n\
\tseq(1) = seq(1) + 1\n\
\tINI_S(seq)\n" <>
    Invoke@@ abbrmods[[1]] <> "\
* " <> prefix <> "END ABBR_S\n\
\tENDTEST(flags, BIT_RESET)\n\n\
* " <> prefix <> "BEGIN ABBR_ANGLE\n\
\tseq(2) = seq(2) + 1\n\
\tINI_ANGLE(seq)\n" <>
    Invoke@@ abbrmods[[2]] <> "\
* " <> prefix <> "END ABBR_ANGLE\n\n\
* " <> prefix <> "BEGIN RES_INI\n\
\tresult(1) = 0\n\
\tresult(2) = 0\n\
* " <> prefix <> "END RES_INI\n\n\
* " <> prefix <> "BEGIN HEL_LOOPS\n" <>
    ({"\tLOOP_HEL(", #1, ")\n\
\tTEST(helicities, BIT_HEL(", #2, "))\n\n"}&)@@@ hels <> "\
\tEXEC(" <> $SymbolPrefix <> "SquaredMEHel, result, flags)\n" <>
    ({"\n\tENDTEST(helicities, BIT_HEL(", #2, "))\n\
\tENDLOOP_HEL(", #1, ")\n"}&)@@@ Reverse[hels] <> "\
* " <> prefix <> "END HEL_LOOPS\n\n\
\tSYNC(result)\n\
\tDEINI(seq)\n\n\
#ifdef CHECK" <>
    ({"\n\tprint *, '", #, " =', ", #}&)/@
      (ToString[#1]&)@@@ invs <> "\n\
\tprint *, 'tree =', result(1)\n\
\tprint *, 'loop =', result(2)\n\
\tstop\n\
#endif\n\n\
* " <> prefix <> "END SQUAREDME\n\
\tend\n\n\
"];
)


writeSquaredMEC[hh_] := (
  (* a) declarations *)
  WriteString[hh, "\
//#define CHECK\n\n" <>
    Hdr["assembly of squared matrix element"] <> "\
#include <math.h>\n\
#include \"" <> prefix <> "vars.h\"\n" <>
    fincl[[2]] <> "\n\n\
#if NOUNDERSCORE\n\
#define " <> $SymbolPrefix <> "SquaredME " <> $SymbolPrefix <> "squaredme\n\
#else\n\
#define " <> $SymbolPrefix <> "SquaredME " <> $SymbolPrefix <> "squaredme_\n\
#endif\n\n" <>
    ({"void ", $SymbolPrefix, #, "(void);\n"}&)/@
      Flatten[{abbrmods, ffmods}] <>
    "\n" <>
    Cond[ Length[mats] > 0,
      "struct formfactors " <> $SymbolPrefix <> "formfactors = {\n",
      varSeq[{"  .", ToString[Head[#]],
       {"[0 ... ", ToString[# - 1], "]"}&/@ List@@ Reverse[#],
       " = NAN"}&/@ mats, ",\n"], " };\n\n" ] <> "\
/**********************************************************************/\n\n\
static RealType " <> $SymbolPrefix <> "sumup(" <>
      sloop[[6]] <> "ComplexType CC" <> sloop[[3]] <> ",\n  " <>
      stree[[6]] <> "ComplexType CC" <> stree[[3]] <> ") {\n\n\
#include \"" <> prefix <> "vars.h\"\n" <>
    sloop[[7]] <>
    stree[[7]] <> "\n\
  ComplexType s = 0;\n" <>
    stree[[9]] <> "\n\
  ComplexType m = 0;" <>
    sloop[[9]] <> "\n\
  m += CC" <> sloop[[4]] <> "*" <>
      ToCode[Inner[MatType, sloop[[1]], stree[[1]], Times]] <> ";" <>
    sloop[[10]] <> "\n\
  s += Re(Conjugate(CC" <> stree[[4]] <> ")*m);" <>
    stree[[10]] <> "\n\
  return s;\n\
}\n\n\
/**********************************************************************/\n\n\
void " <> $SymbolPrefix <> "SquaredMEHel(RealType *result, int *flags) {\n\n\
#include \"" <> prefix <> "vars.h\"\n" <>
    (#8 &)@@@ loops <> "\n\
// " <> prefix <> "BEGIN ABBR_HEL\n" <>
    Invoke@@ abbrmods[[3]] <> "\
// " <> prefix <> "END ABBR_HEL\n\n\
// " <> prefix <> "BEGIN FF_INI" <>
    ({#9, "\n  C", #2, " = 0;", #10, "\n"}&)@@@ loops <> "\
// " <> prefix <> "END FF_INI\n" <>
    ({"\n\
// ", prefix, "BEGIN FF_TREE\n",
      CallDecl[ToDoLoops[ffmods[[1]], Indices]],
      Cond[ TrueQ[treesq],
        "\n  result[0] += ", $SymbolPrefix, "sumup(",
          #5, $SymbolPrefix, "formfactors.Ctree, ",
          #5, $SymbolPrefix, "formfactors.Ctree);" ], "\n\
// ", prefix, "END FF_TREE\n"}&)@@ ntree <>
    ({"\n\
  TEST(flags, BIT_LOOP)\n\
// ", prefix, "BEGIN FF_LOOP\n",
      CallDecl[ToDoLoops[ffmods[[2]], Indices]],
      Cond[ ntree =!= "",
        "\n  result[1] += 2*", $SymbolPrefix, "sumup(",
          #5, $SymbolPrefix, "formfactors.Cloop, ",
          ntree[[5]], $SymbolPrefix, "formfactors.Ctree);" ],
      Cond[ ntree === "" || TrueQ[loopsq],
        "\n  result[1] += ", $SymbolPrefix, "sumup(",
          #5, $SymbolPrefix, "formfactors.Cloop, ",
          #5, $SymbolPrefix, "formfactors.Cloop);" ], "\n\
// ",  prefix, "END FF_LOOP\n\
  ENDTEST(flags, BIT_LOOP)\n"}&)@@ nloop <> "\
}\n\n\
/**********************************************************************/\n\n\
void " <> $SymbolPrefix <> "SquaredME(RealType *result, long long int *helicities, int *flags) {\n\n\
#include \"" <> prefix <> "vars.h\"\n\n\
// " <> prefix <> "BEGIN VARDECL\n" <>
    VarDecl[hels, "int"] <> "\
// " <> prefix <> "END VARDECL\n\n" <>
    sincl[[2]] <> "\n\
  PREP(seq,seq_end, vec,vec_end, varangle,varangle_end, vars,vars_end);\n"];

  (* b) definitions of the invariants *)
  WriteExpr[hh, {"\
// " <> prefix <> "BEGIN INVARIANTS\n",
    invs, "\
// " <> prefix <> "END INVARIANTS\n\n"}, Newline -> ""];

  (* c) calculation of the abbreviations *)
  WriteString[hh, "\
  TEST(flags, BIT_RESET)\n\
// " <> prefix <> "BEGIN ABBR_S\n\
  ++seq(1);\n\
  INI_S(seq);\n" <>
    Invoke@@ abbrmods[[1]] <> "\
// " <> prefix <> "END ABBR_S\n\
  ENDTEST(flags, BIT_RESET)\n\n\
// " <> prefix <> "BEGIN ABBR_ANGLE\n\
  ++seq(2);\n\
  INI_ANGLE(seq);\n" <>
    Invoke@@ abbrmods[[2]] <> "\
// " <> prefix <> "END ABBR_ANGLE\n\n\
// " <> prefix <> "BEGIN RES_INI\n\
  result[0] = 0;\n\
  result[1] = 0;\n\
// " <> prefix <> "END RES_INI\n\n\
// " <> prefix <> "BEGIN HEL_LOOPS\n" <>
    ({"  LOOP_HEL(Hel(", #2, "))\n\
  TEST(helicities, BIT_HEL(", #2, "))\n\n"}&)@@@ hels <> "\
  EXEC(" <> $SymbolPrefix <> "SquaredMEHel, result, flags);\n" <>
    ({"\n  ENDTEST(helicities, BIT_HEL(", #2, "))\n\
  ENDLOOP_HEL(", #1, ")\n"}&)@@@ Reverse[hels] <> "\
// " <> prefix <> "END HEL_LOOPS\n\n\
  SYNC(result);\n\
  DEINI(seq);\n\n\
#ifdef CHECK\n" <>
    ({"!printf(\"", #, " = %g\\n\", ", #, ");\n"}&)/@
      (ToString[#1]&)@@@ invs <> "\
!printf(\"tree = (%g,%g)\\n\",  Re(result[0]), Im(result[0]));\n\
!printf(\"loop = (%g,%g)\\n\",  Re(result[1]), Im(result[1]));\n\
!exit(1);\n\
#endif\n\n\
// " <> prefix <> "END SQUAREDME\n\
}\n"];
)


(* renormalization constants *)

rcpatt[_[_[_[h_[___]]], _]] := HoldPattern[_h]

rcpatt[_[_[_[h_]], _]] := HoldPattern[h]

RCPattern[other___] := Alt[{other,
  Union[rcpatt/@ DownValues[RenConst]]}]


Attributes[WithOpt] = {HoldFirst}

WithOpt[foo_, {}] := (Needs["FeynArts`"]; foo)

WithOpt[foo_, {opt__}] := (
  Needs["FeynArts`"];
  (Options[InsertFields] = #1; #2)&[
    Options[InsertFields],
    SetOptions[InsertFields, opt]; foo ] )

WithOpt[foo_, opt__] := WithOpt[foo, Flatten[{opt}]]


(* These are special versions of Re and Im where the real and
   imaginary part is taken only of the loop integrals,
   see A. Denner, Forts. Phys. 41 (1993) 307, arXiv:0709.1075. *)

ReTilde[expr_] := expr /. int:LoopIntegral[__] :> Re[int]

ImTilde[expr_] := (expr /. int:LoopIntegral[__] :> Im[int]) -
  (expr /. LoopIntegral[__] -> 0)


	(* Note: it seems weird that the left-handed vector component
	   is taken as the coefficient of DiracChain[6, k]: this is
	   because DiracChain[6, k] = DiracChain[k, 7]. *)

DiracCoeff[expr_, g__] :=
Block[ {e, ec},
  _e = _ec = Sequence[];
  (((# /. DiracChain[g] -> 1) - #) /. _DiracChain -> 0)& @ expr
]

LVectorCoeff[se_] := DiracCoeff[se, 6, _]

RVectorCoeff[se_] := DiracCoeff[se, 7, _]

LScalarCoeff[se_] := DiracCoeff[se, 7]

RScalarCoeff[se_] := DiracCoeff[se, 6]


SEPart[f_, se_] := f[se]


Attributes[OnlyIf] = {HoldRest}

OnlyIf[True, a_, _] := a

OnlyIf[False, _, b_] := b

OnlyIf[other__] := OnlyIfEval[other]

OnlyIfEval[_, a_, a_] := a

OnlyIfEval[cond_, a_, b_] := Thread @ IndexIf[ cond,
  a /. Cases[{cond}, i_ == j_ -> (IndexDelta[i, j] -> 1), Infinity],
  b /. Cases[{cond}, i_ == j_ -> (IndexDelta[i, j] -> 0)] ]


FermionRC[m_, a_, b_, se_] :=
  m (a SEPart[LVectorCoeff, se] + b SEPart[RVectorCoeff, se]) +
     b SEPart[LScalarCoeff, se] + a SEPart[RScalarCoeff, se]

BosonRC[se_] := SEPart[Identity, se]


MassRC[f_F, opt___Rule] :=
  FermionRC[TheMass[f], 1/2, 1/2, ReTilde[SelfEnergy[f, opt]]]

MassRC[f_, opt___Rule] := BosonRC[ReTilde[SelfEnergy[f, opt]]]

MassRC[f_, f_, opt___Rule] := MassRC[f, opt]

MassRC[f1_, f2_, opt___Rule] :=
  1/2 (BosonRC[ReTilde[SelfEnergy[f2 -> f1, TheMass[f1], opt]]] +
       BosonRC[ReTilde[SelfEnergy[f2 -> f1, TheMass[f2], opt]]])


FieldRC[f_F, opt___Rule] :=
Block[ {m = TheMass[f], se},
  se = ReTilde[SelfEnergy[f, opt]];
  -{SEPart[LVectorCoeff, se], SEPart[RVectorCoeff, se]} -
    FermionRC[m, m, m, ReTilde[DSelfEnergy[f]]]
]

FieldRC[f_, opt___Rule] := -BosonRC[ReTilde[DSelfEnergy[f, opt]]]

FieldRC[f_, f_, opt___Rule] := FieldRC[f, opt]

FieldRC[f1_, f2_, opt___Rule] := FieldRC[f1, f2, 0, opt]

FieldRC[f1_, f2_, c_, opt___Rule] := OnlyIf[
  And@@ MapThread[Equal, Flatten[{##}]&@@@ {f1, f2}],
  FieldRC[f1, opt],
  FieldRC2[f1, f2, c, opt]
]

FieldRC2[f1_F, f2_F, c_, opt___Rule] :=
Block[ {m1 = TheMass[f1], m2 = TheMass[f2]},
  2/(m1^2 - m2^2) (FermionRC[m2, {m2, m1}, {m1, m2},
    ReTilde[SelfEnergy[f2 -> f1, m2, opt]]] - c)
]

FieldRC2[f1_, f2_, c_, opt___Rule] :=
Block[ {m1 = TheMass[f1], m2 = TheMass[f2]},
  2/(m1^2 - m2^2) (BosonRC[ReTilde[SelfEnergy[f2 -> f1, m2, opt]]] - c)
]


TadpoleRC[f_, opt___Rule] :=
  -BosonRC[ReTilde[SelfEnergy[f -> {}, Indeterminate, opt]]]


WidthRC[f_F, opt___Rule] :=
  FermionRC[TheMass[f], 1, 1, ImTilde[SelfEnergy[f, opt]]]

WidthRC[f_, opt___Rule] :=
  BosonRC[ImTilde[SelfEnergy[f, opt]]]/TheMass[f]


Attributes[TreeCoupling] = {HoldRest}

TreeCoupling[proc_Rule, opt___Rule] :=
  WithOpt[CalcProcess[0, proc, Options[InsertFields]], {opt}]


Attributes[SelfEnergy] = {HoldRest}

SelfEnergy[proc_Rule, m_, opt___Rule] := (# /. K2 -> m^2)& @
  WithOpt[CalcSelfEnergy[proc, Options[InsertFields]], {opt}]

SelfEnergy[f_, opt___Rule] := SelfEnergy[f -> f, TheMass[f], opt]

SelfEnergy[f_, m__] := SelfEnergy[f -> f, m]


Attributes[DSelfEnergy] = {HoldRest}

DSelfEnergy[proc_Rule, m_, opt___Rule] := (# /. K2 -> m^2)& @
  WithOpt[D[CalcSelfEnergy[proc, Options[InsertFields]], K2], {opt}]

DSelfEnergy[f_, opt___Rule] := DSelfEnergy[f -> f, TheMass[f], opt]

DSelfEnergy[f_, m__] := DSelfEnergy[f -> f, m]


RCSub = Simplify

RCInt = Simplify


CalcProcess[loop_, proc_, opt_] :=
Block[ {amp, Neglect, FormSub = RCSub},
  ClearProcess[];
  amp = InsertFieldsHook[
    CreateTopologies[loop, Length[Flatten[{#}]]&/@ proc,
      ExcludeTopologies -> Internal],
    proc ];
  OptPaint[amp];
  amp = CreateFeynAmpHook[amp, Truncated -> !FreeQ[proc, F]];
  Plus@@ CalcFeynAmp[amp, OnShell -> False, Transverse -> False,
           FermionChains -> Chiral, OPP -> False, FileTag -> "rc"] //.
    Dispatch[Abbr[]] /. Mat -> Identity
]


CalcSelfEnergy[proc_, opt_] := CalcSelfEnergy[proc, opt] =
  CalcProcess[1, proc, opt] /. {
    Pair[_k, _k] -> K2,
    Pair[_e | _ec, _k] -> If[ MatchQ[proc, _V -> _V],
	(* default: take only the transverse part of vector-boson SEs *)
      If[TrueQ[$LongitudinalSE], I Sqrt[K2], 0],
      1 ],
    Pair[_e, _ec] -> -1,
    SUNT[_, _] -> 1,
    SUNT[_, _, 0, 0] -> 1/2 }


InsertFieldsHook[args__] := InsertFields[args]

CreateFeynAmpHook[args__] := CreateFeynAmp[args]


ClearSE[] := (DownValues[CalcSelfEnergy] =
  Select[DownValues[CalcSelfEnergy], #[[1, 1, 1, 0]] === Pattern &];)


OptPaint[ins_] := OptPaint[ins, $PaintSE]

OptPaint[ins_, True, ___] := Paint[ins]

OptPaint[ins_, prefix_String, suffix___String] :=
Block[ {file = ChkExist[prefix, ProcessName[ins] <> suffix <> ".ps"]},
  Paint[ins, DisplayFunction -> (Export[file, #]&)]
]


RenConst::nodef =
"Warning: `` might be renormalization constants, but have no definition."

FindRenConst[expr_] :=
Block[ {test = {expr}, orbit, isym, patt, rcs = {}, new,
SelfEnergy, DSelfEnergy},
  Needs["FeynArts`"];
  If[ $Model === "",
    InitializeModel[
      Model /. Cases[test, HoldPattern[Model -> _], Infinity] /.
        Options[InsertFields],
      GenericModel -> (GenericModel /. Options[InsertFields]),
      Reinitialize -> False ] ];

  patt = RCPattern[];

  Cases[DownValues[Dim], _[_[_[i_]], r_Integer] :> (isym[i] = x)];
  isym[other_] := other;

  orbit[other_] := other;

  While[ 
    Apply[(orbit[#1] = Range[#2])&,
      { Cases[test, SumOver[i_, r_, ___] :> {i, r}, Infinity],
        Cases[test, IndexSum[_, r___] :> r, Infinity] (*,
        Cases[DownValues[Dim], _[_[_[i_]], r_Integer] :> {i, r}]*) }, {2}];
    Length[new = Complement[
      Flatten[Cases[test, rc:patt :> Distribute[orbit/@ rc, List], Infinity]], 
      rcs, SameTest -> (isym/@ #1 === isym/@ #2 &) ]] =!= 0,
    rcs = Flatten[{new, rcs}];
    test = RenConst/@ new ];

  new = Select[ Names["Global`d*"],
    (FreeQ[rcs, #] && !FreeQ[expr, #]&)[
      ToExpression[#, InputForm, HoldPattern] ]& ];
  If[ Length[new] =!= 0, Message[RenConst::nodef, new] ];

  Cases[rcs, patt]
]


CalcRenConst[expr_, opt___Rule] :=
  (# -> WithOpt[
    RenConst[#] /. IndexSum -> psum[1] /.
      psum[n_] -> (Product[Fold[isum[], #1, {##2}], {n}]&),
    Options[kind[#]], opt])&/@
    FindRenConst[expr //. Dispatch[Subexpr[]]] /. Plus -> IntCollect


IntCollect[p__] := Plus[p] /; FreeQ[{p}, Re]

(*IntCollect[p__] := Collect[Plus[p], _Re, RCInt]*)

IntCollect[p__] :=
Block[ {RCInt},
  Replace[
    Collect[ Plus[p],
      First/@ DeleteCases[Split @ Sort @ Cases[Plus[p], _Re, Infinity], {_}],
      RCInt ],
    RCInt[x_] -> x, {1} ]
]


RCMod[rcs_, mod_] :=
Block[ {hh},
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["renormalization constants"] <>
    fincl <> "\n\n" <>
    SubroutineDecl[mod] <>
    sincl[[1]] <> "\n" <>
    VarDecl[Union[Cases[rcs, SumOver[i_, _] -> i, Infinity]], "integer"] <>
    "\n"];
  WriteExpr[hh, {sincl[[2]], rcs},
    Optimize -> True, DebugLines -> True];
  WriteString[hh, SubroutineEnd[]];
  Close[hh];
  mod
]

RCAll[mod_, mods_] :=
Block[ {hh},
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["RC invocations"] <>
    fincl[[1]] <> "\n\n" <>
    SubroutineDecl[mod] <>
    ({"\n\tcall ", $SymbolPrefix, #}&/@ mods) <>
    "\n" <> SubroutineEnd[]];
  Close[hh];
  {mod, mods}
]


Options[WriteRenConst] = {
  Folder -> "renconst",
  ExtraRules -> ExtraRules (* i.e. from WriteSquaredME *),
  FilePrefix -> FilePrefix,
  SymbolPrefix -> SymbolPrefix,
  FileHeader -> FileHeader,
  FileIncludes -> FileIncludes,
  SubroutineIncludes -> SubroutineIncludes }

WriteRenConst::norcs = "Warning: no renormalization constants found."

WriteRenConst[rcs:{___Rule}, dir_, opt___Rule] :=
Block[ {folder, xrules, prefix, $SymbolPrefix, fincl, sincl, header,
ModName, Hdr, rcmods, file, hh},

  {folder, xrules, prefix, $SymbolPrefix, header, fincl, sincl} =
    ParseOpt[WriteRenConst, opt] //. Options[WriteSquaredME];

  fincl = Flatten[{fincl, "", ""}];
  sincl = Flatten[{sincl, "", ""}];

  DefModName[MkDir[dir, folder]];

(* Part 1: CalcRenConst.F *)

  If[ Length[rcs] === 0, Message[WriteRenConst::norcs] ];

  rcmods = FileSplit[ToDoLoops[OnePassOrder[rcs /. xrules]],
    "CalcRenConst", RCMod, RCAll];

(* Part 2: renconst.h *)

  hh = OpenCode[ModName["renconst", ".h"]];

  WriteString[hh,
    Hdr["RC declarations"] <>
    VarDecl[Common["renconst"][
      MaxDims[Map[Dim, First/@ rcs, {2}]], "ComplexType" ]] <>
    "\n"];

  Close[hh];

(* Part 3: the makefile *)

  hh = OpenWrite[ModName["makefile", ""]];

  WriteString[hh, "\
OBJS :=" <> ({" \\\n  $(DIR)/", #, ".o"}&)/@ Flatten[rcmods] <> "\n\n\
$(LIB): $(LIB)($(OBJS))\n\n\
$(LIB)($(OBJS)): " <> $MakeDeps[[1]] <> "\n\n\
LIBS += $(LIB)\n\n\
VPATH := $(VPATH):$(DIR)\n\n"];

  Close[hh];

  Cases[DownValues[ModName], _[_, s_String] -> s]
]

WriteRenConst[expr_, dir_, opt___Rule] :=
  WriteRenConst[CalcRenConst[expr], dir, opt]


(* low-level output functions *)

SetLanguage::badlang = "Unknown language ``."

SetLanguage["C"] := (
  varDecl = varDeclC;
  SubroutineDecl = subroutineDeclC;
  SubroutineEnd[] = "}\n";
  callDecl = callDeclC;
  loopElem = LoopElemC;
  writeSquaredME = writeSquaredMEC;
  OpenCode = OpenC;
  CodeForm = CForm;
  $MakeDeps = {"Cdecl.d Cinline.d", "Cnum.d"};
  $CodeExt = ".c";
  $CodeIndent = "  ";
  $CodeEoln = ";";
  $CodeIf = "  if( ";
  $CodeThen = ") {\n";
  $CodeElseif = "  else if( ";
  $CodeElse = "  else {\n";
  $CodeEndIf = "  }\n";
  $CodeCall = "  ";
  $Code = "C";
)

SetLanguage["Fortran"] = (
  varDecl = varDeclF;
  SubroutineDecl = subroutineDeclF;
  SubroutineEnd[] = "\tend\n";
  callDecl = callDeclF;
  loopElem = LoopElemF;
  writeSquaredME = writeSquaredMEF;
  OpenCode = OpenFortran;
  CodeForm = FortranForm;
  $MakeDeps = {"Fdecl.d Finline.d", "Fnum.d"};
  $CodeExt = ".F";
  $CodeIndent = "";
  $CodeEoln = Sequence[];
  $CodeIf = "if( ";
  $CodeThen = ") then\n";
  $CodeElseIf = "else if( ";
  $CodeElse = "else\n";
  $CodeEndIf = "endif\n";
  $CodeCall = "call ";
  $Code = "Fortran" )

SetLanguage[lang_] := (
  Message[CodeLanguage::badlang, lang];
  OpenCode = Abort;
  $Failed )

OpenC[file_, opt___] := (
  $Code = "C";
  OpenWrite[toc <> Escape[file],
    FormatType -> CForm, opt, PageWidth -> 67] )

OpenFortran[file_, opt___] := (
  $Code = "Fortran";
  OpenWrite[tofortran <> Escape[file],
    FormatType -> FortranForm, opt, PageWidth -> 67] )

tofortran = "!" <> Escape[ToFileName[$FormCalcBin, "ToFortran"]] <> " > "

toc = "!" <> Escape[ToFileName[$FormCalcBin, "ToC"]] <> " > "

Format[addrOf[x_], CForm] := SequenceForm["&", x]

Format[addrOf[x_], FortranForm] := x


TimeStamp[] := TimeStamp@@ Date[]

TimeStamp[y_, m_, d_, H_, M_, _] :=
  ToString[d] <> "-" <>
  {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
   "Sep", "Oct", "Nov", "Dec"}[[m]] <> "-" <>
  ToString[y] <> " " <>
  ToString[H] <> ":" <> StringTake["0" <> ToString[M], -2]


(* The following routines are concerned with breaking a large
   expression into pieces the Fortran compiler will compile.
   This is controlled by two variables:

   - $BlockSize is the maximum LeafCount a single Fortran statement
     may have.  The cutting-up of expressions into such blocks is
     performed by the function WriteExpr.

   - $FileSize is the maximum LeafCount a whole file may have.
     The function which separates large expressions into file-size
     fragments is SizeSplit. *)

Attributes[batch] = {Flat}

Coalesce[(ru:Rule | RuleAdd)[v_, p_Plus], r___] :=
  Level[
    { ReplacePart[
        RuleAdd[v, Plus@@ #]&/@ Flatten[Coalesce@@ p],
        ru, {1, 0} ],
      {r} }, {2}, Coalesce ] /; LeafCount[p] > size

Coalesce[a_, b_, r___] :=
  Coalesce[batch[a, b], r] /; LeafCount[{a, b}] < size

Coalesce[a_, r___] := {batch[a], Coalesce[r]}

Coalesce[] = {}


Attributes[BlockSplit] = {Listable}

BlockSplit[expr_] := expr /; LeafCount[expr] < $BlockSize

BlockSplit[expr_] :=
Block[ {size = $BlockSize},
  List@@@ Flatten[Coalesce[expr]]
]


Attributes[NumberMod] = {Listable}

NumberMod[s__String] := StringJoin[s]


FileSplit[expr_List, mod_, writemod_, ___] :=
  {writemod[expr, mod]} /; LeafCount[expr] < $FileSize

FileSplit[expr_List, mod_, writemod_, writeall_:(#2&)] :=
Block[ {size = $FileSize},
  writeall[mod, MapIndexed[
    writemod[List@@ #1, NumberMod[mod, ToString@@ #2]]&,
    Flatten[Coalesce@@ expr] ]]
]

FileSplit[CodeExpr[vars__, expr_List], mod_,
  writemod_, writeall_:(#2&)] :=
Block[ {size = $FileSize},
  writeall[mod, MapIndexed[
    writemod[CodeExpr[vars, List@@ #1], NumberMod[mod, ToString@@ #2]]&,
    Flatten[Coalesce@@ expr] ]]
]

FileSplit[other_, r__] := FileSplit[{other}, r]


RemoveDups[expr_CodeExpr] = expr

RemoveDups[expr_] :=
  Fold[RemoveDups, #, -Range[3, Depth[#] - 1]]& @ Flatten[{expr}]

RemoveDups[expr_, lev_] :=
Block[ {obj, elem, tmps, new},
  _obj = {};
  elem[{p_, r___}] := (obj[#] = {obj[#], p})& @ expr[[p, r]];
  elem/@ Position[ expr /. {_DoLoop -> 1, _IndexIf -> 1},
    p_Plus /; LeafCount[N[p]] > minleaf,
    {lev}, Heads -> False ];
  tmps = Cases[ DownValues[obj],
    _[_[_[x_]], p_ /; Depth[p] > 3] :> {x -> dup[++dupc], Min[p]} ];
  new = Block[{Plus}, Apply[Set, tmps, {2}]; expr];
  new = #1 /. Flatten[#2] & @@
    Reap[new /. r:(_dup -> _dup) :> (Sow[r]; {})];
  Fold[ InsertDef, new, Reverse[tmps] ] //Flatten
]

InsertDef[expr_, {var_ -> tmp_, pos_}] :=
  MapAt[{tmp -> var, #}&, expr, pos]


Attributes[TmpList] = {HoldFirst}

TmpList[expr_] := Reverse[Reap[expr]]

ToTmp[expr_] := (Sow[# -> expr]; #)& @ tmp[]


Attributes[SplitExpr] = {Listable}

SplitExpr[(ru:Rule | RuleAdd)[var_, expr_]] :=
  BlockSplit @ TmpList[ ru[var, Replace[expr,
    Plus -> (If[LeafCount[#] > $BlockSize, ToTmp[#], #]&[Plus[##]]&),
    {2, Infinity}, Heads -> True]] ]

SplitExpr[other_] := other


Attributes[RhsMap] = {Listable}

RhsMap[foo_, (ru:Rule | RuleAdd)[var_, expr_]] := ru[var, foo[expr]]

RhsMap[_, other_] := other


Options[PrepareExpr] = {
  Optimize -> False,
  Expensive -> {},
  MinLeafCount -> 10,
  DebugLines -> False,
  FinalTouch -> Identity,
  ResetNumbering -> True }

PrepareExpr[expr_, opt___Rule] :=
Block[ {optim, expen, minleaf, debug, final, reset,
process, doloop, new, vars, tmps, tmp, dup, dupc = 0},
  {optim, expen, minleaf, debug, final, reset} =
    ParseOpt[PrepareExpr, opt];
  process = RhsMap[final, Flatten[SplitExpr @ Prep[#]]] &;
  If[ TrueQ[optim], process = process /. p_Prep :> RemoveDups[p] ];
  If[ reset, SymbolNumber["dup"] = SymbolNumber["tmp"] = 0 ];
  tmp[] := NewSymbol["tmp", 0];
  doloop = Hoist@@ Flatten[{expen}];
  new = Flatten[{expr}];
  vars = Cases[new, (Rule | RuleAdd)[var_, _] -> var, Infinity];
  If[ debug =!= False, new = AddDebug[new, debug] ];
  new = process[new];
  dup[c_] := dup[c] = NewSymbol["dup", 0];
  new = new /. CodeExpr[_, t_, x_] :>
    (vars = DeleteCases[vars, Alt@@ t]; x);
  tmps = Cases[new, (ru:Rule | RuleAdd)[var_, _] -> var, Infinity];
  CodeExpr[MaxDims[vars], Complement[tmps, vars], Flatten[new]]
]


NFirst[Nbb] = bb0;
NFirst[Ncc] = cc0;
NFirst[Ndd] = dd0;
NFirst[Nee] = ee0;
NFirst[Nff] = ff0;
_NFirst = 1


Attributes[AddDebug] = {Listable}

AddDebug[s_String, _] := s

AddDebug[DoLoop[expr_, ind__], tag_] :=
  DoLoop[{DebugLine[#1&@@@ {ind}, tag], AddDebug[expr, tag]}, ind]

AddDebug[i_IndexIf, tag_] := MapIf[AddDebug[#, tag]&, i]

AddDebug[(ru:Rule | RuleAdd)[var_, expr_], tag_] :=
  {ru[var, expr], DebugLine[var, tag]}

DebugLine[var_, True] := DebugLine[var /. NDim -> NFirst]


Attributes[Prep] = {Listable}

Prep[DoLoop[expr_, ind__]] := doloop[process[expr], ind]

Prep[(ru:Rule | RuleAdd)[var_, i_IndexIf]] :=
  MapIf[process[ru[var, #]]&, i]

Prep[(ru:Rule | RuleAdd)[var_, expr_List]] := Prep @
  MapIndexed[IniLHS[ru, var], ToDoLoops[expr] /. _SumOver -> 1]

Prep[(ru:Rule | RuleAdd)[var_, expr_]] := Prep @
  ru[var, SplitSums[expr]] /; !FreeQ[expr, SumOver]

Prep[(ru:Rule | RuleAdd)[var_, expr_]] := Prep @
  TmpList[ru[var, expr /. i_IndexIf :> ToTmp[i]]] /;
  !FreeQ[expr, IndexIf]

Prep[other_] := other


IniLHS[Rule, lhs_][DoLoop[rhs_, ind__], {1}] :=
  {lhs -> 0, DoLoop[RuleAdd[lhs, Plus@@ rhs], ind]}

IniLHS[Rule, lhs_][rhs_, {1}] := lhs -> rhs

IniLHS[_, lhs_][DoLoop[rhs_, ind__], _] :=
  DoLoop[RuleAdd[lhs, Plus@@ rhs], ind]

IniLHS[_, lhs_][rhs_, _] := RuleAdd[lhs, rhs]


Hoist[] = DoLoop

Hoist[a_, b__] := Hoist[a | b]

Hoist[patt_][expr_, i__] :=
Block[ {veto, abb, got = {}, c},
  veto = Alt@@ Union[
    DoIndex/@ Cases[expr, DoLoop[_, j__] :> j, Infinity],
    Cases[expr, SumOver[j_, ___] -> j, Infinity] ];
  abb = Cases[expr, x:patt /; FreeQ[x, veto], Infinity] //Union;
  abb = (got = Join[got, c = Complement[#2, got]]; #1 -> c)&@@@
    Sort[hsel[abb]/@ {i}, Length[ #1[[2]] ] < Length[ #2[[2]] ]&];
  abb[[1, 2]] = {};
  abb = (#1 -> (tmp[] -> # &)/@ #2)&@@@ abb;
  Fold[hdo, expr /. Reverse/@ Flatten[Last/@ abb], abb]
]


hsel[abb_][it:{s_, ___}] := it -> Select[abb, !FreeQ[#, s]&]

hsel[abb_][it_] := it -> Select[abb, !FreeQ[#, it]&]


hdo[DoLoop[expr_, j__], i_ -> {}] := DoLoop[expr, i, j]

hdo[expr_, i_ -> {}] := DoLoop[expr, i]

hdo[expr_, i_ -> {t__}] := DoLoop[{t, expr}, i]


SplitSums[li_List, wrap___] := SplitSums[Plus@@ li, wrap]

SplitSums[x_, wrap_:Identity] := {wrap[x]} /; FreeQ[x, SumOver]

SplitSums[x_, wrap_:Identity] :=
Block[ {term},
  term[_] = 0;
  assign[Expand[x, SumOver]];
  term[_] =.;
  #[[1, 1, 1]] wrap[Plus@@ Flatten[ #[[2]] ]]&/@ DownValues[term]
]

assign[p_Plus] := assign/@ p

assign[t_Times] := (term[#1] = {term[#1], #2})&@@ cull/@ t

assign[other_] := term[1] = {term[1], other}

cull[o_SumOver] := {o, 1}

cull[other_] := {1, other}


FindIndices[var_ -> _] := Union[Cases[var, _Symbol]]

FindIndices[t_Times] := Cases[t, SumOver[i__] -> {i}]

_FindIndices = {}

ToDoLoops[li:{__}, indices_:FindIndices] :=
Block[ {do, si},
  _do = {};
  Scan[(do[#1] = {do[#1], Tag[##]})&[indices[#], #]&, li];
  si = Flatten[Cases[DownValues[do], _[_[_[{___}]], a_] -> a]];
  DoLoop[ Last/@ #, Sequence@@ #[[1, 1]] ]&/@ 
    Split[OnePassOrder[si], #1[[1]] === #2[[1]] &]
]

(* old version:
ToDoLoops[li_List, indices_:FindIndices] :=
Block[ {do},
  do[_] = {};
  Scan[(do[#1] = {do[#1], #2})&[indices[#], #]&, li];
  Cases[ DownValues[do],
    _[_[_[{ind___}]], a_] :> DoLoop[Flatten[a], ind] ]
]
*)

ToDoLoops[other_, ___] := other


DoLoop[{a___}] := a

DoLoop[a_] := a


DoIndex[{i_, ___}] := i

DoIndex[i_] := i


Options[WriteExpr] = {
  HornerStyle -> True,
  FinalCollect -> True,
  FinalFunction -> Identity,
  Type -> False,
  TmpType -> (*Type*) "ComplexType",
  IndexType -> False,
  RealArgs -> Level[{PaVeIntegral, CutIntegral, CutMasters,
    Bput, Cput, Dput, Eput, Fput, Log, Sqrt}, {-1}],
  Newline -> "\n" }

WriteExpr[_, _[], ___] = {}

WriteExpr[hh_, CodeExpr[vars_, tmpvars_, expr_], opt___Rule] :=
Block[ {horner, fcoll, ffun, type, tmptype, indtype, rargs, newline,
var, block = 0},
  {horner, fcoll, ffun, type, tmptype, indtype, rargs, newline} =
    ParseOpt[WriteExpr, opt];
  rargs = Alt[rargs];
  horner = If[ horner =!= True, {},
    p_Plus :> (HornerForm[p] /. HornerForm -> Identity) /; Depth[p] < 6 ];
  fcoll = If[ TrueQ[fcoll], TermCollect, Identity ];
  _var = {};
  VarType[vars, type];
  VarType[tmpvars, tmptype /. Type -> type];
  VarType[Union[DoIndex/@
    Cases[expr, DoLoop[_, i__] -> i, Infinity]], indtype];
  _var =.;
  WriteString[hh,
    VarDecl[ Flatten[#2], #1[[1, 1]] ]&@@@ DownValues[var] <> "\n"];
  Flatten[{WriteBlock[hh, expr]}]
]

WriteExpr[hh_, expr_, opt___Rule] := WriteExpr[hh,
  PrepareExpr[expr, FilterOptions[PrepareExpr, opt]],
  FilterOptions[WriteExpr, opt]]


VarType[_, False] = 0

VarType[v:{__}, s:_String | {__String}] := var[s] = {var[s], v}

VarType[v_List, r:_Rule | {__Rule}] :=
  MapThread[VarSet, {v, Replace[v, r, {1}]}]

VarType[v_List, f_] := MapThread[VarSet, {v, f/@ v}]

VarSet[v_, v_] = 0

VarSet[v_, s:_String | {__String}] := var[s] = {var[s], v}


$DebugCmd = {"#ifdef DEBUG\n", "#endif\n", "DEB(", ")"}

DebugStatement[var_, tag_, {bline_:"", eline_:"", bcmd_, ecmd_}] :=
  bline <>
  "!" <> bcmd <> "\"" <> ({ToString[#], ": "}&)/@ tag <>
    var <> " =\", " <> var <> ecmd <> $CodeEoln <> "\n" <>
  eline


Attributes[WriteBlock] = {Listable}

WriteBlock[hh_, s_String] := (WriteString[hh, s <> newline]; s)

WriteBlock[hh_, DebugLine[var_, tag___]] :=
  WriteBlock[hh, DebugStatement[ToCode[var], {tag}, $DebugCmd]]

WriteBlock[hh_, DoLoop[expr_, ind__]] :=
  WriteBlock[hh, {#1, expr, #2}]&@@ DoDecl[ind]

WriteBlock[hh_, IndexIf[cond_, a_, r___]] := WriteBlock[hh, {
  $CodeIf <> ToCode[cond] <> $CodeThen,
  a, ElseIf[r],
  $CodeEndIf }]

ElseIf[] = ElseIf[RuleAdd[_, 0]] = {}

ElseIf[a_] := {$CodeElse, a}

ElseIf[cond_, a_, r___] := {
  $CodeElseIf <> ToCode[cond] <> $CodeThen,
  a, ElseIf[r] }

WriteBlock[hh_, ru_[var_, {sub__, expr_}]] :=
  WriteBlock[hh, {sub, ru[var, expr]}]

WriteBlock[_, RuleAdd[_, 0]] := Sequence[]

WriteBlock[hh_, RuleAdd[var_, expr_]] := (
  Write[hh, $CodeIndent, var -> ffun[var + FExpr[expr]], $CodeEoln];
  WriteString[hh, newline];
  var -> var + expr
)

WriteBlock[hh_, var_ -> subr_Function] := (
  Write[hh, $CodeCall,
    ffun[FExpr[subr[addrOf[var] /. NDim -> NFirst]]], $CodeEoln];
  WriteString[hh, newline];
  var -> expr
)

WriteBlock[hh_, var_ -> expr_] := (
  Write[hh, $CodeIndent, var -> ffun[FExpr[expr]], $CodeEoln];
  WriteString[hh, newline];
  var -> expr
)


FExpr[expr_] := fcoll[ expr /.
    Complex[a_, b_] -> a + cI b /.
    Dminus4 -> -2/Divergence /.
    E^x_ :> exp[x] /.
    f:rargs[__] :> RealArgs[f] /.
    p_Integer^r_Rational :> (HoldForm[#]^r &)[ N[p] ] /.
    (* p:_Integer^_Rational :> N[p] /. *)
    Den[p_, m_, d___] :> (p - m)^-d /.
    IGram[x_] :> 1/x /.
    horner ] /.
  Times -> OptTimes /.
  WeylChain -> SplitChain


Scan[(RealArgs[#[mi_, hel_, r_, n4_, ne_, pm__]] :=
  #[mi, hel, r, n4, ne, NArgs[pm]])&, CutIntegral]

RealArgs[f_] := NArgs/@ f


NArgs[i_Integer] := N[i]

NArgs[other_] := other

NArgs[i__] := Sequence@@ NArgs/@ {i}


Attributes[NZap] = {HoldAll}

NZap[sym_] := (NValues[sym] = {}; sym) /;
  Length[NValues[sym]] > 0

_NZap = {}

NClear[patt_String:"Global`*"] := Flatten[
  ToExpression[#, InputForm, NZap]&/@ Names[patt]
] //Union


Unprotect[Rule, Rational, Power]

( Format[a_ -> b_, #] := SequenceForm[a, " = ", b];
  Format[Rational[a_, b_], #] := HoldForm[a]/b;
  Format[a_^n_Integer /; n < -1, #] := (1/HoldForm[#] &)[ a^-n ];
)&/@ {FortranForm, CForm}

Format[s_Symbol^2, CForm] := HoldForm[s s];
Format[s_Symbol^-2, CForm] := 1/HoldForm[s s]

Protect[Rule, Rational, Power]


OptTimes[t__] :=
Block[ {p = Position[N[{t}], _Real, 1, Heads -> False]},
  OptNum[Times@@ Extract[{t}, p], Times@@ Delete[{t}, p]]
]

OptNum[const_, 1] = const

OptNum[const_Integer, var_] := const var

OptNum[n_?Negative r_., var_] := -OptNum[-n r, var]

OptNum[const_, var_] := HoldForm[HoldForm[const] var]

End[]


Format[ _Continuation ] = "    "
  (* eliminate those `>' in front of continuation lines so one can cut
     and paste more easily *)

FormSetup = "\
#-\n\
#:SmallSize 5000000\n\
#:LargeSize 20000000\n\
#:WorkSpace 50000000\n\
#:MaxTermSize 300000\n\
#:TermsInSmall 30000\n\
off stats;\n\
format 75;\n\n"

$BlockSize = 700

$FileSize = 30 $BlockSize

$MaxFortranName = 30

$RecursionLimit = 1024

$PaintSE = False

EndPackage[]


(* global definitions for specific models *)

powQ[1] = powQ[-1] = False

powQ[other_] := IntegerQ[other]

Sq/: (Sq[v_] = v2_) := (
  v^(n_?EvenQ) ^:= v2^(n/2);
  (v^(n_?powQ) ^:= v2^((n - Sign[n])/2) #^Sign[n])&[ v /. Pattern -> (#1&) ];
  Square[v] = v2
)

Sq/: (Sq[v_] =.) := (
  v/: v^(n_?EvenQ) =.;
  v/: v^(n_?powQ) =.;
  Square[v] =.
)


(* definitions for the Standard Model *)

Sq[EL] = 4 Pi Alfa;
Sq[Alfa] = Alfa2;
Alfa2/: Alfa2/Alfa = Alfa

Sq[GS] = 4 Pi Alfas;
Sq[Alfas] = Alfas2;
Alfas2/: Alfas2/Alfas = Alfas

Sq[SW] = SW2;
Sq[CW] = CW2;
CW2/: CW2 - 1 = -SW2;
SW2/: SW2 - 1 = -CW2;
CW2/: CW2 + SW2 = 1

Sq[MZ] = MZ2;
Sq[MW] = MW2;
Sq[MH] = MH2;

Sq[ME] = ME2;  Sq[MM] = MM2;  Sq[ML] = ML2;
Sq[MU] = MU2;  Sq[MC] = MC2;  Sq[MT] = MT2;
Sq[MD] = MD2;  Sq[MS] = MS2;  Sq[MB] = MB2

MLE[a__] := Mf[2,a];
MQU[a__] := Mf[3,a];
MQD[a__] := Mf[4,a];
Sq[Mf[a__]] = Mf2[a]

Mf[2,1] = ME;  Mf[2,2] = MM;  Mf[2,3] = ML;
Mf[3,1] = MU;  Mf[3,2] = MC;  Mf[3,3] = MT;
Mf[4,1] = MD;  Mf[4,2] = MS;  Mf[4,3] = MB

Mf2[2,1] = ME2;  Mf2[2,2] = MM2;  Mf2[2,3] = ML2;
Mf2[3,1] = MU2;  Mf2[3,2] = MC2;  Mf2[3,3] = MT2;
Mf2[4,1] = MD2;  Mf2[4,2] = MS2;  Mf2[4,3] = MB2

Conjugate[CKM[a__]] ^:= CKMC[a];
Conjugate[CKMC[a__]] ^:= CKM[a]

SUNN = 3

(* these symbols represent real quantities, i.e. Conjugate[sym] = sym
   for any of these.  Thinking e.g. of complex masses this looks
   dangerous but then again it's easy to remove any such definition.
   The function that really needs this is SquaredME. *)

Scan[ (RealQ[#] = True)&,
  { EL, Alfa, Alfa2, GS, Alfas, Alfas2,
    SW, CW, SW2, CW2,
    MW, MW2, MZ, MZ2,
    MH, MH2, MG0, MG02, MGp, MGp2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, MT, MT2,
    MD, MD2, MS, MS2, MB, MB2, _Mf, _Mf2 } ]

(* Model parameters which are defined using the parameter statement in
   Fortran (i.e. as numeric constants; see model.h) are given some
   numeric value here.  Using this information, the OptTimes function can
   significantly optimize the generated Fortran code.  The idea is to put
   everything that is known as constant at compile time in one place,
   i.e. rearrange products such that they are of the form (const)*(vars),
   then the compiler will usually collect all of these constants into one
   number. *)

Scan[ (N[#] = Random[])&,
  { cI, Alfa, Alfa2, SW2, CW, CW2,
    MW, MW2, MZ, MZ2,
    ME, ME2, MM, MM2, ML, ML2,
    MU, MU2, MC, MC2, MT, MT2,
    MD, MD2, MS, MS2, MB, MB2 } ]

SMReduce[foo_][expr_, r___] := foo[expr /.
  SW2 -> 1 - CW2 /.
  {CW -> MW/MZ, CW2 -> MW2/MZ2}, r]

SMShorten[foo_][x__] := SMReduce[foo][x] /.
  MW2 - MZ2 -> -SW2 MZ2 /.
  {MZ/MW -> 1/CW, MZ2/MW2 -> 1/CW2,
   MW/MZ -> CW, MW2/MZ2 -> CW2}

SMSimplify = SMShorten[Simplify];
SMFullSimplify = SMShorten[FullSimplify]

(* definitions for the MSSM *)

SetOptions[ CalcFeynAmp,
  NoExpand -> {USf, USfC, UASf, UASfC,
    UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC,
    UHiggs, UHiggsC, ZHiggs, ZHiggsC},
  NoBracket -> {CKM, CKMC, USf, USfC, UASf, UASfC,
    UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC,
    UHiggs, UHiggsC, ZHiggs, ZHiggsC} ]

Af[t_, g_] := Af[t, g, g]

USf[t_, g_][a_, b_] := USf[a, b, t, g];
UASf[t_][a_, b_] := UASf[a, b, t]

Conjugate[USf[a__]] ^:= USfC[a];
Conjugate[USfC[a__]] ^:= USf[a]

Conjugate[UASf[a__]] ^:= UASfC[a];
Conjugate[UASfC[a__]] ^:= UASf[a]

Conjugate[UCha[a__]] ^:= UChaC[a];
Conjugate[UChaC[a__]] ^:= UCha[a]

Conjugate[VCha[a__]] ^:= VChaC[a];
Conjugate[VChaC[a__]] ^:= VCha[a]

Conjugate[ZNeu[a__]] ^:= ZNeuC[a];
Conjugate[ZNeuC[a__]] ^:= ZNeu[a]

Conjugate[UHiggs[a__]] ^:= UHiggsC[a];
Conjugate[UHiggsC[a__]] ^:= UHiggs[a]

Conjugate[ZHiggs[a__]] ^:= ZHiggsC[a];
Conjugate[ZHiggsC[a__]] ^:= ZHiggs[a]

Conjugate[Af[a__]] ^:= AfC[a];
Conjugate[AfC[a__]] ^:= Af[a]

Conjugate[MUE] ^:= MUEC;
Conjugate[MUEC] ^:= MUE

Conjugate[Mino3] ^:= Mino3C;
Conjugate[Mino3C] ^:= Mino3

Conjugate[SqrtEGl] ^:= SqrtEGlC;
Conjugate[SqrtEGlC] ^:= SqrtEGl

SqrtEGl/: SqrtEGl SqrtEGlC = 1

Sq[SqrtEGl] = Mino3/MGl;
Sq[SqrtEGlC] = Mino3C/MGl

Sq[SA] = SA2;
Sq[CA] = CA2;
CA2/: CA2 + SA2 = 1

Sq[TB] = TB2;
Sq[SB] = SB2;
Sq[CB] = CB2;
CB2/: CB2 + SB2 = 1;
CB2/: CB2 TB2 = SB2;
CB/: CB TB = SB

Sq[CBA] = CBA2;
Sq[SBA] = SBA2;
CBA2/: CBA2 + SBA2 = 1

Sq[MGl] = MGl2;
Sq[MSf[a__]] = MSf2[a];
Sq[MASf[a__]] = MASf2[a];
Sq[MCha[a__]] = MCha2[a];
Sq[MNeu[a__]] = MNeu2[a]

Sq[MHiggs[a__]] = MHiggs2[a];
Sq[MHiggstree[a__]] = MHiggstree2[a]

Sq[Mh0] = Mh02;  Sq[Mh0tree] = Mh0tree2;
Sq[MHH] = MHH2;  Sq[MHHtree] = MHHtree2;
Sq[MA0] = MA02;  Sq[MA0tree] = MA0tree2;
Sq[MHp] = MHp2;  Sq[MHptree] = MHptree2

Scan[ (RealQ[#] = True)&,
  { TB, CB, SB, TB2, CB2, SB2, C2B, S2B,
    CA, SA, CA2, SA2, C2A, S2A,
    CAB, SAB, CBA, SBA, CBA2, SBA2,
    Mh0, Mh02, MHH, MHH2, MA0, MA02, MHp, MHp2, MGl,
    _MSf, _MSf2, _MCha, _MCha2, _MNeu, _MNeu2,
    _MHiggs, _MHiggstree } ]

MSSMReduce[foo_, red_:SMReduce][expr_, r___] :=
  red[foo][expr /. SBA2 -> 1 - CBA2, r]

MSSMShorten[foo_, red_:SMSimplify][x__] :=
  MSSMReduce[foo, red][x] /. CBA2 - 1 -> -SBA2

MSSMSimplify = MSSMShorten[Simplify];
MSSMFullSimplify = MSSMShorten[FullSimplify]


SUSYTrigExpand[expr_] := expr /. {
  SB -> sb, CB -> cb, SB2 -> sb^2, CB2 -> cb^2,
  TB -> sb/cb, TB2 -> sb^2/cb^2,
  SA -> sa, CA -> ca, SA2 -> sa^2, CA2 -> ca^2,
  C2A -> ca^2 - sa^2, S2A -> 2 ca sa,
  C2B -> cb^2 - sb^2, S2B -> 2 cb sb,
  CAB -> ca cb - sa sb, SAB -> cb sa + ca sb,
  CBA -> ca cb + sa sb, SBA -> ca sb - cb sa,
  CBA2 -> (ca cb + sa sb)^2, SBA2 -> (ca sb - cb sa)^2 }

SUSYTrigReduce[expr_] := expr /.
  {ca -> CA, sa -> SA, cb -> CB, sb -> SB}

SUSYTrigSimplify[expr_, simp_:Simplify] :=
  SUSYTrigReduce[simp[SUSYTrigExpand[expr]]]


MassDim0 = { EL, Alfa, Alfa2, GS, Alfas, Alfas2,
  SW, CW, SW2, CW2,
  TB, CB, SB, TB2, CB2, SB2, C2B, S2B,
  CA, SA, CA2, SA2, C2A, S2A,
  CAB, SAB, CBA, SBA, CBA2, SBA2,
  SqrtEGl, SqrtEGlC,
  _CKM, _CKMC,
  _USf, _USfC, _UASf, _UASfC,
  _UCha, _UChaC, _VCha, _VChaC, _ZNeu, _ZNeuC,
  _UHiggs, _UHiggsC, _ZHiggs, _ZHiggsC }

MassDim1 = { MW, MZ, MH, MG0, MGp,
  ME, MM, ML, MU, MC, MT, MD, MS, MB, _Mf,
  Mh0, MHH, MA0, MHp, MGl,
  MUE, MUEC, Mino3, Mino3C,
  _MSf, _MASf, _MCha, _MNeu, _Af, _AfC }

MassDim2 = { MW2, MZ2, MH2, MG02, MGp2,
  ME2, MM2, ML2, MU2, MC2, MT2, MD2, MS2, MB2, _Mf2,
  Mh02, MHH2, MA02, MHp2, MGl2,
  _MSf2, _MASf2, _MCha2, _MNeu2 }

MassDim[expr_] := expr /.
  (# -> Random[] &)/@ MassDim0 /.
  (# -> Mass Random[] &)/@ MassDim1 /.
  (# -> Mass^2 Random[] &)/@ MassDim2


(* make Needs["FeynArts`"] work after loading FormCalc *)

If[ !NameQ["FeynArts`$FeynArts"],
  Unprotect[$Packages];
  $Packages = DeleteCases[$Packages, "FeynArts`"];
  Protect[$Packages] ]

Null

