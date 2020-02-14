(*

This is FormCalc, Version 9.8
Copyright by Thomas Hahn 1996-2019
last modified 22 Apr 19 by Thomas Hahn

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


(* symbols from FeynArts *)

BeginPackage["FeynArts`"]

{ Topology, TopologyList, FeynAmp, FeynAmpList,
  Process, GraphID, FieldPoint, Incoming, Outgoing,
  Generic, Classes, Particles,
  Insertions, G, Mass, GaugeXi, VertexFunction,
  PropagatorDenominator, FeynAmpDenominator,
  FourMomentum, Internal, External, TheMass,
  Index, IndexDelta, IndexEps, IndexSum, SumOver,
  MatrixTrace, FermionChain, NonCommutative, LeviCivita,
  CreateTopologies, ExcludeTopologies, Tadpoles,
  InitializeModel, $Model, Model, GenericModel, Reinitialize,
  InsertFields, InsertionLevel, AmplitudeLevel,
  ExcludeParticles, ExcludeFieldPoints, LastSelections,
  Restrictions, CreateFeynAmp, Truncated,
  RenConst, MassShift, Paint, DiagramGrouping }

P$Field = (F | S | V | T | U | _Mix | _Rev)[__]

EndPackage[]


(* symbols from LoopTools *)

BeginPackage["LoopTools`"]

PaVeIntegral = A0i | B0i | C0i | D0i | E0i | F0i |
  A0 | A00 | B0 | B1 | B00 | B11 | B001 | B111 |
  DB0 | DB1 | DB00 | DB11 | C0 | D0 | E0 | F0

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

Epsi::usage =
"Epsi[id] used inside a tensor-coefficient array allows the generated 
code to address a loop integral's eps^{0,-1,-2} coefficients."


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

{Mbb, Mcc, Mdd, Mee, Mff, njcoeff}

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

{ MuTilde, dm4M, qfM, qcM, numM, intM, paveM, cutM, extM,
  dotM, abbM, fermM, sunM, helM, powM, addM, mulM, subM,
  root$$, d$$, e$$, i$$, dummy$$, g5M, g6M, g7M, dirM,
  iM, sM, EiKi, njM, tnj, xnj, bnj, b0nj, b1nj, b2nj,
  vTnj, v0nj, v1nj, v2nj, v3nj, v4nj }

FormLoopMomenta = {q1, q2}

FormInd = {Ind1, Ind2, Ind3, Ind4, Ind5, Ind6, Ind7, Ind8, Ind9}

EndPackage[]


(* symbols from the model files live in Global` *)

{ DiracMatrix, DiracSlash, ChiralityProjector,
  DiracSpinor, MajoranaSpinor, SpinorType, DiracObject,
  PolarizationVector, PolarizationTensor,
  MetricTensor, FourVector, ScalarProduct,
  Lorentz, Lorentz4, EpsilonScalar,
  SUNT, SUNF, SUNTSum, SUNEps, Colour, Gluon }

TreeCoupling::usage =
"TreeCoupling[from -> to] calculates the tree-level contribution to
the process from -> to.  TreeCoupling[..., opt] specifies options for
CreateTopologies, InsertFields, and CreateFeynAmp to be used in the
computation."

VertexFunc::usage =
"VertexFunc[from -> to] calculates the one-loop contribution to the
process from -> to.  VertexFunc[..., opt] specifies options for
CreateTopologies, InsertFields, and CreateFeynAmp to be used in the
computation."

SelfEnergy::usage =
"SelfEnergy[from -> to, mass] calculates the self-energy with incoming
particle from and outgoing particle to, taken at k^2 = mass^2. 
SelfEnergy[f] calculates the self-energy of particle f on its mass
shell.  SelfEnergy[..., opt] specifies options for CreateTopologies,
InsertFields, and CreateFeynAmp to be used in the computation."

DSelfEnergy::usage =
"DSelfEnergy[from -> to, mass] calculates the derivative with respect to
k^2 of the self-energy with incoming particle from and outgoing particle
to, taken at k^2 = mass^2.  DSelfEnergy[f] calculates the self-energy of
particle f on its mass shell.  DSelfEnergy[..., opt] specifies options
for CreateTopologies, InsertFields, and CreateFeynAmp to be used in the
computation."

K2::usage =
"K2 represents the momentum squared in SelfEnergy and DSelfEnergy."

SEHook::usage =
"SEHook[se, amp, k2 -> m2] provides a hook on the final step on the
computation of (the derivative of) a self-energy.  It should replace
k2 -> m2 in amp.  Take care that evaluating the first argument unchanged
results in a recursion."

$RCTop::usage = $RCIns::usage = $RCAmp::usage = $RCRes::usage =
"$RCTop, $RCIns, $RCAmp, and $RCRes respectively contain the output of
CreateTopologies, InsertFields, CreateFeynAmp, and CalcFeynAmp produced
by the last call to SelfEnergy, DSelfEnergy, VertexFunc, or
TreeCoupling for inspection and debugging."

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
"MassRC[f] computes the one-loop mass renormalization constant of field f. 
MassRC[f1, f2] computes the one-loop mass renormalization constant for
the f1-f2 transition. 
Field specifications of the form f @ m use mass m instead of TheMass[f]. 
For fermions the output is a list {left-handed RC, right-handed RC}. 
MassRC[..., opt] specifies options for CreateTopologies, InsertFields,
and CreateFeynAmp to be used in the computation."

FieldRC::usage =
"FieldRC[f] computes the one-loop field renormalization constant of
field f.  FieldRC[f1, f2] computes the one-loop field renormalization
constant for the f1-f2 transition.  FieldRC[f1, f2, c] subtracts c from
the self-energy entering into the calculation. 
Field specifications of the form f @ m use mass m instead of TheMass[f]. 
For fermions the output is a list {left-handed RC, right-handed RC}. 
FieldRC[..., opt] specifies options for CreateTopologies, InsertFields,
and CreateFeynAmp to be used in the computation."

TadpoleRC::usage =
"TadpoleRC[f] computes the one-loop tadpole renormalization constant of
field f.  TadpoleRC[..., opt] specifies options for CreateTopologies,
InsertFields, and CreateFeynAmp to be used in the computation."

WidthRC::usage =
"WidthRC[f] computes the one-loop width of field f. 
A field specification of the form f @ m uses mass m instead of TheMass[f]. 
WidthRC[..., opt] specifies options for CreateTopologies, InsertFields,
and CreateFeynAmp to be used in the computation."

$RenConst::usage =
"$RenConst lists the functions under which renormalization constants are
given."


BeginPackage["FormCalc`",
  {"FeynArts`", "LoopTools`", "OPP`", "Form`", "Global`"}]

(* some internal symbols must be visible for FORM/ReadForm *)

{ SUNSum, ReadForm, ReadFormDebug, FormExpr }

(* some internal symbols made visible for debugging *)

LoopIntegral = Join[PaVeIntegral, Blank/@ CutIntegral]

{ FormKins, FormProcs, KinFunc, InvSum, PairRules,
  CurrentProc, LastAmps, DenList, DenMatch,
  FormSetup, ToFPlus, UnitarityDebug, SUNObjs,
  DenyNoExp, DenyHide }

$Dminus4MaxPower::usage =
"$Dminus4MaxPower gives the maximum power of the Dminus4 symbol needed."

$Dminus4MaxPower = 1	(* for one-loop *)


(* symbols appearing in the output *)

Amp::usage =
"Amp[proc][expr1, expr2, ...] is the result of the calculation of
diagrams of the process proc.  The result is divided into parts expr1,
expr2, ..., so that index sums (marked by SumOver) apply to the whole
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

Pol::usage =
"Pol[e[i], mu, nu] is the polarization tensor of external particle i. 
Pol[ec[i], mu, nu] is its conjugate.
Pol[e1, e2, ..., mu, nu] denotes the product tensor e1*e2*... with
overall tensor indices mu and nu and Pol[..., 0, 0] its trace.
The default assumption is that the polarization tensor factorizes into
two polarization vectors, which is a specific choice for massless
particles, and can be turned off with Clear[Pol]."

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
FeynArts conventions so that the amplitude can be processed by
DeclareProcess and CalcFeynAmp."

DeclareProcess::usage =
"DeclareProcess[amps] sets up internal definitions for the computation
of the amplitudes amps."

OnShell::usage =
"OnShell is an option of DeclareProcess.  It specifies whether FORM
should put the external particles on their mass shell, i.e. apply
ki^2 = mi^2.  The special value ExceptDirac omits application of the
Dirac equation to on-shell momenta."

ExceptDirac::usage =
"ExceptDirac is a possible value for the OnShell option of DeclareProcess. 
It omits application of the Dirac equation to on-shell momenta."

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

MomElim::usage =
"MomElim is an option of DeclareProcess.  It controls in which way
momentum conservation is used to eliminate momenta.  False performs no
elimination, an integer between 1 and the number of legs substitutes
the specified momentum in favour of the others, and Automatic tries all
substitutions and chooses the one resulting in the fewest terms."

DotExpand::usage =
"DotExpand is an option of DeclareProcess.  It controls whether the
terms collected for momentum elimination are expanded again.  This
prevents kinematical invariants from appearing in the abbreviations but
typically leads to poorer simplification of the amplitude."

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
"FermionOrder is an option of CalcFeynAmp. 
For Dirac spinor chains (FermionChains -> Chiral or VA) it determines
the ordering of the external fermions within the chains.  Choices are
None, Fierz, Automatic, Colour, or an explicit ordering, e.g. {2,1,4,3}. 
None applies no reordering.  Fierz applies the Fierz identities twice,
thus simplifying the chains but keeping the original order.  Colour
applies the ordering of the external colour indices to the spinors. 
Automatic chooses a lexicographical ordering.
For Weyl spinor chains (FermionChains -> Weyl) it determines the
structuring of the amplitude with respect to the fermion chains:
Setting FermionOrder -> Mat wraps the Weyl fermion chains in Mat, too,
so that they end up at the outermost level of the amplitude and their
coefficients can be read off easily."

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

PaVeReduce::usage =
"PaVeReduce is an option of CalcFeynAmp.  False retains the one-loop
tensor-coefficient functions.  LoopTools reduces all tensors not
available in LoopTools.  True and Raw both reduce the tensors all the
way down to scalar integrals but Raw keeps the Gram determinants in
the denominator in terms of dot products while True simplifies them
using invariants."

LoopTools::usage =
"LoopTools is a value for the PaVeReduce option of CalcFeynAmp.
It specifies that only the tensor-coefficient functions not available
in LoopTools are reduced."

SortDen::usage =
"SortDen is an option of CalcFeynAmp.  It determines whether the
denominators of loop integrals shall be sorted.  This is usually done to
reduce the number of loop integrals appearing in an amplitude."

CancelQ2::usage =
"CancelQ2 is an option of CalcFeynAmp.  It controls cancellation of
terms involving the integration momentum squared.  If set to True,
powers of q1.q1 in the numerator are cancelled by a denominator, except
for OPP integrals, where the denominators are controlled by the
CombineDen option."

CombineDen::usage =
"CombineDen is an option of CalcFeynAmp.  It determines whether to
combine OPP integrals with common denominators, as in:
N2/(D0 D1) + N3/(D0 D1 D2) -> (N2 D2 + N3)/(D0 D1 D2).
True/False turns combining integrals on/off, Automatic combines integrals
only if the rank of the combined numerator is not larger than the sum
of the individual numerator ranks (which is faster)."

OPP::usage =
"OPP is an option of CalcFeynAmp.  It specifies an integer N starting
from which an N-point function is treated with OPP methods.  For
example, OPP -> 4 means that A, B, C functions are reduced with
Passarino-Veltman and D and up with OPP."

OPPMethod::usage =
"OPPMethod is an option of CalcFeynAmp.  It can take the values NumRat,
AnaRat, or Ninja, which determine how OPP integrals are treated:
NumRat assumes that the OPP library will numerically reconstruct the
rational terms and to this end provides it with the (D-4)-dimensional
part of the numerator.  AnaRat adds the rational terms analytically and
thus zeroes the (D-4)-dimensional part of the numerator.  Ninja applies
the Ninja expansion to the numerators."

NumRat::usage =
"NumRat is a possible value for the CalcFeynAmp option OPPMethod. 
It arranges for numerical evaluation of the rational terms by the OPP
library, by providing the (D-4)-dimensional part of the numerator as
the second argument to the numerator function."

AnaRat::usage =
"AnaRat is a possible value for the CalcFeynAmp option OPPMethod.
It arranges for analytical calculation of the rational terms and
consequently provides the OPP library with the (D-4)-dimensional part
of the numerator (second argument of the numerator function) set to
zero."

Ninja::usage =
"Ninja is a possible value for the CalcFeynAmp option OPPMethod.
It arranges for evaluation of the loop integrals with the Ninja library
by applying the Ninja expansion to the numerators."

OPPQSlash::usage =
"OPPQSlash is an option of CalcFeynAmp.  While the integration momentum
q1 is conceptually a D-dimensional object, it may be be treated
4-dimensionally in OPP because the numerator function actually computes
q1.q1 - MuTildeSq for every q1.q1 and the OPP libraries can reconstruct
the rational terms from the dimensionful scale MuTildeSq.
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
"NoBracket is an option of CalcFeynAmp and PolarizationSum. 
NoBracket -> {sym1, sym2, ...} specifies that sym1, sym2, ... are not
collected inside a multiplication bracket during the FORM calculation."

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

Unabbr::usage =
"Unabbr[expr] expands all abbreviations and subexpressions in expr. 
Unabbr[expr, patt] expands only those free of patt.
Unabbr[expr, !f [, patt]] does not expand inside objects that match f."

UnAbbr = Unabbr

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
"RegisterAbbr[abbr] registers a list of abbreviations so that
future invocations of CalcFeynAmp will make use of them.  Note that
abbreviations introduced for different processes are in general not
compatible."

Abbreviate::usage =
"Abbreviate[expr, f] introduces abbreviations for subexpressions
in expr.  If f is an integer, abbreviations are introduced for sums
from level f downward.  Otherwise, f is taken as a function where
f[subexpr] is invoked for the subexpressions of expr and returns
True if an abbreviation shall be introduced for it, False if not,
or the part of subexpr for which an abbreviation shall be introduced."

AbbrevSet::usage =
"AbbrevSet[expr] sets up the AbbrevDo function that introduces
abbreviations for subexpressions.  The expression given here is not
itself abbreviated but used for determining the summation indices."

AbbrevDo::usage =
"AbbrevDo[expr, lev], where lev is an integer, introduces abbreviations
for subexpressions appearing at level lev of expr or deeper.
AbbrevDo[expr, fun], where fun is a function, introduces abbreviations
for all subexpressions of expr for which fun[subexpr] returns True.
AbbrevDo must first be defined by AbbrevSet."

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

SubExpr = Subexpr

ClearSubexpr::usage =
"ClearSubexpr[] clears the internal definitions of the subexpressions
introduced by Abbreviate."

RegisterSubexpr::usage =
"RegisterSubexpr[subexpr] registers a list of subexpressions so that
future invocations of Abbreviate will make use of them."

OptimizeAbbr::usage =
"OptimizeAbbr[abbr, simp] optimizes the abbreviations in abbr by
eliminating common subexpressions.  The function simp is applied to
new abbreviations introduced for common subexpressions and defaults
to Simplify."

$OptPrefix::usage =
"$OptPrefix specifies the prefix for additional abbreviations introduced
by OptimizeAbbr, i.e. the Opt in Opt123."

SubstAbbr::usage =
"SubstAbbr[exprlist, patt] removes abbreviations matching patt and 
substitutes them back into exprlist.
SubstAbbr[exprlist, patt, deny] does not substitute variables matching deny."

SubstSimpleAbbr::usage =
"SubstSimpleAbbr[exprlist] removes (substitutes back) abbreviations which
match P$SimpleAbbr.
SubstSimpleAbbr[exprlist, deny] does not substitute variables matching deny."

P$SimpleAbbr::usage =
"P$SimpleAbbr is a pattern which matches a `simple' abbreviation of the
form var -> s where s is a number or a symbol times a number."

ExtractInt::usage =
"ExtractInt[expr] extracts the loop integrals from expr for evaluation
with LoopTools/OPP.  It output is a list {loopint, abbrexpr}, where
loopint is the list of loop integrals and abbrexpr is expr with the
loop integrals substituted by the variables on the l.h.s. of loopint."

VarSort::usage =
"VarSort[list] sorts the list of variable definitions (var -> val),
taking into account explicit priorities (var -> Prio[p][val]) and the
running numbers of variable names, e.g. abb100 is sorted after abb11."

Prio::usage =
"var -> Prio[p][val] tells VarSort that the definition var -> val is to
be sorted with priority p, i.e. after definitions with priority p - 1 and
before definitions with priority p + 1, where 100 is the default priority.
VarSort strips the Prio[p] part from the definition."

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

ToArrayRules::usage =
"ToArrayRules[expr, vars] gives the rules used for the substitution
in ToArray[expr, vars]."

Renumber::usage =
"Renumber[expr, var1, var2, ...] renumbers all var1[n], var2[n], ... in
expr."

RenumberRules::usage =
"RenumberRules[expr, vars] gives the rules used for the substitution
in Renumber[expr, vars]."

Enum::usage =
"Enum[ind] associates the indices ind with integers which are used for
determining array dimensions while the indices themselves are kept in
symbolic form.  The syntax is similar to C, e.g. Enum[a, b, c -> 5, d]
assigns a -> 1, b -> 2, c -> 5, d -> 6."

ClearEnum::usage =
"ClearEnum[] clears all Enum values."

MaxDims::usage =
"MaxDims[args] returns a list of all distinct functions in args with the
highest indices that appear, e.g. MaxDims[foo[1, 2], foo[2, 1]] returns
{foo[2, 2]}."

Keep::usage =
"Keep[expr, (name), (path/)] loads path/name.m if that file exists,
otherwise it evaluates expr and stores the result (together with the
Abbr[] and Subexpr[]) in that file.  If name is omitted, Hash[expr]
is used as filename.  path is optional and defaults to $KeepDir. 
Keep[lhs = rhs] is short for lhs = Keep[rhs, \"lhs\"]."

$KeepDir::usage =
"$KeepDir specifies the default directory for storing intermediate
expressions with Keep."

$KeepAbbr::usage =
"$KeepAbbr = False instructs Keep to not store the Subexpr[] and
Abbr[] necessary to compute abbreviated expressions.  Note that
Keep may be unable to fully recover saved expressions this way."

AbbrExpr::usage =
"AbbrExpr[expr, subexpr, abbr] contains an expression including its
subexpressions and abbreviations.  This format is used by Keep for
file storage of abbreviated expressions."


(* miscellaneous functions *)

FCPrint::usage =
"FCPrint[v, s] prints s if v <= $FCVerbose."

$FCVerbose::usage =
"$FCVerbose determines the extent of run-time messages in FormCalc. 
It ranges from 0 (no messages) to 3 (all messages)."

$FormAbbrDepth::usage =
"$FormAbbrDepth gives the minimum depth of an expression to be
considered for abbreviationing on return from FORM."

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

FormQF::usage =
"FormQF is a function applied to the loop-momentum-dependent parts
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
"TermCollect[expr] pulls common factors out of sums.
TermCollect[expr, wrap] applies wrap to terms from which a factor has
been pulled out."

Profile::usage =
"Profile[f][expr] if functionally equivalent to f[expr] but prints
information on timing and size of input and output.
Profile[f, tag][expr] prints tag instead of the function name in front
of the profiling information."

ProfIn::usage =
"ProfIn[n] contains the input of the n-th call to Profile."

ProfOut::usage =
"ProfOut[n] contains the output of the n-th call to Profile."

$ProfLine::usage =
"$ProfLine counts the calls to Profile."

SplitTerms::usage =
"SplitTerms[f, expr, n] applies f to expr, n terms at a time."

ExprHeads::usage =
"ExprHeads[expr] returns all non-system symbols and heads in expr."

ExprParts::usage =
"ExprParts[expr, hmust, (hmay)] returns all subexpressions of expr
which must contain the symbols and functions in hmust and may
contain the symbols and functions in hmay."

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
fast.  Unlike Simplify, it does not modify b and c. 
Pool[expr, wrap] applies wrap to the (b + c) part."

MapOnly::usage =
"MapOnly[f, h, patt][expr] maps f onto subexpressions of head h in expr
which must contain each of the items in patt and no symbols other than
the ones appearing in patt.
For example, MapOnly[f, h, a1|a2, b][expr] maps f onto all
h-subexpressions which contain b and a1 or a2 and no other symbols."

Creep::usage =
"Creep[f, patt][expr] applies f to subexpressions of expr which contain
only the patterns in patt."

OnSize::usage =
"OnSize[n1, f1, n2, f2, ..., fdef][expr] returns f1[expr] if
LeafCount[expr] < n1, f2[expr] if LeafCount[expr] < n2, etc., and
fdef[expr] if the expression is still larger.  fdef can take the
special value Map which means that fdef[expr] recursively applies
the entire OnSize function to the parts of expr.  If omitted, fdef
defaults to Identity."

ApplyUnitarity::usage =
"ApplyUnitarity[expr, mat, d, (wrap)] simplifies expr by exploiting the
unitarity of the d-dimensional matrix mat.  The optional argument wrap
specifies the function wrapped around linear combinations of products of
matrix elements for simplification."

ColumnIndex::usage =
"ColumnIndex is an option of ApplyUnitarity and determines whether to
simplify products of the form U[i, j] U^*[k, j]."

RowIndex::usage =
"RowIndex is an option of ApplyUnitarity and determines whether to
simplify products of the form U[j, i] U^*[j, k]."

ComplementFunction::usage =
"ComplementFunction is an option of ApplyUnitarity and specifies its
central decision function.  This function is called with two arguments,
a set of indices {a} = {a1,...,an} and its complement {Ca}.  
If it returns True, sums of the form
  Usq[a][1, i,j] := U[i,a1] U^*[j,a1] + ... + U[i,aN] U^*[j,aN] and
  Usq[a][2, i,j] := U[a1,i] U^*[a1,j] + ... + U[aN,i] U^*[aN,j]
are replaced by KroneckerDelta[i, j] - Usq[Ca][x, i,j]."

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
"MultiplyDiagrams[amp, func] multiplies the diagrams in amp with the
factor returned by the function func.  The latter is invoked for each
diagram either as func[n, id, amplitude] (for a fully inserted diagram),
or as func[n, id, generic amplitude, insertion]."

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
"IndexIf[cond, a, b] is identical to If[cond, a, b] except that a
and b are not held unevaluated.  IndexIf[cond1, a1, cond2, a2, ...]
represents a sequence of if-statements (like Which) and is equivalent
to IndexIf[cond1, a1, IndexIf[cond2, a2, ...]].  IndexIf[cond, a, b]
is converted to the Which-like syntax IndexIf[cond, a, True, b]."

MapIf::usage =
"MapIf[f, i] maps f over the expression parts only (i.e. not over the
conditions) if i is an IndexIf expression.  For all other types of
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

HoldCode::usage =
"HoldCode[expr] preserves expr until write-out to a Fortran or C file. 
Unlike Hold or HoldForm it does not inhibit the evaluation of expr."

ColourSimplify::usage =
"ColourSimplify[expr] simplifies the colour objects in expr.
ColourSimplify[tree, loop] simplifies the colour objects in
(tree^* loop)."

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
"HelicityME[tree, loop] calculates the helicity matrix elements for all
combinations of spinor chains that appear in the expression
(tree^* loop).  Terms of this kind arise in the calculation of the
squared matrix element.  The arguments do not necessarily have to be
amplitudes since they are only used to determine which spinor chains
to select from the abbreviations.  The symbol All can be used to
select all spinor chains currently defined in the abbreviations."

WeylME::usage =
"WeylME[tree, loop] arranges the Weyl chains appearing in the
expression (tree^* loop) in matrix elements suitable for evaluation
of amplitudes that were processed with FermionOrder -> Mat."

MatFactor::usage =
"MatFactor is an option of WeylME.  It specifies a function fac which
computes a correction factor fac[Fi, Fj] that is multiplied with each
matrix element Mat[Fi, Fj]."

ColourME::usage =
"ColourME[tree, loop] calculates the colour matrix elements.  ColourME
is very similar to HelicityME, except that it computes the matrix
elements for SU(N) objects, not for spinor chains."

All::usage =
"All as an argument of HelicityME, ColourME, and WeylME indicates that
all spinor chains or SUNT objects currently defined in the abbreviations
should be used instead of just those appearing in the argument."

Hel::usage =
"Hel[i] is the helicity of the ith external particle.  It can take the
values +1, 0, -1, where 0 stands for an unpolarized particle."

s::usage =
"s[i] is the ith helicity reference vector."

Mat::usage =
"Mat[Fi SUNi] (one argument) is a matrix element in an amplitude, i.e. 
an amplitude is a linear combination of Mat objects.\n\
Mat[Fi, Fj] (two arguments) appears in the squared matrix element and 
stands for the product of the two arguments, Fi Fj^*.  Such expressions 
are calculated by HelicityME, ColourME, and WeylME."

Lor::usage =
"Lor[i] is a contracted Lorentz index in a product of Dirac chains."

SquaredME::usage =
"SquaredME[tree, loop] returns the matrix element (tree^* loop). 
This performs a nontrivial task only for fermionic amplitudes: the
product of two fermionic amplitudes\n\
    M1 = a1 F1 + a2 F2 + ... and\n\
    M2 = b1 F1 + b2 F2 + ... is returned as\n\
    M1 M2^* = FF[F1][a1] FFC[F1][b1^*] Mat[F1, F1] + \n\
              FF[F2][a2] FFC[F1][b1^*] Mat[F2, F1] + ...\n\
The special case |tree|^2 can be written as SquaredME[tree] which
is of course equivalent to SquaredME[tree, tree]."

FF::usage =
"FF[mat][expr] denotes expr as the form factor of matrix element mat."

FFC::usage =
"FFC[mat][expr] denotes expr as the complex conjugated form factor of
matrix element mat."

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

SumLegs::usage =
"SumLegs is an option of PolarizationSum.  It specifies which of
the external legs to include in the polarization sum, or All for
summation over all external vector bosons."

GaugeTerms::usage =
"GaugeTerms is an option of PolarizationSum.  It controls the
treatment of the gauge-dependent terms in the polarization sum of
massless vector bosons, which should eventually cancel in
gauge-invariant subsets of diagrams. 
GaugeTerms -> True keeps the gauge-dependent terms. 
GaugeTerms -> False inserts the gauge-dependent terms, to let
potential cancellations happen, and removes remaining eta vectors. 
GaugeTerms -> Off omits the gauge-dependent terms completely."

$HelicityME::usage =
"$HelicityME contains any helicity matrix elements implicitly computed
by PolarizationSum."

$ColourME::usage =
"$ColourME contains any colour matrix elements implicitly computed
by PolarizationSum."

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

$TreeSquare::usage =
"$TreeSquare contains the default value for the TreeSquare options of
HelicityME, ColourME, WeylME, and WriteSquaredME.  Explicitly setting
the TreeSquare option overwrites $TreeSquare and thereby sets the
default value for the others."

$TreeSquare = True

LoopSquare::usage =
"LoopSquare is an option of WriteSquaredME.  It specifies whether to
add the square of the 1-loop amplitude, |M_1|^2, to the result in the
SquaredME subroutine.  This term is of order alpha^2 with respect to the
tree-level contribution, |M_0|^2.  Usually one takes into account only
the interference term, 2 Re M_0^* M_1, which is of order alpha."

$LoopSquare::usage =
"$LoopSquare contains the default value for the LoopSquare options of
HelicityME, ColourME, WeylME, and WriteSquaredME.  Explicitly setting
the LoopSquare option overwrites $LoopSquare and thereby sets the
default value for the others."

$LoopSquare = False

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
As for SubroutineIncludes, it is admissible to put a list of strings,
in which case the first two elements will be included at the beginning
of the file."

SubroutineIncludes::usage =
"SubroutineIncludes is an option of WriteSquaredME and WriteRenConst.
It specifies per-subroutine #include statements (or other declarations).
As for FileIncludes, it is admissible to put a list of strings, in which
case the first element will be included before, the second after local
variable declarations, and the third at the end of the routine."

FileHeader::usage =
"FileHeader is an option of WriteSquaredME and WriteRenConst and
specifies a file header.  This string may contain %f, %d, and %t, which
are substituted at file creation by file name, description, and time
stamp, respectively."


(* renormalization constants *)

FindRenConst::usage =
"FindRenConst[expr] returns a list of all renormalization constants
found in expr, recursively and including subexpressions.
FindRenConst[expr, h] finds only those of type h."

RenConstList::usage =
"RenConstList[h][rcs] denotes a list of renormalization constants rcs
of type h."

CalcRenConst::usage =
"CalcRenConst[expr] calculates the renormalization constants appearing 
in expr.  CalcRenConst[expr, h] calculates only those of type h."

WriteRenConst::usage =
"WriteRenConst[expr, dir] calculates the renormalization constants
appearing in expr and generates code from the results.  The resulting 
files (the code itself and the corresponding declarations) are written
to the directory dir.
WriteRenConst[expr, h, dir] writes only those of type h."

CreateTopologiesHook::usage =
"CreateTopologiesHook[l, i -> o] is the function called by SelfEnergy
and DSelfEnergy to create l-loop topologies with i incoming and
o outgoing legs.  It is normally equivalent to CreateTopologies, but
may be redefined to change the way the topologies are created."

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

RenConstHook::usage =
"RenConstHook[rc, expr] is the function called to compute expr.  It
normally returns rc -> expr but may be redefined to inspect or modify
the computation."

ClearSE::usage =
"ClearSE[] clears the internal definitions of already calculated
self-energies."

ProcName::usage =
"ProcName[amp] constructs a string suitable as symbol or filename for
the inserted topology or amplitude list amp which is unique to the model
and particle selection."

PaintSE::usage =
"PaintSE[ins, True] invokes Paint[ins]. 
PaintSE[ins, pre] invokes Paint[ins], saving the graphics in files with
names constructed from ProcName[ins], prefixed with the string pre which
may include a directory. 
PaintSE[ins, pre, suf] further appends suf to the filename. 
PaintSE[ins] executes PaintSE[ins, $PaintSE]."

$PaintSE::usage =
"$PaintSE is the default second argument of PaintSE and controls painting
of diagrams in SelfEnergy, DSelfEnergy, VertexFunc, and TreeCoupling.
Admissible values are True, False, pre, or {pre, suf}, where pre and suf
are strings used to construct the filename."

PutSE::usage =
"PutSE[{top, ins, amp, res}, pre] saves the outputs respectively of
{CreateTopologies, InsertFields, CreateFeynAmp, CalcFeynAmp} to
{file.top, file.ins, file.amp, file.res}, where file is constructed
from ProcName[ins], prefixed with the string pre which may include a
directory.
PutSE[{...}, pre, suf] further appends suf to the filename.
PutSE[{...}] executes PutSE[{...}, $PutSE]."

$PutSE::usage =
"$PutSE is the default second argument of PutSE and controls saving
of results in SelfEnergy, DSelfEnergy, VertexFunc, and TreeCoupling.
Admissible values are pre or {pre, suf}, where pre and suf are strings
used to construct the filename."

$LongitudinalSE::usage =
"$LongitudinalSE specifies that the longitudinal rather than the
transverse part of the vector-boson self-energies is taken in SelfEnergy
and DSelfEnergy."


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

ToDef::usage =
"ToDef[expr] returns the Fortran or C form of expr as a string suitable
for use in a preprocessor #define."

ToFortran::usage =
"ToFortran has been superseded by ToCode."

ToFortran = ToCode

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
subexpressions each of which has a LeafCount less than $BlockSize."

FileSplit::usage =
"FileSplit[exprlist, mod, writemod, (writeall)] splits exprlist into
batches with LeafCount less than $FileSize.  If there is only one
batch, writemod[batch, mod] is invoked to write it to file.  Otherwise,
writemod[batch, modN] is invoked on each batch, where modN is mod
suffixed by a running number, and in the end writeall[mod, res] is
called, where res is the list of writemod return values.  The optional
writeall function can be used e.g. to write out a master subroutine
which invokes the individual modules."

FunctionNames::usage =
"FunctionNames[base, ind] constructs two names out of base and ind,
where the latter are typically indices.  The first is the direct
concatenation of the elements, separated by underscores, and is for
use in file names and similar uncritical places.  The second is
limited to a maximum of $MaxFunctionName characters and should be
used for Fortran symbols to comply with compiler limits.  Truncation
is done by leaving out delimiting underscores and if that is not
enough, contracting index names down to 2 characters."

$MaxFunctionName::usage =
"$MaxFunctionName specifies the maximum length of a symbol name in
Fortran."

RuleAdd::usage =
"RuleAdd[var, expr] is written out by WriteExpr as var -> var + expr."

Call::usage =
"arr[dim] -> Call[subroutine[args]] is written out by WriteExpr as
a call to subroutine[args], which is assumed to fill the array arr."

AddrOf::usage =
"AddrOf[var] is written out by WriteExpr as the address-of operator
applied to var, e.g. &var in C."

ToVars::usage =
"ToVars[patt, symname][exprlist] introduces variables for all
subexpressions in exprlist matching patt.  The names for the
variables are determined by the function symname which receives the
expression being abbreviated and must return a symbol name for it.
If symname is a string, NewSymbol[symname, 0] is taken as naming
function, i.e. the variable names will be symname1, symname2, ...
The numbering is consecutive across ToVars calls but can be reset
by assigning SymbolNumber[symname] = 0."

$TmpPrefix::usage =
"$TmpPrefix specifies the prefix for temporary variables introduced
by PrepareExpr and WriteExpr, i.e. the tmp in tmp123."

$DupPrefix::usage =
"$DupPrefix specifies the prefix for variables introduced through
Optimize -> True in PrepareExpr and WriteExpr, i.e. the dup in dup123."

PrepareExpr::usage =
"PrepareExpr[{var1 -> expr1, var2 -> expr2, ...}] prepares a list of
variable assignments for code generation.  Expressions with a LeafCount
larger than $BlockSize are split into several pieces, as in\n
\tvar = part1\n\
\tvar = var + part2\n\
\t...\n
thereby possibly introducing temporary variables for subexpressions. 
The output is a CodeExpr[vars, tmpvars, exprlist] object, where vars
are the original and tmpvars the temporary variables introduced by
PrepareExpr."

$WriteExprDebug::usage =
"$WriteExprDebug contains the output of PrepareExpr generated during
the last invocation of WriteExpr."

WriteExpr::usage =
"WriteExpr[file, exprlist] writes a list of variable assignments to
file.  The exprlist can either be a CodeExpr object or a list of
expressions of the form {var1 -> expr1, var2 -> expr2, ...}, which is
first converted to a CodeExpr object using PrepareExpr.  WriteExpr
returns a list of the subexpressions that were actually written."

HornerStyle::usage =
"HornerStyle is an option of WriteExpr.  It specifies whether polynomial
expressions are written out in Horner form."

FinalCollect::usage =
"FinalCollect is an option of WriteExpr.  It specifies whether common
factors are collected in the final expression, just before write-out to
Fortran or C."

FinalFunction::usage =
"FinalFunction is an option of WriteExpr.  It specifies a function to be
applied to the final expressions, just before write-out to Fortran or C."

Type::usage =
"Type is an option of WriteExpr.  Admissible values are a string, rules,
a function, or False. 
If a string is given, e.g. \"RealType\", declarations of that type are
written out for all variables.
If rules are given, e.g. {x -> \"RealType\", y -> \"ComplexType\"}, they
are applied to the list of variables and any string in the result is 
taken as the type of the corresponding variable. 
If a function f is given, f[v] is formed for each variable v and if it
evaluates to a string, that string is taken as the type of variable v.  
Otherwise f[v][defs] is formed, where defs are the definitions involving
v and if the output is a string, that is taken as the type of v.
No declarations are produced for all f[v][defs] that do not evaluate to
a string, or if an explicit False is given."

TmpType::usage =
"TmpType is an option of WriteExpr.  It is the counterpart of Type for
the temporary variables.  TmpType -> Type uses the settings of the Type
option."

IndexType::usage =
"IndexType is an option of WriteExpr.  It is the counterpart of Type
for do-loop indices.  IndexType -> Type uses the settings of the Type
option."

DeclIf::usage =
"DeclIf is an option of WriteExpr.  With DeclIf -> var, where var is a
string suitable for a preprocessor variable, preprocessor statements of
the form
	#ifndef var
	#define var
	[declarations]
	#else
	[code]
	#endif
are generated to separate declarations and code, i.e. a file so generated
is supposed to be included once in the declarations section and once in
the code part.  This can be necessary in particular in Fortran, if e.g.
non-declaration statements such as statement functions need to be placed
between declarations and code."

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
"Newline is an option of WriteExpr.  It initializes the value of
$Newline, which is the string printed after each code statement."

$Newline::usage =
"$Newline is the string printed after each code statement by WriteExpr.
For finer control conditional expressions may be used which depend on
variables set by Hold[var = ...] insertions in the output stream, such
as $IndexIf and $DoLoop."

$IndexIf::usage = $DoLoop::usage =
"$IndexIf and $DoLoop are integers defined inside WriteExpr which
indicate the expression depth in IndexIf and DoLoop, respectively."

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
and/or checking statements are generated for each variable.  Admissible
values are
0 = no statements are generated,
1 = debugging statements are generated,
2 = checking statements are generated,
3 = debugging and checking statements are generated.
Debugging messages are usually generated for the expressions specified
by the user only.  To cover intermediate variables (e.g. the ones
introduced for optimization), too, specify the negative of the values
above.  The actual statements of type i (1 = debug, -2 = check-pre,
2 = check-post) are constructed from the strings in $DebugCmd[i]."

DebugLabel::usage =
"DebugLabel is an option of PrepareExpr.  It specifies with which label
debugging/checking statements are printed.  False disables printing,
True prints the variable name, and a string prefixes the variable name.
Any other value is understood as a function which is queried for each
variable assignment and its output, True, False, or a string,
individually determines generation of the debug statement."

MakeTmp::usage =
"MakeTmp is an option of PrepareExpr.  It specifies a function for
introducing user-defined temporary variables, e.g. ToVars."

Declarations::usage =
"Declarations is an option of PrepareExpr.  It specifies a pattern suitable
for Cases (i.e. patt, patt -> rhs, or patt :> rhs) that selects all objects
to be declared as variables."

FinalTouch::usage =
"FinalTouch is an option of PrepareExpr.  It specifies a function which
is applied to each final subexpression, just before write-out to file."

ResetNumbering::usage =
"ResetNumbering is an option of PrepareExpr.  It restarts numbering
of variable names for temporary and duplicate expressions at zero."

NoDebug::usage =
"NoDebug[var -> expr] does not generate a debug statement in PrepareExpr."

CodeExpr::usage =
"CodeExpr[vars, tmpvars, exprlist] is the output of PrepareExpr and
contains a list of expressions ready to be written to a code file,
where vars are the original variables and tmpvars are temporary
variables introduced in order to shrink individual expressions to a
small-enough size."

DebugLine::usage =
"DebugLine[i, var] emits a debugging statement (print-out of variable
var) of type i (1 = debug, -2 = check-pre, 2 = check-post) when written
out with WriteExpr.
DebugLine[i, var, tag] prefixes the debugging message with tag."

$DebugCmd::usage =
"$DebugCmd[i] specifies the actual debugging/checking statement of
type i (1 = debug, -2 = check-pre, 2 = check-post) in Fortran or C. 
It is given as StringForm format for three arguments: a tag prefix
(string), the variable name (string) and value (number)."

$DebugPre::usage = $DebugPost::usage =
"$DebugPre[i] and $DebugPost[i] specify commands issued before and
after the a debugging/checking statement of type i (1 = debug,
-2 = check-pre, 2 = check-post) in Fortran or C."

$DebugFF::usage =
"$DebugFF is the debug level for the form-factor modules."

$DebugAbbr::usage =
"$DebugAbbr[cat] is the debug level for the abbreviation modules,
where cat is the abbreviation category (\"s\", \"a\", or \"h\")."

$DebugNum::usage =
"$DebugNum is the debug level for the numerator modules."

$DebugRC::usage =
"$DebugRC is the debug level for the renormalization-constant modules."

SplitSums::usage =
"SplitSums[expr] splits expr into a list of expressions so that index
sums (marked by SumOver) always apply to the whole of each part. 
SplitSums[expr, wrap] applies wrap to the coefficients of the SumOver."

ToDoLoops::usage =
"ToDoLoops[list, (ifunc)] splits list into patches which must be summed
over the same set of indices.  ifunc is an optional argument, where
ifunc[expr] must return the indices occurring in expr."

DoLoop::usage =
"DoLoop[expr, ind] is a symbol introduced by ToDoLoops indicating that
expr is to be summed over the indices ind."

Dim::usage =
"Dim[i] returns the highest value the index i takes on.
A manual assignment Dim[i] = n generates correct array dimensions for
index i.  To have do-loops generated, too, assign to DoDim instead."

DoDim::usage =
"DoDim[i] returns the highest value the index i takes on, for all
indices collected during amplitude evaluation from SumOver statements.
A manual assignment DoDim[i] = n generates both correct array dimensions
and a loop over index i."

MoveDepsRight::usage =
"MoveDepsRight[r1, ..., rn] shuffles variable definitions (var -> value)
among the lists of rules ri so that the definitions in each list do
not depend on definitions in ri further to the left.  For example,
MoveDepsRight[{a -> b}, {b -> 5}] produces {{}, {b -> 5, a -> b}}, i.e.
it moves a -> b to the right list because that depends on b."

MoveDepsLeft::usage =
"MoveDepsLeft[r1, ..., rn] shuffles variable definitions (var -> value)
among the lists of rules ri so that the definitions in each list do
not depend on definitions in ri further to the right.  For example,
MoveDepsLeft[{a -> b}, {b -> 5}] produces {{b -> 5, a -> b}, {}}, i.e.
it moves b -> 5 to the left list because that depends on b."

JoinDeps::usage =
"JoinDeps is the function with which the dependent abbreviations
are joined with the given lists in MoveDepsRight and MoveDepsLeft."

OnePassOrder::usage =
"OnePassOrder[r] orders a list of interdependent rules so that the
definition of each item (item -> ...) comes before its use in the
right-hand sides of other rules."

$OnePassDebug::usage =
"When OnePassOrder detects a recursion among the definitions of a list,
it deposits the offending rules in an internal format in $OnePassDebug
as debugging hints."

FindDeps::usage =
"FindDeps[list, patt] finds all variables in the list of definitions
(var -> val) whose r.h.s. directly or indirectly depend on patt."

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
"DoDecl[v, m] returns two strings with the declarations for a loop over
v from 1 to m.  DoDecl[v, {a, b}] returns the same for a loop from a to b. 
DoDecl[v] invokes DoDim[v] to determine the upper bound on v."

CallDecl::usage =
"CallDecl[names] returns a string with the invocations of the subroutines
names, taking into account possible loops indicated by DoLoop."

$SymbolPrefix::usage =
"$SymbolPrefix is a string prepended to all externally visible symbols
in the generated code to avoid symbol collisions."


(* symbols used in the generated code *)

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
"exp[x] is the exponential function in Fortran and C."

cI::usage =
"cI represents the imaginary unit in Fortran and C."

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
"ChainXn[sL,epsL, vec..., epsR,sR] is the representation of the spinor
chain <sL|epsL|vec...|epsR|sR> in Fortran and C.  The spinors sL and
sR are multiplied by the spinor metric if epsL,R = 1, respectively. 
ChainVn starts with a sigma, ChainBn with a sigma-bar."

MomEncoding::usage =
"MomEncoding[f, i] is the encoded version of momentum i with (integer)
prefactor f."

HelDim::usage =
HelAll::usage =
"HelDim[i] and HelAll[i] are wrappers around the first index of a
HelType array.  They allow to insert additional dimensions/indices
for vectorization in Fortran."

HxH::usage =
"HxH[v1, v2] is the product of vectors v1 and v2."

IxH::usage =
"IxH[v1, v2] is the product of integer vector v1 and vector v2."

SxH::usage =
"SxH[s, v] is the product of scalar s and vector v."

SxI::usage =
"SxI[s, v] is the product of scalar s and integer vector v."

StoH::usage =
"StoH[s] vectorizes the scalar s."

ItoH::usage =
"ItoH[s] vectorizes the scalar integer s."

ConjugateH::usage =
"ConjugateH[v] conjugates the vector v."

ReH::usage =
"ReH[v] takes the real part of the vector v."

Hel0::usage = Pair0::usage = Eps0::usage = k0::usage = s0::usage =
"Hel0, Pair0, Eps0, k0, s0 are unvectorized versions respectively of
Hel, Pair, Eps, k, s."


(* system variables *)

$Editor::usage =
"$Editor specifies the command line used for editing FORM code in a
detached window."

$EditorModal::usage =
"$Editor specifies the command line used for editing FORM code in a
non-detached (modal) window."

$FormCalc::usage =
"$FormCalc gives the FormCalc version as integers {major, minor}."

$FormCalcVersionNumber::usage =
"$FormCalcVersionNumber gives the FormCalc version as a real number."

$FormCalcVersion::usage =
"$FormCalcVersion gives the FormCalc version as human-readable string."

$FormCalcDir::usage =
"$FormCalcDir is the directory from which FormCalc was loaded."

$FormCalcSrc::usage =
"$FormCalcSrc is the directory containing the FormCalc source files."

$FormCalcBin::usage =
"$FormCalcBin is the directory containing the FormCalc binary files."

$ReadForm::usage =
"$ReadForm contains the location of the ReadForm executable."

$ReadFormHandle::usage =
"$ReadFormHandle contains the ReadForm MathLink handle."

$FormCmd::usage =
"$FormCmd specifies the invocation of the FORM executable.  Arguments
are separated by a vertical bar (|)."

$DriversDir::usage =
"$DriversDir is the path where the driver programs for the generated
code are located."

$BlockSize::usage =
"$BlockSize is the maximum LeafCount a single code statement written
out by WriteExpr may have.  Any expression with LeafCount > $BlockSize
will be chopped up before being written to the code file."

$FileSize::usage =
"$FileSize gives the maximum LeafCount the expressions in a single
code file may have.  If the expressions grow larger than $FileSize,
the file is split into several pieces."


Begin["`Private`"]

$FormCalc = {9, 8}

$FormCalcVersionNumber = 9.8

$FormCalcVersion = "FormCalc 9.8 (22 Apr 2019)"

$FormCalcDir = DirectoryName[$InputFileName /.
  $HoldPattern[$InputFileName] :>
    (File /. FileInformation[System`Private`FindFile[$Input]])]

$FormCalcSrc = ToFileName[{$FormCalcDir, "FormCalc"}]

$FormCalcBin = ToFileName[{$FormCalcDir, $SystemID}]

$DriversDir = ToFileName[{$FormCalcDir, "drivers"}]

Print[""];
Print[$FormCalcVersion];
Print["by Thomas Hahn"];

$ReadForm = ToFileName[$FormCalcBin, "ReadForm"];

Check[
  $ReadFormHandle = Install[$ReadForm],
  ReadForm::notcompiled = "The ReadForm executable `` could not be \
installed.  Did you run the compile script first?";
  Message[ReadForm::notcompiled, $ReadForm];
  Abort[] ]

atexit[_[_, cmd_], ___] := $Epilog := (Uninstall[$ReadFormHandle]; cmd);
atexit[] := $Epilog := Uninstall[$ReadFormHandle];
atexit@@ OwnValues[$Epilog]

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


$FCVerbose = 1

Attributes[FCPrint] = {HoldRest}

FCPrint[v_Integer, s__] := Print[s] /; v <= $FCVerbose


$NumberMarks = False

Off[General::spell1, General::spell, Pattern::patv]


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


(* generic functions *)

FilterOpt[f_, opt___] := Sequence@@
  Cases[Flatten[{opt}], _[Alternatives@@ First/@ Options[f], _]]

ParseOpt[f_, opt___] := (
  Message[f::optx, #, f]&/@ Complement[First/@ {opt}, #];
  # //. Level[{{opt}, Options[f]}, {2}]
)&[ First/@ Options[f] ]


ToStr[x__] := StringJoin[ToString/@ Flatten[{x}]]

ToSymbol[x__] := ToExpression[ToStr[x]]

RPad[s_, n_] := s <> Table[" ", {n - StringLength[s]}]


_SymbolNumber = 0

next[stub_[h_]] :=
  ToString[h] <> ToString[stub] <> ToString[++SymbolNumber[stub[h]]]

next[stub_] := ToString[stub] <> ToString[++SymbolNumber[stub]]

NewSymbol[stub_[h_], i___] := NewSymbol[stub, i] /; Context[h] === "System`"

NewSymbol[stub_, 0] := ToExpression[next[stub]]

NewSymbol[stub_] := ToExpression[NestWhile[next, next[stub], NameQ]]


Kind[Tag[___, x_]] := Kind[x]

Kind[h_ -> _] := kind[h]

_Kind = Sequence[]

kind[h_[___], ___] := h

kind[h_, ___] := h


VarSort[li_List] := Last/@ Sort[VarPrio/@ li]

VarPrio[elem_] :=
Block[ {prio = 100, sortas},
  {prio, sortas, #}& @ VarOrd[elem]
]

VarOrd[Tag[t___, r_]] := Tag[t, VarOrd[r]]

VarOrd[lhs_ -> Prio[p_][rhs_]] := (prio = p; VarOrd[lhs -> rhs])

VarOrd[lhs_ -> rhs_] := (sortas = SymSplit[kind[lhs]]; lhs -> rhs)

VarOrd[other_] := sortas = other


Attributes[SymSplit] = {HoldRest}

SymSplit[sym_, cmd___] :=
Block[ {s = ToString[sym], i},
  i = Position[DigitQ/@ Characters[s], False][[-1,1]];
  {s, i} = StringTake[s, {{1, i}, {i + 1, -1}}];
  cmd;
  ToExpression[s][ToExpression[i]]
]


Attributes[ToArray] = {Listable}

ToArray[sym_Symbol] := SymSplit[sym]

ToArray[other_] := other

ToArray[expr_, vars__] := expr /. Dispatch[ToArrayRules[expr, vars]]

ToArrayRules[expr_, vars__] :=
  Flatten[ToArraySym[ToString/@ Flatten[{vars}]]/@ Symbols[expr]]

ToArraySym[vars_][sym_] := SymSplit[sym, If[!MemberQ[vars, s], Return[{}]]]


Renumber[expr_, vars__] := expr /. Dispatch[RenumberRules[expr, vars]]

RenumberRules[expr_, vars__] :=
Block[ {v = Flatten[{vars}], old, new},
  old = Union[Cases[expr, #[__], Infinity]]&/@ v;
  new = MapThread[Array, {v, Length/@ old}];
  Thread[Flatten[old] -> Flatten[new]]
]


ExpandSums[x_, ___] := x /; FreeQ[x, SumOver]

ExpandSums[a_. IndexDelta[i_, j_] SumOver[i_, _], h___] :=
  ExpandSums[a /. i -> j, h]

ExpandSums[a_. s__SumOver, h_:Sum] :=
  h[ExpandSums[a], Evaluate[Sequence@@ xran@@@ {s}]]

ExpandSums[other_, h___] := ExpandSums[#, h]&/@ other

xran[i_] := i

xran[i_, n_, ___] := {i, n}



Alt[l_List] := Alt@@ Flatten[l]

Alt[s_] := s

Alt[s___] := Alternatives[s]


ClearEnum[] := (Clear[enum]; enum[i_] := i)

ClearEnum[]

Enum[ind__] :=
Block[ {enumc = 0},
  enumSet/@ Flatten[{ind}]
]

enumSet[i_ -> c_] := i -> (enum[i] = enumc = c)

enumSet[i_] := i -> (enum[i] = ++enumc)


MaxDims[args__] := hRan/@
  Split[Union[Flatten[{args}]], Head[#1] === Head[#2]&]

hRan[{s__Symbol}] := s

hRan[{s__String}] := s

hRan[x_] := x[[1,0]]@@ MapThread[iRan, List@@@ x]


iRan[i__HelAll] := HelDim[Level[{i}, {2}, iRan]]

iRan[i___] := Ran@@ SortBy[Flatten[Unran/@ {i}], enum]


Ran[i:Except[1], ___, j_] := i ;; j

Ran[___, m_] := m


Unran[i_ ;; j_] := {i, j}

Unran[i_] := {1, i}



RanOff[i_ ;; _] := i

_RanOff = 1


RanDim[i_ ;; j_] := j - i + 1

RanDim[i_] := i


Attributes[Keep] = {HoldFirst}

Keep[lhs_ = rhs_] := lhs = Keep[rhs, Block[{lhs}, ToString[lhs]], $KeepDir]

Keep[lhs_ = rhs_, other__] := lhs = Keep[rhs, other]

Keep[cmd_] := Keep[cmd, ToString[Hash[Hold[cmd]]], $KeepDir]

Keep[cmd_, name_String] := Keep[cmd, name, $KeepDir]

Keep[cmd_, name_String, prefix_String] :=
Block[ {file = ChkExist[prefix, name <> ".m"], save, res, subexpr, abbr},
  If[ FileType[file] === File,
    FCPrint[1, "loading ", file];
    KeepGet[Get[file]],
  (* else *)
    save = res = cmd;
    If[ $KeepAbbr =!= False,
      subexpr = Subexpr[];
      If[ FreeQ[res, Alt[Kind/@ subexpr]], subexpr = {} ];
      abbr = Abbr[];
      If[ FreeQ[{res, subexpr}, Alt[Kind/@ abbr]], abbr = {} ];
      save = AbbrExpr[res, subexpr, abbr];
    ];
    Put[save, file];
    FCPrint[1, "saved in ", file];
    res ]
]

$KeepDir = "keep/"


AbbrExpr[expr_, {}, {}] := expr


KeepGet[AbbrExpr[expr_, subexpr_, abbr_]] := (
  RegisterAbbr[abbr];
  RegisterSubexpr[subexpr];
  expr )

KeepGet[expr_] := expr


Attributes[ToForm] = {Listable}

ToForm[x_String] := x

ToForm[x_] := ToString[x, InputForm]


ToList[p_Plus] := List@@ p

ToList[other_] := {other}


ToSeq[li:{__List}, opt___] := StringReplace[
  StringTake[ToString[li, CForm, opt], {11, -3}],
  "),List(" -> ",\n  "]

ToSeq[li_List, opt___] :=
  StringTake[ToString[li, CForm, opt], {6, -2}]

ToSeq[x_, ___] := ToForm[x]


ToBool[True] = "1"

ToBool[___] = "0"


ToCode[x_String] := x

ToCode[x_List] := StringTake[toCode[x], {6, -2}]

ToCode[x_] := toCode[x]

ToDef[x_] := StringReplace[ToCode[x], {" " -> "", "\"" -> ""}]


ToCat[n_, {}] := Table[{}, {n}]

ToCat[_, li_] := Flatten/@ Transpose[li]


Symbols[expr_] := Union[Cases[expr, _Symbol, {-1}]]


otest[f_, h_, s_][t___] := Catch[
  Scan[ If[FreeQ[{t}, #], Throw[h[t]]]&, s ];
  f[h[t]] ]

omap[f_, sym_][sym_, expr_] := f[expr]

_omap[_, expr_] := expr

MapOnly[f_, h_, s__][expr_] :=
  expr /. h -> otest[omap[f, Symbols[{s}]], h, {s}]


Creep[f_, patt_, p__] := Creep[f, patt | p]

Creep[_, patt_][expr_] := expr /; FreeQ[expr, patt]

Creep[f_, patt_][expr_] := f[expr] /; NumberQ[expr /. patt :> Random[]]

c_Creep[expr_] := c/@ expr


FromPlus[h_, p_Plus] := h@@ p

FromPlus[_, other_] := other


nobrk = Alternatives[]

DotSimplify[f_, _][expr_] := f[expr] /; FreeQ[expr, nobrk]

DotSimplify[_, f_][expr_] := f[expr]


Attributes[FPlus] = {Flat}

ToFPlus := Plus :>
  (If[ TrueQ[Re[ (4711 #1 _)[[1]] ] < 0], -FPlus@@ -{##}, FPlus[##]]&)

TermCollect[x_, wrap_:Identity] := x /. ToFPlus //.
  { FPlus[a_ b_, a_ c_] :> a FPlus[b, c],
    FPlus[a_, a_ b_] :> a FPlus[1, b] } /.
  FPlus :> (wrap[Plus[##]]&) //.
  { a_ b_ + a_ c_ :> a wrap[b + c] (*,
    a_ + a_ b_ :> a wrap[1 + b] *) }


SplitTerms[f_, p_Plus, n_Integer] := Plus@@ f/@
  Plus@@@ Partition[List@@ p, n, n, {1, 1}, {}]

SplitTerms[f_, other_, _] := f[other]


ExprHeads[expr_] := Union @ Cases[expr,
  s_Symbol /; Context[s] =!= "System`", Infinity, Heads -> True]

ExprParts[expr_, hmust_, hmay___] :=
Block[ {cond},
  defcond[Flatten[{hmust}], Flatten[{hmust, hmay}]];
  Union @ Cases[{expr}, p_Plus /; cond[ExprHeads[p]], Infinity]
]

defcond[{}, hmay_] := cond =
  Complement[#, hmay, SameTest -> MatchQ] === {} &

defcond[hmust_, hmay_] := cond =
  Complement[#, hmay, SameTest -> MatchQ] === {} && 
  Complement[hmust, #, SameTest -> (MatchQ[#2, #1]&)] === {} &


DenCancel[x_] := x //. a_ d_Den + b_ d_Den :> (a + b) d /. {
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
  Plus@@ (wrap[#2]/#1[[1,1]] &)@@@ DownValues[numer]
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

usq/: usq[a__][p_, i_, i_] + usq[i_][q_, k_, k_] :=
  usq[a, k][p, i, i] /; p =!= q

usq/: usq[a_][p_, i_, i_] + usq[a_][p_, k_, k_] :=
  usq[i, k][3 - p, a, a]


usqplus[simp_][p__] := simp[Plus[p]] /;
  Length[Position[{p}, usq[__][__], 2, 1]] > 0

usqplus[_][p__] := Plus[p]


uxdef[i___][c___] := (ux[i] = KroneckerDelta[##2] - usq[c][##] &) /;
  comp[{i}, {c}]

uxdef[i___][___] := ux[i] = usq[i]


Options[ApplyUnitarity] = {
  RowIndex -> True,
  ColumnIndex -> True,
  ComplementFunction -> (Length[#1] > Length[#2] &)
}

ApplyUnitarity[expr_, U_, dim_Integer,
  wrap:Except[_Rule]:FullSimplify, opt___Rule] :=
Block[ {rows, cols, comp, ux, uexpr},
  {rows, cols, comp} = ParseOpt[ApplyUnitarity, opt];

  Evaluate[ux@@ Range[dim]] = KroneckerDelta[##2] &;
  ux[i__] := uxdef[i]@@ Complement[Range[dim], {i}];

  UnitarityDebug = uexpr = TermCollect[expr //. Flatten[{
    If[ cols === False, {},
      (u:U[i_, j_])^n_. (uc:Conjugate[U[k_, j_]])^nc_. :>
        (usq[j][1, i, k]^# u^(n - #) uc^(nc - #) &) @ Min[n, nc] ],
    If[ rows === False, {},
      (u:U[j_, i_])^n_. (uc:Conjugate[U[j_, k_]])^nc_. :>
        (usq[j][2, i, k]^# u^(n - #) uc^(nc - #) &) @ Min[n, nc] ]
  }]] /. usq -> ux;

  uexpr /. Plus -> usqplus[wrap] /. usq -> ux /. {
    usq[a__][1, i_, k_] :> Plus@@ (U[i, #] Conjugate[U[k, #]] &)/@ {a},
    usq[a__][2, i_, k_] :> Plus@@ (U[#, i] Conjugate[U[#, k]] &)/@ {a} }
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

IndexIf[if__, else_] := IndexIf[if, True, else] /; EvenQ[Length[{if}]]

IndexIf[_, a_, True, a_] := a

IndexIf[cond_, {}, True, a_] := IndexIf[!cond, a]

IndexIf[cond1_, IndexIf[cond2_, a_]] := IndexIf[cond1 && cond2, a]

IndexIf[a___, i_IndexIf] := Level[{{a}, i}, {2}, IndexIf]


MapIf[f_, i_IndexIf] := MapAt[f, i, Table[{j}, {j, 2, Length[i], 2}]]

MapIf[f_, other_] := f/@ other


Off[Optional::opdef]

ToIndexIf[expr_, s_String] := ToIndexIf[expr, Alt[Names[s]]]

ToIndexIf[expr_, patt_:_] :=
  Fold[ singleIf, expr, Union @ Cases[expr,
    (IndexDelta | IndexDiff)[i:patt..] :> {i}, Infinity] ]

singleIf[expr_, {i__}] := expr /.
  { IndexDelta[i] -> suck[1, 1],
    IndexDiff[i] -> suck[2, 1] } /.
  a_. suck[1, x_] + a_. suck[2, y_] :> a IndexIf[Equal[i], x, y] /.
  suck[h_, x_] :> IndexIf[{Equal, Unequal}[[h]][i], x]

suck/: r_ suck[h_, x_] := suck[h, r x] /; FreeQ[r, Mat (*| SumOver*)]

suck/: suck[h_, x_] + suck[h_, y_] := suck[h, x + y]

(* suck/: suck[h_, x_] + y_ := suck[3 - h, y] /; x + y == 0 *)


(* preparations for FORM *)

ExtWF = {e, ec, z, zc}

ConjWF = ({#1 -> #2, #2 -> #1}&)@@@ Partition[ExtWF, 2] //Flatten

	(* note: KinVecs determines the order in which vectors
	   appear in FORM and hence ultimately in functions like
	   Pair and Eps *)
KinVecs = {eta, e, ec, z, zc, k}


Attributes[KinFunc] = {HoldAll}

KinFunc[args__] := Function[Evaluate[KinVecs], args]

(* (KinFunc[args__] := Function[#, args])&[ KinVecs ] *)

FromFormRules = Outer[ ToSymbol["FormCalc`", #2, #1] -> #2[#1] &,
  Range[8], KinVecs ]

FormKins = Apply[#1&, FromFormRules, {2}]

FromFormRules = Flatten[FromFormRules]


MomThread[f_][i_Index] := f[i]

MomThread[f_][p_Symbol] := f[p]

MomThread[f_][p_] := Replace[MomReduce[p], k_Symbol :> f[k], {-1}]


MomReduce[p_Plus] := Fewest[p, p + MomSum, p - MomSum]

MomReduce[p_] := p


Fewest[a_, b_, r___] := Fewest[a, r] /; Length[a] <= Length[b]

Fewest[_, b__] := Fewest[b]

Fewest[a_] := a


fvec[p_] := (vecs = {vecs, Symbols[p]}; MomReduce[p])

fvec[p_, mu_] := (vecs = {vecs, Symbols[p]}; MomThread[#[mu]&][p])


iname[type_, n_] := iname[type, n] =
  ToSymbol[StringTake[ToString[type], 3], n]


Attributes[idelta] = {Orderless}

idelta[c1:Index[Colour, _], c2_] := SUNT[c1, c2]

idelta[g1:Index[Gluon, _], g2_] := 2 SUNT[g1, g2, 0, 0]

idelta[x__] := IndexDelta[x]


ieps[c__] := Block[{eps = SUNEps[c]}, eps /; !FreeQ[eps, Index[Colour, _]]]

ieps[x__] := IndexEps[x]


isum[expr_, j__, i_] := isum[isum[expr, i], j]

isum[expr_, {i_, f_:1, t_, s_:1}] :=
  Floor[(t - f + 1)/s] expr /; FreeQ[expr, i]

isum[expr_, {i_, 1, t_}] := isum[expr, {i, t}]

isum[x__] := IndexSum@@ Flatten[{x}]


KinFunc[
  pvt[_, k, mu:Index[EpsilonScalar, _]] := z[mu];
  pvtc[_, k, mu:Index[EpsilonScalar, _]] := zc[mu];
  pvt[_, k, mu_] := e[mu];
  pvtc[_, k, mu_] := ec[mu];
  pvt[_, k, mu__] := Pol[e, mu];
  pvtc[_, k, mu__] := Pol[ec, mu]
]@@@ FormKins

pvt[fi_, -k_, mu__] := pvtc[fi, k, mu]

Conjugate[pvt] ^= pvtc


	(* this assumes factorization of the polarization tensor *)
Pol[e_, mu__] := Times@@ e/@ {mu}


Attributes[scalar] = {Orderless}

scalar[0, _] = 0

scalar[a_Symbol, b_Symbol] := a . b

scalar[a_, p:_[__]] := MomThread[scalar[a, #]&][p]


prop[0, m_Symbol, d___] := -m^(-2 d)

prop[p_, m__] := prop[-p, m] /; !FreeQ[p, -q1]

prop[p_, m_, d___] := Den[p, m^2, d]


(*
loop[a___, Den[p_, m_], b___, Den[p_, x_ m_], c___] :=
  loop[a, Den[p, m] - Den[p, x m], b, c] Den[m, 0]/(1 - x)

loop[a___, Den[p_, m1_], b___, Den[p_, m2_], c___] :=
  loop[a, Den[p, m1] - Den[p, m2], b, c]/(m1 - m2) (* *Den[m1, m2]*) /;
  m1 =!= m2
*)

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
Block[ {c = 96, newlhs, newrhs = rhs, patt},
  newlhs = lhs /.
    {Blank -> FormArg, BlankSequence | BlankNullSequence -> FormArgs} /.
    Pattern -> ((newrhs = newrhs /. #1 -> #2; #2)&);
  (newlhs /. patt[x_] :> x <> "?") -> (newrhs /. patt[x_] :> x)
]

FormArg[h_:Identity] := h[patt["sM" <> FromCharacterCode[++c]]]

FormArgs[h_:Identity] := h["?" <> FromCharacterCode[++c]]


OrdSq[r:_[_, rhs_]] := {{}, r} /; VectorQ[lhs, FreeQ[rhs, #]&]

OrdSq[r_] := {r, {}}

SortSq[dv_] :=
Block[ {lhs = #[[1,1,1]]&/@ dv},
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

CurrentProc = CurrentOptions = Sequence[]

Options[DeclareProcess] = {
  OnShell -> True,
  Invariants -> True,
  Transverse -> True,
  Normalized -> True,
  InvSimplify -> True,
  MomElim -> Automatic,
  DotExpand -> False,
  Antisymmetrize -> True }

Attributes[DeclareProcess] = {Listable}

DeclareProcess::incomp =
"Calculation of incompatible process(es) attempted. 
Need ClearProcess[] between processes with different external kinematics."

DeclareProcess::syntax =
"Wrong syntax: DeclareProcess expects FeynAmpList objects as arguments."

DeclareProcess[fal:FeynAmpList[__][___].., opt___Rule] := (
  DeclP[#, ParseOpt[DeclareProcess, opt], CurrentOptions];
  Level[{fal}, {2}, FormAmp[#]]
)& @ ChkProc[{LastAmps = fal, CurrentProc}, DeclareProcess, Abort[]]

DeclareProcess[l__List, opt___Rule] := DeclareProcess@@ Flatten[{l, opt}]

_DeclareProcess := (Message[DeclareProcess::syntax]; Abort[])


DeclP[_, opt_, opt_] := Null

DeclP[proc_, opt:{onshell_, inv_, transv_, norm_,
                  invsimp_, momelim_, dotexp_, antisymm_}, ___] :=
Block[ {signs, masses, kins, dot, n, neglect, square, id, momrange,
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
  FormVectors = Transpose[KinFunc[Evaluate[KinVecs]]@@@ kins];

  Attributes[dot] = {Orderless, Listable};

  kikj = dot[Moms, Moms];
  If[ onshell =!= False && Length[kins] === Length[masses],
    kikj = MapThread[Set, {kikj, masses}] ];

  If[ transv, eiki = KinFunc[{EiKi[e, k], EiKi[ec, k]}]@@@ kins ];

  If[ norm, eiei = KinFunc[e.ec -> -1]@@@ kins ];

  Switch[ Length[masses],
    0 | 1,
      Null,
    2,
      (*dot[ Moms[[2]], Moms[[2]] ] =.;*)
      MomSubst = {s_Spinor :> s, Moms[[2]] -> Moms[[1]]};
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

        If[ invsimp && onshell =!= False && Legs =!= 3,
          invproc = ToForm[MapIndexed[
            { "id `foo'(sM?) = `foo'(sM, sM*replace_(",
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

  kikj = (#1[[1]] -> #2)&@@@ (DownValues[dot] /. dot -> Dot);
  FormSymbols = Flatten[{FormSymbols, Last/@ kikj}];

  momrange = Range[Legs];
  If[ momelim > 0 && momelim <= Legs,
    momrange = Append[Delete[momrange, {momelim}], momelim] ];

  FormProcs = "\n\
#define Legs \"" <> ToString[Legs] <> "\"\n\
#define Invariants \"" <> ToSeq[Cases[invs, _Symbol]] <> "\"\n\
#define OnShell \"" <> ToString[onshell] <> "\"\n\
#define MomElim \"" <> ToString[momelim && Length[MomSubst] === 0] <> "\"\n\
#define MomRange \"" <> ToSeq[momrange] <> "\"\n\
#define DotExpand \"" <> ToBool[dotexp] <> "\"\n\
#define Antisymmetrize \"" <> ToBool[antisymm] <> "\"\n\
#define Dminus4MaxPower \"" <> ToString[$Dminus4MaxPower] <> "\"\n\n" <>
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
    Level[ {{{{1}}}, Select[FromFormRules, !FreeQ[ {amps}, #[[1]] ]&]},
      {4}, Max];
  MapIndexed[ToAmp, FeynAmpList[Process -> proc][amps] /.
    Reverse/@ FromFormRules /.
    PolarizationVector | PolarizationTensor -> pvt /.
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
  AmpName[Select[id, FreeQ[#, Number]&]] ->
    {FLines[TrivialSums[ampden[amp /. toGen]]], {}}

LevelSelect[Automatic][r__, gm_ -> ins_] := LevelSelect[
  Which[
    !FreeQ[ins, Particles], Particles,
    !FreeQ[ins, Classes], Classes,
    True, Generic ]
][r, gm -> ins]

LevelSelect[Generic][id_, _, gen_, ___] :=
  AmpName[Select[id, FreeQ[#, Classes | Particles | Number]&]] ->
    {FLines[ampden[gen /. toGen]], {}}

LevelSelect[lev_][id_, _, gen_, gm_ -> levins_] :=
Block[ {ins, red, amp, pc},
  _pc = 0;
  ins = TrivialSums/@ Cases[{levins}, Insertions[lev][r__] :> r, Infinity];
  ins = ampden[GetDen[gen], gm]/@ ins;
  red = Thread[Flatten[ReduceIns[gm, Transpose[ins]]], Rule];
  amp = gen /. red[[1]] /. {d_Den :> (d /. small[m_] :> m), _small -> 0};
  If[ Length[ red[[2]] ] === 0, amp *= Length[ins] ];
  AmpName[id] -> {FLines[amp], red[[2]]}
]


FormIns[name_ -> {}] := name

FormIns[name_ -> ins_] := "i" <> name -> name *
  Plus@@ (Level[#, {2}, "replace_"]&)/@ Transpose[Thread/@ ins]


AmpDen[amp_] := amp GetDen[amp]

AmpDen[1, _] = Identity

AmpDen[den_, gm_] := Block[ {ins = #1},
  ins[[-1]] *= den /. Thread[gm -> #1];
  ins ]&


IdDen[amp_] := amp

IdDen[_, _] = Identity


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


Attributes[sumover] = {HoldAll}

sumover[sum_, _][i:Index[Colour | Gluon, _], r__] := sum[i, r]

sumover[_, sum_][i__] := sum[i]


SUNSum[i_, _, External] := SUNSum[i]


TrivialSums[ins_ -> _] := TrivialSums[ins]

TrivialSums[ins_] := ins /; FreeQ[ins, SumOver]

TrivialSums[ins_] := ins /. SumOver -> CarryOut[ins /. _SumOver -> 1]


CarryOut[test_][i_, n_, r___] :=
  sumover[SUNSum, (ranges = {i -> i == n, ranges}; SumOver)][i, n, r] /;
  !FreeQ[test, i]

_CarryOut[_, v_, External] := Sqrt[v]

_CarryOut[_, v_, ___] := v


Attributes[OffShell] = {Listable}

OffShell[fal:FeynAmpList[__][___].., extmass___Rule] :=
Block[ {off},
  offdef@@@ {extmass};
  Sequence@@ ampoff/@ {fal}
]

offdef[n_, f_Function] := off[n] = f

offdef[n_, m_] := off[n] = m &

ampoff[amp_] :=
Block[ {c = 0},
  (#1 /. Flatten[#2])&@@ Reap[amp /. (Process -> proc_) :>
    (Process -> Map[partoff[#, off[++c]]&, proc, {2}])]
]

partoff[x_, _off] := x

partoff[{fi_, p_, m_, qn__}, f_] := (
  Sow[(s:DiracSpinor | MajoranaSpinor | _SpinorType)[k:_. p, _] :> s[k, #]];
  {fi, p, #, qn}
)& @ f[m]


Combine[fal:FeynAmpList[__][___]..] := (
  ChkProc[#, Combine];
  Level[#, {1}, #[[1,0]]] )& @ {fal}

Combine[amp:Amp[_][___]..] :=
Block[ {comp},
  ChkProc[#, Combine];
  _comp = 0;
  Map[Component, #, {2}];
  _comp =.;
  #[[1,0]]@@ Cases[DownValues[comp], _[_[_[x___]], r_] :> r x]
]& @ {amp}

Component[r_ s__SumOver] := comp[s] += r

Component[other_] := comp[] += other


MultiplyDiagrams[diags_, f_] :=
Block[ {DiagNo = 0},
  MulDiags[f][diags]
]

MulDiags[f_][fal:FeynAmpList[h__][___]] :=
  MulDiags[f, Count[AmplitudeLevel /. {h}, Generic]]/@ fal

MulDiags[f_][fal:FormAmp[_][___]] := MulDiags[f, 1]/@ fal

MulDiags[f_, gc_:0][FeynAmp[id_, q_, gen_, gm_ -> ins_]] := (
  DiagNo += gc;
  FeynAmp[id, q, gen, gm -> MulIns[f, id, gen, gm]/@ ins] )

MulDiags[f_, ___][FeynAmp[id_, q_, amp_]] :=
  FeynAmp[id, q, f[++DiagNo, id, amp] amp]

i_MulIns[ins_ -> more_] := i[ins] -> i/@ more

MulIns[f_, id_, gen_, gm_][{r__, fac_}] :=
  {r, f[++DiagNo, id, gen, Thread[gm -> {r, fac}]] fac}


TagDiagrams[diags_, tag_:Diagram] := MultiplyDiagrams[diags, tag[#1]&]


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

DenyNoExp = {_String, Spinor, Den, intM, PowerOf,
  (*IndexDelta,*) IndexSum, SumOver, Sequence@@ SUNObjs}

DenyHide = Level[{SumOver, PowerOf, IndexDelta, IndexEps, SUNObjs},
  {-1}, Alternatives]

FinalFormRules = {Power -> FormPower, Complex[re_, im_] :> re + "i_" im}


Attributes[CalcFeynAmp] = {Listable}

Options[CalcFeynAmp] = {
  CalcLevel -> Automatic,
  Dimension -> D,
  NoCostly -> False,
  FermionChains -> Weyl,
  FermionOrder -> Automatic,
  Evanescent -> False,
  InsertionPolicy -> Default,
  SortDen -> True,
  CombineDen -> Automatic,
  PaVeReduce -> False,
  CancelQ2 -> True,
  OPP -> 6,
  OPPMethod -> Ninja,
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

CalcFeynAmp::ddim = "Warning: `` implies the use of Fierz identities \
which are in general not applicable in D dimensions.  \
You know what you are doing."

CalcFeynAmp[FormAmp[proc_][amp___], opt___Rule] :=
Block[ {lev, dim, nocost, fchain, forder, evanes,
inspol, sortden, combden, pavered, cancelq2, opp, oppmeth, oppqsl,
g5test, g5eps, noexp, nobrk, momrul, pre, post, tag, edit, retain,
uniq, vecs, kc = 0, hd, ic, inssym, mmains,
indices = {}, ranges = {}, haveferm = False, indsym, momrange,
intmax, extmax = 0, ampden, vars, hh, amps, res, traces = 0},

  { lev, dim, nocost, fchain, forder, evanes, inspol,
    sortden, combden, pavered, cancelq2,
    opp, oppmeth, oppqsl,
    g5test, g5eps, noexp, nobrk, momrul,
    pre, post, tag, edit, retain } = ParseOpt[CalcFeynAmp, opt];

  opp = opp /. {True -> 1, Rational -> -1, False -> 100};

  If[ dim === 0,
    If[ Length[#] > 0, Message[CalcFeynAmp::ddim, #] ]&[{
      If[ forder =!= None, FermionOrder -> forder, {} ],
      If[ fchain === Weyl, FermionChains -> fchain, {} ]
    } //Flatten] ];

  NoExpandRule = {
    p_Plus :> (TermCollect[p] /. Plus -> addM) /; !FreeQ[p, #1] && FreeQ[p, #2],
    p_Plus^n_?Negative :> (TermCollect[p] /. Plus -> addM)^n /; FreeQ[p, #2]
  }&[ Alt[noexp],
      Alt[{FormLoopMomenta, FormVectors, DenyNoExp}] ];

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
    { g:G[_][_][__][_] :> g,
      IndexDelta -> idelta,
      IndexEps -> ieps,
      IndexSum -> isum } /.
    FormMom[proc] /.
    FourMomentum[Internal, i_] :> FormLoopMomenta[[i]] /.
    PolarizationVector | PolarizationTensor -> pvt /.
    Prepend[momrul,
      (DiracSpinor | MajoranaSpinor | _SpinorType)[s_. p_Symbol, m_] :>
        Spinor[p, Neglect[m], s]] /.
    { LeviCivita -> Eps,
      ScalarProduct -> scalar,
      PropagatorDenominator -> prop,
      FeynAmpDenominator -> loop,
      FourVector -> fvec,
      g:G[_][_][__][_] :> g,	(* for gn *)
      ChiralityProjector[c_] :> ga[(13 - c)/2],
      (DiracSlash | DiracMatrix) -> MomThread[ga],
      NonCommutative -> noncomm,
      MetricTensor[mu_, _]^2 :> "d_"[mu, mu],
      MetricTensor -> "d_" };

  ampden = If[ combden === False, IdDen,
    DenList = Reverse[SortBy[DenList, Length]];
    AmpDen ];
  amps = LevelSelect[lev]@@@ amps;
  Apply[ (amps[[##,0]] = DenExtend)&,
    Sort[Thread[-Length@@@ Extract[amps, #] -> #]& @
      Position[amps, _ExtDen]], {2} ];

  intmax = Max[0, Cases[amps, i_intM :> Length[i], Infinity]];

  vecs = Complement[Flatten[vecs], Flatten[{FormLoopMomenta, FormVectors}]];
  FormVectors = DeleteCases[Append[FormVectors, vecs], {}];

  indsym[type_, n_] := (
    indices = {indices, #};
    ranges = {IndexDim[type, #], ranges};
    indsym[type, n] = #
  )& @ iname[type, n];

  hd = Amp[Apply[{#1, k[++kc], ##3}&, proc /. Index -> indsym, {2}]];

  amps = DeleteCases[amps, {_ -> {0, _}}];
  If[ Length[amps] === 0, Return[hd[0]] ];

  mmains = Cases[DownValues[inssym], _[_[_[ins_]], s_Symbol] :> s == ins];

  amps = ToCat[2, Thread/@ amps] /. NoExpandRule;
  If[ Position[amps[[2]], _Rule, {3}, 1] === {}, inspol = Begin ];

  amps = amps /. Index -> indsym /. FinalFormRules;
  mmains = mmains /. Index -> indsym;
  ranges = Append[Union[Flatten[ranges]] /. Index -> indsym, i_ :> i == 0];

  indices = Union[Flatten[ {indices, iM,
    Cases[DownValues[indsym], _[_, s_Symbol] :> s],
    Cases[amps, SumOver[i_, ___] | IndexSum[_, i_, ___] |
                IndexDelta[i__] | IndexEps[i__] :> i, Infinity],
	(* possibly some colour or gluon indices have been fixed by
	   the user in FeynArts; these would by default be declared
	   as symbols and not be seen by FORM: *)
    Cases[Cases[amps, SUNObjs[__], Infinity] /. SUNSum -> (#1 &),
      _Symbol, {-1}]} ]];

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
#define Dim \"" <> ToString[dim] <> "\"\n\
#define InsertionPolicy \"" <> ToForm[inspol] <> "\"\n\
#define IntMax \"" <> ToForm[intmax] <> "\"\n\
#define ExtMax \"" <> ToForm[Max[intmax, extmax]] <> "\"\n\
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
#define CancelQ2 \"" <> ToBool[cancelq2] <> "\"\n\
#define OPP \"" <> ToForm[Max[2, opp]] <> "\"\n\
#define OPPMethod \"" <> ToForm[oppmeth] <> "\"\n\n" <>
    vars[[1]] <> "\n\
table HEL(1:" <> ToString[Legs] <> ");\n" <>
    Array[{"fill HEL(", ToString[#], ") = ", ToForm[Hel[#]], ";\n"}&,
      Legs] <> "\n\n"];

  FormWrite[hh, amps[[1]] /. MomSubst];

  WriteString[hh,
    FormProcs <>
    FormConst[vars, nobrk] <>
    "#procedure Insertions\n"];
  res = Plus@@ FormWrite[ hh, FormIns/@ amps[[2]] ];
  WriteString[hh, ".sort\ndrop;\n\n"];
  FormWrite[hh, FormCalc`Result -> res];
  WriteString[hh, "#endprocedure\n\n" <>
    FormCode["Common.frm"] <>
    FormCode["CalcFeynAmp.frm"]];
  Close[hh];

  nobrk = Alt[nobrk];	(* for DotSimplify *)
  FormPre[amps];

  hd@@ post/@ FormOutput[mmains][edit, retain][[1]]
]

CalcFeynAmp[amps__, opt___Rule] := CalcFeynAmp[
  DeclareProcess[amps, FilterOpt[DeclareProcess, opt]],
  FilterOpt[CalcFeynAmp, opt] ]


(* FORM interfacing *)

trace[g__] := (traces = Max[traces, ++fline];
  NonCommutativeMultiply[ga[], g] /. {
    a_. ga[li__] ** ga[6] + a_. ga[li__] ** ga[7] :> a ga[li],
    a_. ga[6] + a_. ga[7] :> a
  } /. ga -> om)

chain[g1_, g___] := (fline = Max[fline + 1, 100];
  haveferm = True;
  NonCommutativeMultiply[g1, ga[], g] /. ga -> om)

dobj[g1_, g___][di__] :=
Block[ {fline = iM},
  dirM[NonCommutativeMultiply[g1, ga[], g] /. ga -> om, di]
]


om[1] := "g_"[fline]

om[5] := "g5_"[fline]

om[6] := "g6_"[fline]/2

om[7] := "g7_"[fline]/2

om[q1] := "g_"[fline, q1] + I "g5_"[fline] MuTilde /; oppqsl

om[li___] := "g_"[fline, li]


FLines[expr_] :=
Block[ {fline = 0},
  expr /.
    MatrixTrace -> trace /.
    FermionChain -> chain /.
    DiracObject -> dobj /.
    NonCommutativeMultiply[a_] :> a
]


Attributes[FormWrite] = {Listable}

FormWrite[hh_, lhs_ -> rhs_] :=
Block[ {fline = 0},
  Write[hh, "L ", lhs, " = ", rhs, ";"];
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

FormId[lhs_ -> rhs_] :=
  {"id ", ToForm[lhs], " = ", ToForm[rhs], ";\n"}

FormId[fun_] :=
  {"#call ", ToForm[fun], "\n"}

FormSq[_[lhs_, h_[x___]]] :=
  {"id ", #1, "^2 = ", #2, ";\n",
   "id ", #1, "^-2 = 1/(", #2, ");\n"}&[ ToForm[lhs], ToForm[h[x]] ] /;
  Context[h] === "System`"

FormSq[_[lhs_, rhs_]] :=
  {"id ", #1, "^2 = ", #2, ";\n",
   "id ", #1, "^-2 = ", #2, "^-1;\n"}&[ ToForm[lhs], ToForm[rhs] ]


NCFuncs = Spinor | g5M | g6M | g7M

FormVars[dim_, expr_, inds_:{}, ranges_:{}] :=
Block[ {theexpr, vars, func},
  theexpr = {expr, dim, FormSymbols};

  vars = Complement[Symbols[theexpr],
    Flatten[{inds, FormLoopMomenta, FormVectors}]];

  func = Complement[
    Cases[Head/@ Level[theexpr, {1, -2}, Heads -> True], _Symbol],
    Flatten[{FormLoopMomenta, FormVectors, ExtWF, Pol,
      ga, MatrixTrace, FermionChain, NonCommutativeMultiply,
      Rule, Equal, Plus, Times, Power, Dot, Rational}] ];

  { { FormDecl["s ", vars],
      FormDecl["cf ", DeleteCases[func, NCFuncs]],
      FormDecl["f ", Cases[func, NCFuncs]],
      If[ dim === 4, "", "d D;\n"],
      FormDecl["i ", Replace[inds, ranges, 1]], "\
t Pol;\n\n\
#define LoopMomenta \"", ToSeq[FormLoopMomenta], "\"\n\
#define Vectors \"`LoopMomenta',\n\
  ", ToSeq[FormVectors], "\"\n\
v `Vectors';\n" },
    inds, vars, func }
]


FormConst[vars_, nobrk_] := {
  "#procedure ConstBracket\n",
  FormDecl["ab ", DeleteCases[
    Complement[Flatten @ vars[[{3, 4}]], Flatten[{D, Dminus4, nobrk}]],
      (x_ /; MemberQ[{"FeynArts`", "FormCalc`", "LoopTools`", "Form`"},
        Context[x]]) | SUNObjs | Conjugate ]],
  "#endprocedure\n\n" }


tempnum = 1000 $KernelID + 1 /. HoldPattern[$KernelID] -> 0

FormFile[stub_, n_] :=
  ToFileName[Directory[], stub <> ToString[n] <> ".frm"]

OpenForm[stub_:"fc"] :=
Block[ {hh},
  While[ FileType[tempfile = FormFile[stub, tempnum]] =!= None,
    ++tempnum ];
  FCPrint[1, "\npreparing FORM code in " <> tempfile];
  hh = OpenWrite[toform <> Escape[tempfile],
    FormatType -> InputForm, PageWidth -> 73];
  WriteString[hh, FormSetup];
  hh
]

toform = "!" <> Escape[ToFileName[$FormCalcBin, "ToForm"]] <> " > "


Attributes[FormExpr] = {HoldAll}

FormExpr[x__] := Block[{addM = Plus, mulM = Times, powM = Power, subM}, {x}]


	(* simplification settings *)

$FormAbbrDepth = 3

FormPre = AbbrevSet[#, Preprocess -> FormMat]&

FormSub = AbbrevDo

FormDot = DotSimplify[Simplify, TermCollect]

FormMat = TermCollect

FormQC = FormNum =
  TermCollect[DenCancel[#]] (* //. a_ p_Pair + b_ p_Pair :> (a + b) p *) &

FormQF = Simplify


Attributes[Profile] = {HoldFirst}

Profile[f_Symbol] := prof[f, Block[{f}, ToString[f]]]

Profile[other__] := prof[other]

$ProfLine = 0

prof[f_, tag___][expr_, r___] :=
Block[ {in, n = ++$ProfLine},
  ProfIn[n] = expr;
  in = tag <> RPad[ToStr["  ProfIn[", n, "]:", LeafCount[expr]], 23];
  If[ $Notebooks === False, WriteString["stdout", in, "\r"] ];
  ( WriteString["stdout", in,
      RPad[ToStr["  ProfOut[", n, "]:", LeafCount[#2]], 24],
      "  Time:", InputForm[#1], "\n"];
    ProfOut[n] = #2
  )&@@ Timing[f[expr, r]]
]


Attributes[FormExec] = {HoldAll}

FormExec[s__Set] := Block[{s}, ReadForm[$FormCmd, tempfile]]


FormOutput[r___][edit_, retain_] :=
Block[ {res},
  Switch[ edit,
    True, Pause[1]; Run[StringForm[$Editor, tempfile]]; Pause[3],
    Modal, Pause[1]; Run[StringForm[$EditorModal, tempfile]] ];
  FCPrint[1, "running FORM... "];
  If[ Check[
        res = Set@@@ Level[
          {r, FromFormRules, {Dot -> p$$}},
          {2}, FormExec ];
        True,
        False] &&
    !retain, DeleteFile[tempfile] ];
  FCPrint[1, "ok\n"];
  res
]


FormPower[Rational[x_, y_], Rational[n_, d_]] := ("root_"[d, x y^(d - 1)]/y)^n

FormPower[y_?NumberQ, Rational[n_?Negative, d_]] := ("root_"[d, y^(d - 1)]/y)^-n

FormPower[x_?NumberQ, Rational[n_, d_]] := "root_"[d, x^n]

FormPower[f:_[__], n_?Negative] := powM[f, n]

FormPower[x_, n_Integer] := x^n

FormPower[other__] := powM[other]


(* things to do when the amplitude comes back from FORM *)

root$$[d_, x_] := x^(1/d)

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

fermM[0] = 0

fermM[x_] := fermM[x] = NewSymbol["F"]


pair[x_] := pair[x] = NewSymbol["Pair"]

eps[x_] := eps[x] = NewSymbol["Eps"]

pol[x_] := pol[x] = NewSymbol["Pol"]


abb[x_] := x /; Depth[x] < $FormAbbrDepth

abb[n_?NumberQ x_] := n abb[x]

(*abb[x_?AtomQ] = x*)

abb[x_^n_Integer] := abb[x]^n

abb[x_] := abb[x] = NewSymbol["Abb"]


abbM[x_] := x /; FreeQ[x, Pair | Eps]

abbM[x_] := abb[x /. {p_Pair :> pair[p], e_Eps :> eps[e], t_Pol :> pol[t]}]

dotM[x_] := abbM[FormDot[x]]


sunM[x_] := sunM[x] = NewSymbol["SUN"]


numM[n_, x_, i___] := tonum[n, qcount@@ #1, i][##]&@@
  CoefficientList[DenCancel[x], {njM, dm4M}]

tonum[_, r_, i___][{n4_, neps_:0, ___}] :=
  { r,
    num[Num[FormNum[n4]], i],
    num[Num[FormNum[neps]], i] }

tonum[n_, r_, i___][{n4_, ___}, {mu_, ___}, {t3_, ___}, {t2_, ___}] :=
Block[ {mindeg = r - n + 3, ct, cm, numstub},
  { r,
    numstub = "Num";
      num[Num[FormNum[n4]], i],
    numstub = "MuExp";
      num[Num[{0} -> FormNum[Coefficient[mu, tnj, 0]]], i],
    numstub = "T3Exp";
      num[Num@@ Flatten @
        Table[ct = Coefficient[t3, tnj, r - d];
          Table[{d, m} -> FormNum[Coefficient[ct, MuTildeSq, m]],
            {m, 0, d/2}],
          {d, 0, mindeg}], i],
    numstub = "T2Exp";
      num[Num@@ Flatten @
        Table[ct = Coefficient[t2, tnj, r - d];
          Table[cm = Coefficient[ct, xnj, x];
            Table[{d, x, m} -> FormNum[Coefficient[cm, MuTildeSq, m]],
              {m, 0, (d - x)/2}],
            {x, 0, d}],
          {d, 0, mindeg}], i] }
]


numstub = "Num"

num[n_[], ___] := None	(* keep 'n_' here for dv *)

num[n_] := num[n] = NewSymbol[numstub]

num[n_, i_List] := Level[{{n}, Select[i, !FreeQ[n, #]&]}, {2}, num]

num[n_, i__] := num[n, i] = NewSymbol[numstub][i]


qcount[p_Plus, ___] := Max[qcount/@ List@@ p]

qcount[t_Times, ___] := Plus@@ qcount/@ List@@ t

qcount[x_, ___] := Count[x, q1, {-1}]


qfM[x_] := Collect[x, {njM, dm4M}, FormQF]
 (* /. {p_Pair :> pair[p] /; FreeQ[p, q1],
        e_Eps :> eps[e] /; FreeQ[e, q1]} *)


qcM[x__, i___List] :=
  qc[FormQC[Times[x] /. {p_Pair :> pair[p], e_Eps :> eps[e]}], i]

qc[x_, ___List] := x /; LeafCount[x] < 10

qc[n_?NumberQ x_, i___List] := n qc[x, i]

qc[p_Plus, i___List] := -qc[-p, i] /; MatchQ[p[[1]], _?Negative _.]

qc[x_] := qc[x] = NewSymbol["QC"]

qc[x_, i_List] := Level[{{x}, Select[i, !FreeQ[x, #]&]}, {2}, qc]

qc[x_, i__] := qc[x, i] = NewSymbol["QC"][i]


addpatt[s_] := Pattern[#, _]&/@ s

unpatt[expr_] := expr /. Pattern -> (#1 &)


Attributes[dv] = {HoldAll}

dv[x_, _] := {} /; !FreeQ[x, Pattern]

dv[_, _?NumberQ] = {}

dv[_[_[x_, ___]], s_] := addpatt[s] -> x


Abbr[] := Flatten @ Apply[dv, DownValues/@
  {fermM, pair, eps, pol, abb, sunM, num, qc, xi}, {2}]

Abbr[patt__] :=
Block[ {all = Abbr[], need = Flatten[{patt}], omit},
  omit = Cases[need, !x_ :> x];
  need = DeleteCases[need, !_];
  FixedPoint[ abbsel[First/@ #]&, abbsel[need, omit] ]
]

abbsel[{need__}, {omit__}] := Select[all,
  !FreeQ[#, Alternatives[need]] && FreeQ[#, Alternatives[omit]]&]

abbsel[{need__}, ___] := Select[all, !FreeQ[#, Alternatives[need]]&]


Unabbr[expr_, args___] := expr //.
  Dispatch[abbrList[Flatten[{Subexpr[], Abbr[]}]][args]]

abbrList[abbr__][] := Flatten[{abbr}]

abbrList[abbr__][!patt_, r___] := abbrList[x:patt :> x, abbr][r]

abbrList[pre___, abbr_][patt_, r___] :=
  abbrList[pre, Select[abbr, FreeQ[#, patt]&]][r]


IndexHeader[h_, {}] := h

IndexHeader[h_, {i__}] := h[i]

IndexHeader[h_[i___], expr__] := IndexHeader[h, Flatten[{i,
  Select[Union @ Symbols[Level[{expr}, {-2}]],
    Head[Dim[#]] =!= DoDim &]}]]


MomEncode[other_] := other /; FreeQ[other, k]

MomEncode[p_Plus] := MomEncode/@ p

MomEncode[f_. k[i_]] := MomEncoding[f, i]


(* AbbrevInt introduces abbreviations for the loop integrals.
   They fall into two categories:
   1. A0, A00,
   2. B0i, C0i, D0i, E0i, F0i.
   For the latter the LoopTools functions [BCDEF]put are used to
   compute all tensor coefficients at once. *)

IntSym[f_, i___] := ToSymbol[f, i, ++intc[f]]

Function[ {X0i, Xput, Nxx},
  AbbrevInt[X0i[i_, args__]] :=
  Block[ {abb},
    abb[id_] = IndexHeader[IntSym[X0i][id], args];
    abbint[X0i[id_, args]] = abb[Epsi[id]];
    lint = {lint, Rule@@ {abb[Nxx_], Call[Xput[AddrOf[abb[1]], args]]}};
    abb[Epsi[i]]
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
    abbint[Xcut[args]] = abb = IndexHeader[IntSym[Xcut][], args];
    lint = {lint, abb -> MomEncode/@ Xcut[args]};
    abb
  ];
  AbbrevInt[Xmas[args__]] :=
  Block[ {abb},
    abb[id_] = IndexHeader[IntSym[Xmas][id], args];
    lint = {lint, Rule@@ {abb[Mxx_], Call[Xmas[AddrOf[abb[1]], args]]}};
    abbint[Xmas[args]] = abb[1]
  ] ]@@@
{ {Bcut, Bmas, Mbb},
  {Ccut, Cmas, Mcc},
  {Dcut, Dmas, Mdd},
  {Ecut, Emas, Mee},
  {Fcut, Fmas, Mff} }

AbbrevInt[func_] :=
Block[ {abb = IndexHeader[IntSym[Head[func], "i"][], func]},
  lint = {lint, abb -> func};
  abbint[func] = abb
]


ExtractInt[expr_] :=
Block[ {abbint, intc, new, lint = {}},
  abbint[f_] := AbbrevInt[f];
  _intc = 0;
  new = expr /. int:LoopIntegral[__] :> abbint[int];
  {Flatten[lint], new}
]


toGen = x:(G[_][_][__][__] | _Mass | _GaugeXi | VertexFunction[_][__]) :>
  gn[x /. Internal | External | Loop :> Sequence[]]

gn[G[_][cto_][fi__][kin__]] := gn[G[_][cto][fi][kin]] =
  Gsym["G", Gfi[fi], Gkin[kin]]

gn[Mass[fi_]] := gn[Mass[fi]] =
  Gsym["M", Gfi[fi]]

gn[GaugeXi[fi_]] := gn[GaugeXi[fi]] =
  Gsym["X", Gfi[fi]]

gn[VertexFunction[cto_][fi__]] := gn[VertexFunction[cto][fi]] =
  Gsym["V", cto, Gfi[fi]]


Gsym[h_, r___] := ToSymbol[h,
  DeleteCases[Level[{r}, {-1}, Heads -> True],
    s_Symbol /; Context[s] == "System`"]]


Gfi[fi__] := gfi/@ {fi} /. i_Index :> FromCharacterCode[Mod[Hash[i], 26] + 97]

gfi[(s:2 | -2) sv_[i__]] := gfi[s/2 (sv /. {SV -> VS, VS -> SV})[i]]

gfi[-fi_[i__]] := fi["b", i]

gfi[fi_] := fi


Gkin[kin__] := {kin} /. {
  KI1 | KI2 | KI3 | KI4 | KI5 | KI6 |
  SI1 | SI2 | SI3 | SI4 | SI5 | SI6 -> Identity,
  -Mom[i_] :> StringTake["KPQRUV", {i}],
  Mom[i_] :> StringTake["kpqruv", {i}],
  MetricTensor -> "g",
  DiracSlash | DiracMatrix -> "ga",
  ChiralityProjector[-1] -> "L",
  ChiralityProjector[+1] -> "R",
  FourVector | ScalarProduct | NonCommutative -> Dot }


GenericList[] := Flatten[dv@@@ DownValues[gn]]


Options[ClearProcess] = {
  ZapFunction -> Remove
}

ClearProcess[opt___Rule] :=
Block[ {zapf},
  {zapf} = ParseOpt[ClearProcess, opt];
  Apply[zap, DownValues/@ {fermM, pair, eps, pol, abb, sunM,
    num, qc, coup, mass, xi, vertex}, {2}];
  CurrentProc = CurrentOptions = Sequence[];
  DenList = {};
  Clear[SymbolNumber];
  _SymbolNumber = 0;
]


Attributes[zap] = {HoldAll}

zap[p_, s_Symbol] := (Unset[p]; zapf[s]) /; FreeQ[p, Pattern]

zap[p_, s_Symbol[__]] := (Unset[p]; zapf[s]) /; FreeQ[p, Pattern]


subdef[minleaf_, deny_, fuse_, pre_, ind_] := (
  DownValues[subabb] = {
    HoldPattern[subabb[x_]] :> (abbsub@@ Select[ind, !FreeQ[x, #]&])[pre[x]]
  } /. HoldPattern[abbsub@@ Select[{}, _]] -> abbsub[];

  DownValues[subcond] = {
    HoldPattern[subcond[x_]] :> LeafCount[x] > minleaf && FreeQ[x, deny]
  } /. HoldPattern[x_ && FreeQ[_, _[]]] :> x;

  DownValues[subfun] = {
    HoldPattern[subfun[x_]] :>
      Block[{try = subtry[subF[x], x]}, try /; try =!= False],
    HoldPattern[subfun[s_SumOver x_]] :> s subfun[x],
    HoldPattern[subfun[h_[x_]]] :> h[subfun[x]],
    HoldPattern[subfun[h_[x__]]] :> h[
      subfun[Select[h[x], subF[#] =!= True &], h],
      subfuse[Select[h[x], subF], h]
    ] /; Length[Intersection[Attributes[h], {Flat, Orderless}]] === 2,
    HoldPattern[subfun[x_]] :> subfun/@ x,
    HoldPattern[subfun[h_[x___], h_]] :> subfun/@ h[x],
    HoldPattern[subfun[x_, _]] :> subfun[x]
  } /. HoldPattern[subfuse[x_, _]] :> subadd[x] /; fuse;

  DownValues[AbbrevDo] = {
    HoldPattern[AbbrevDo[expr_, lev_Integer:0]] :>
      Replace[expr, Plus -> sublev, {lev, Infinity}, Heads -> True],
    HoldPattern[AbbrevDo[expr_, f_]] :>
      Block[{subF = f}, subfun[expr]]
  } /. sublev :> (pre[Plus[##]] /. Plus -> sublev &) /; pre =!= Identity
)


AbbrevSet[expr_, opt___Rule] :=
Block[ {minleaf, deny, fuse, pre},
  {minleaf, deny, fuse, pre} = ParseOpt[Abbreviate, opt];
  subdef[minleaf, Alt[deny], fuse, pre,
    Union[Cases[expr, SumOver[i_, ___] :> i, Infinity],
      LoopInd[], FormInd]];
]


AbbrevDo::unset = "Run AbbrevSet before AbbrevDo."

_AbbrevDo := (Message[AbbrevDo::unset]; Abort[])


Options[Abbreviate] = {
  MinLeafCount -> 10,
  Deny -> {k, q1},
  Fuse -> True,
  Preprocess -> Identity }

Abbreviate[expr_, x_:2, opt___Rule] :=
Block[ {subabb, subcond, subfun, AbbrevDo},
  AbbrevSet[expr, opt];
  AbbrevDo[expr, x]
]


sublev[x:(_Integer | _Rational) _., r__] :=
 (#1 subadd[Block[{plus = Plus}, #2]]&)@@
    FactorTermsList[Plus@@ ({x, r} /. Plus -> plus)]

sublev[p___] := subadd[Plus[p]]


subfuse[h_[x___], h_] := subadd/@ h[x]

subfuse[x_, _] := subadd[x]


subadd[x_^2 - y_^2] := subadd[x - y] subadd[x + y]

subadd[x_] := subabb[x] /; subcond[x]

subadd[x_] := x


subtry[False, _] := False

subtry[True, x_] := subabb[x] /; subcond[x]

subtry[f_[___], _] := False /; f === subF

subtry[fx_, x_] := (x /. fx -> subabb[fx]) /; subcond[fx]

_subtry = False


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
  _[_[_[x_]], s_] :> (addpatt[s] -> x) /; FreeQ[x, Pattern] ]


Options[ClearSubexpr] = {
  ZapFunction -> Remove
}

ClearSubexpr[opt___Rule] :=
Block[ {zapf},
  {zapf} = ParseOpt[ClearProcess, opt];
  zap@@@ SubValues[abbsub];
]


likeQ[patt_][s_] := StringMatchQ[ToString[s], patt]


abbcat[p_Pair] := {1, p, 1, 1, 1}

abbcat[e_Eps] := {1, 1, e, 1, 1}

abbcat[t_Pol] := {1, 1, 1, t, 1}

abbcat[d_DiracChain] := {1, 1, 1, 1, d}

abbcat[w_WeylChain] := {1, 1, 1, 1, w}

abbcat[other_] := {other, 1, 1, 1, 1}

regabb[s_ -> rhs_Times] :=
  Level[{{s}, Times@@@ ToCat[5, abbcat/@ List@@ rhs]}, {2}, regabb]

regabb[s_ -> rhs_] := Level[{{s}, abbcat[rhs]}, {2}, regabb]

regabb[s_, x_, 1, 1, 1, 1] := regabb[s, x]

  regabb[s_?(likeQ["SUN*"]), x_] := setabb[sunM, x, s]

  regabb[s_?(likeQ["Abb*"]), x_] := setabb[abb, x, s]

  regabb[s_?(likeQ["QC*"]), x_] := setabb[qc, x, s]

regabb[s_, 1, p_, 1, 1, 1] := setabb[pair, p, s]

regabb[s_, 1, 1, e_, 1, 1] := setabb[eps, e, s]

regabb[s_, 1, 1, 1, t_, 1] := setabb[pol, t, s]

regabb[s_, x_, 1, 1, 1, f_] := setabb[fermM, x f, s]

regabb[s_, rhs___] :=
  (Message[RegisterAbbr::unknown, #]; #)&[ s -> Times[rhs] ]


General::unknown = "Don't know how to register ``."

General::conflict = "Conflict of definitions for ``."

Attributes[RegisterAbbr] = {Listable}

RegisterAbbr[s_ -> x_Num] := setabb[num, x, s]

RegisterAbbr[s_ -> other_] := regabb[s -> 1 other]

RegisterAbbr[other_] := (Message[RegisterAbbr::unknown, other]; other)


setabb[h_, val_, s_] :=
  setval[h, val, s, Cases[DownValues[h], _[_[_[v_]], s] :> v]]

setsub[h_[i___], val_, s_] :=
  setval[h[i], val, s, Cases[SubValues[h], _[_[_[i][v_]], s] :> v]]

setval[h_, val_, s_, {}] := h[val] = s

setval[h_, val_, s_, {val_}] := s

setval[h_, val_, s_, _] := (Message[RegisterAbbr::conflict, s]; $Failed)


Attributes[RegisterSubexpr] = {Listable}

RegisterSubexpr[s_Symbol[i__] -> x_] :=
  setsub[abbsub[i], x, unpatt[s[i]]]

RegisterSubexpr[s_Symbol -> x_] :=
  setsub[abbsub[], x, s]

RegisterSubexpr[other_] :=
  (Message[RegisterSubexpr::unknown, other]; other)


$OptPrefix = "Opt"

newvar[] := NewSymbol[$OptPrefix]

newvar[i__] := NewSymbol[$OptPrefix][i]


Attributes[set] = {Flat, Orderless}


defT[x_, t_] := defL[x, unpatt[x], lot@@ t]

defP[x_, p_] := (
  defL[-x, -#, lot@@ -p];
  defL[x, #, lot@@ p] )& @ unpatt[x]

defL[xp_, x_, l:lot[t__]] := {l = x; xp, set[t], {}}

defL[xp_, x_, _Integer x_] := {{}, {}, xp -> 0}

defL[xp_, _, other_] := {{}, {}, xp -> other}


Attributes[comT] = Attributes[comP] = {Listable}

comT[a_set, b_set] := Intersection[a, b]

comP[a_set, b_set] :=
Block[ {cp, cm},
  cp = Intersection[a, b];
  If[ Length[a] > 2 Length[cp],
    cm = Intersection[a, Thread[-b, set]];
    If[ Length[cm] > Length[cp], cp = cm ] ];
  cp
]

_comT = _comP = set[]


rulP[a_, b_] := {a -> b, Thread[-a, set] -> -b}


Overlap[] = set[]

Overlap[x__] := Intersection[x]


_CSE[] = {}

CSE[h_, defH_, comH_, rulH_, simp_][abbr__] :=
Block[ {lot, alias, var, def, com, tmp, new = {}},
  Attributes[lot] = {Flat, Orderless};
  {var, def, alias} = ToCat[3, defH@@@ SortBy[{abbr}, Length[ #[[2]] ]&]];
  def = def /. alias;

  Do[
    While[ Length[com = comH@@ Take[def, {i, i + 1}]] > 3,
      tmp = Ceiling[Length[com]/2];
      tmp = Overlap@@ Select[ comH[com, Drop[def, {i, i + 1}]],
        Length[#] > tmp & ];
      If[ Length[tmp] > 3, com = tmp ];
      tmp = Position[def, com, 1, 1, Heads -> False];
      If[ Length[tmp] > 0,
        tmp = tmp[[1,1]];
        def[[tmp,0]] = List;
        def = def /. rulH[com, var[[tmp]] //unpatt];
        def[[tmp,0]] = set,
      (* else *)
        tmp = newvar@@ Select[ Union[Level[Take[var, {i, i + 1}], {2}]],
          !FreeQ[com, unpatt[#]]& ];
        new = {new, tmp -> h@@ com};
        def = def /. rulH[com, unpatt[tmp]] ]
    ],
  {i, Length[def] - 1}];

  {simp[new], Thread[var -> h@@@ def], alias}
]


AbbrCat[rul:_[_, _Plus]] := {{}, {}, rul}

AbbrCat[rul:_[_, t_Times]] := {{}, rul, {}} /; FreeQ[t, DiracChain]

AbbrCat[rul_] := {rul, {}, {}}


OptimizeAbbr[{}, ___] = {}

OptimizeAbbr[rul:{__Rule}, simp_:Simplify] := Flatten @
  MapThread[Apply, {
    { List,
      CSE[Times, defT, comT, Rule, simp],
      CSE[Plus, defP, comP, rulP, simp] },
    ToCat[3, AbbrCat/@ rul] }]


oneSubst[exprlist_] :=
Block[ {subst, new, subs},
  new = select[exprlist];
  subs = #[[1,1,1]]&/@ DownValues[subst];
  If[ Length[subs] === 0, Throw[new] ];
  allsubs = {allsubs, subs};
  new //. subs
]

nestSubst[exprlist_] := Catch[Nest[oneSubst, exprlist, 10]]

SubstAbbr[exprlist_, patt_, deny___] :=
Block[ {select, glob = Null, rm, substAbbr, allsubs = {}, subvars, pos},
  Attributes[select] = {Listable};
  Scan[(select[r:(# -> _)] := r)&, {deny}];
  select[r:(var_ -> _)] := r /; !FreeQ[glob, var];
  select[r:patt] := subst[r] = Sequence[];
  select[CodeExpr[vars_, tmps_, expr_]] :=
    CodeExpr[rm[vars], rm[tmps], select[expr]];
  select[i_IndexIf] := MapIf[substAbbr, i];
  select[other_] := other;
  subvars := subvars = Union[First/@ Flatten[allsubs]];
  NestWhile[
    ( glob = ReplacePart[#, Null &, pos] /. CodeExpr -> (#3 &);
      ReplacePart[#, nestSubst, pos] )&,
    nestSubst[exprlist],
    Length[pos = Position[#, substAbbr, Infinity, 1]] > 0 &
  ] /. rm[vars_] :> Complement[vars, subvars]
]

P$SimpleAbbr =
  (_ -> _Symbol | _?NumberQ | _?NumberQ _Symbol) |
  ((lhs_ -> rhs_) /; LeafCount[lhs] >= LeafCount[rhs])

SubstSimpleAbbr[exprlist_, deny___] :=
  SubstAbbr[exprlist, P$SimpleAbbr, deny]


(* UV and IR finiteness checks *)

ToNewBRules = {
  B0[args__] :> B0i[bb0, args],
  B1[args__] :> B0i[bb1, args],
  B00[args__] :> B0i[bb00, args],
  B11[args__] :> B0i[bb11, args],
  B001[args__] :> B0i[bb001, args],
  B111[args__] :> B0i[bb111, args],
  DB0[args__] :> B0i[dbb0, args],
  DB1[args__] :> B0i[dbb1, args],
  DB00[args__] :> B0i[dbb00, args],
  DB11[args__] :> B0i[dbb11, args],
  C0[args__] :> C0i[cc0, args],
  D0[args__] :> D0i[dd0, args],
  E0[args__] :> E0i[ee0, args],
  F0[args__] :> F0i[ff0, args] }

ToOldBRules = {
  B0i[bb0, args__] :> B0[args],
  B0i[bb1, args__] :> B1[args],
  B0i[bb00, args__] :> B00[args],
  B0i[bb11, args__] :> B11[args],
  B0i[bb001, args__] :> B001[args],
  B0i[bb111, args__] :> B111[args],
  B0i[dbb0, args__] :> DB0[args],
  B0i[dbb1, args__] :> DB1[args],
  B0i[dbb00, args__] :> DB00[args],
  B0i[dbb11, args__] :> DB11[args] }


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
Block[ {div = RCPattern[RenConst, Divergence], foo = f, Finite = fin},
  FindDiv[expr /. ToNewBRules /.
    {int:PaVeIntegral[__] :> UVDiv[int] + fin int, D -> Dminus4 + 4} /.
    Dminus4^(n_?Negative) :> (-2/Divergence)^n]
]


UVDivergentPart[expr_] := MapDiv[SubDiv, expr, 0]


UVSeries[expr_, pow_] := Coefficient[UVSeries[expr], Dminus4, pow]

UVSeries[expr_] := MapDiv[SerDiv, expr, 1]


UVDiv[A0[m_]] := m Divergence

UVDiv[A00[m_]] := m^2/4 Divergence

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


(* matrix elements *)

SelectAbbr[HelicityME][abbr_, patt_, loop_, 1] :=
  abbsq[HelicityME] @ inamp[loop] @ Select[abbr, !FreeQ[#, patt]&]

SelectAbbr[h_][abbr_, patt_, loop_, tree_] := abbsq[h]@@
  Through[{inamp[loop], inamp[tree]} @ Select[abbr, !FreeQ[#, patt]&]]

inamp[All] = Identity

inamp[amp_][abbr_] := Select[abbr, !FreeQ[amp, First[#]]&]


General::nomat = "Warning: No matrix elements to compute."

abbsq[h_][{}] := (Message[h::nomat]; {})

_abbsq[loop_] := {loop, {1}}

abbsq[h_][{}, {}] := (Message[h::nomat]; {})

_abbsq[{}, tree_] := {tree, tree}

_abbsq[loop_, {}] := {loop, loop}

abbsq[_][amp__] :=
  Thread[If[{$TreeSquare, $LoopSquare}, Hold[amp], {amp}]] /. Hold -> Union


ConjF[f_ -> expr_, {1, _}] := f[ToString[f] -> expr]

ConjF[f_ -> expr_, {2, _}] := f[ToString[f] <> "C" ->
  expr /. DiracChain -> ConjChain /. ep_Eps :> -ep /. ConjWF]


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
ConjOm[_, 1] = {1, 1};
ConjOm[_, -1] = {1, -1}


FormLor[n_] := FormLor[n] = "N" <> ToString[n] <> "_?"


HelicitySq[f_[s_ -> _], 1] := Mat[f][ToStr["Mat", f] -> s]

HelicitySq[f_[s_ -> _], fc_[sc_ -> _]] :=
  Mat[f, fc][ToStr["Mat", f, fc] -> s <> "*" <> sc]


HelTab[Hel[i_]] := helM[i, e]

HelTab[h_] := h + e


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
  TreeSquare :> $TreeSquare,
  LoopSquare :> $LoopSquare,
  Dimension -> 4,
  EditCode -> False,
  RetainFile -> False }

HelicityME::noprocess = "No process defined so far.  \
HelicityME works only after DeclareProcess or CalcFeynAmp."

HelicityME::weyl = "Warning: HelicityME does not work on WeylChains.  \
CalcFeynAmp uses DiracChains with the option FermionChains -> Chiral or VA."

HelicityME[tree_, opt___?OptionQ] := HelicityME[tree, tree, opt]

HelicityME[tree_, loop_, opt___?OptionQ] :=
Block[ {dim, edit, retain,
abbr, ferm, vars, hels, helM, hh, e, mat},
  If[ CurrentProc === {},
    Message[HelicityME::noprocess];
    Abort[] ];
  {$TreeSquare, $LoopSquare, dim, edit, retain} =
     ParseOpt[HelicityME, opt] /. Options[CalcFeynAmp];

  abbr = Abbr[];
  If[ !FreeQ[abbr, WeylChain], Message[HelicityME::weyl] ];

  abbr = Check[
    SelectAbbr[HelicityME][abbr, DiracChain[_Spinor, ___], loop, tree],
    Return[{}] ];

  ferm = Union[Cases[abbr, Spinor[k[i_], __] :> i, Infinity]];

  abbr = MapIndexed[ConjF, abbr, {2}] /. Lor -> FormLor /.
    Reverse/@ FromFormRules /. Eps -> "e_" /.
    FinalFormRules;

  mat = Outer[HelicitySq, Sequence@@ abbr] //Flatten;

  FCPrint[1, "> ", Length[mat], " helicity matrix elements"];

  hels = HelTab/@ Hel/@ ferm;
  dim = If[dim === 0, D, 4];
  abbr = Level[abbr, {3}];
  vars = FormVars[dim, {Last/@ abbr, hels}];

  hh = OpenForm["fc-hel-"];
  WriteString[hh, "\
#define FermionChains \"None\"\n\n" <>
    vars[[1]] <> "\n\
table HEL(" <> ToString[Min[ferm]] <> ":" <>
               ToString[Max[ferm]] <> ", e?);\n" <>
    MapThread[{"fill HEL(", ToForm[#1], ") = ", ToForm[#2], ";\n"}&,
      {ferm, hels}] <> "\n" <>
    FormProcs <>
    FormCode["Common.frm"] <>
    FormCode["HelicityME.frm"]];

  FormWrite[hh, abbr];

  WriteString[hh, ".sort\ndrop;\n\n"];

  FormWrite[hh, Level[mat, {2}]];

  WriteString[hh, "#call Emit\n"];
  Close[hh];

  (e[#] = s[#])&/@ ferm;
  helM[i_, s_] := Hel[i] + s;

  FormPre[mat];

  Thread[Head/@ mat -> Plus@@@ FormOutput[][edit, retain]]
]


WeylSq[fac_][f_ -> _, fc_ -> _] :=
  Mat[f, fc] -> f Conjugate[fc] fac[f, fc]


Options[WeylME] = {
  TreeSquare :> $TreeSquare,
  LoopSquare :> $LoopSquare,
  MatFactor -> (1 &)
}

WeylME::weyl = "Warning: WeylME does not work on DiracChains."

WeylME[tree_, opt___?OptionQ] := WeylME[tree, tree, opt]

WeylME[tree_, loop_, opt___?OptionQ] :=
Block[ {fac, abbr},
  ChkProc[{tree, loop}, WeylME];

  {$TreeSquare, $LoopSquare, fac} = ParseOpt[WeylME, opt];

  abbr = Abbr[];
  If[ !FreeQ[abbr, DiracChain], Message[WeylME::dirac] ];
  abbr = Check[
    SelectAbbr[WeylME][abbr, WeylChain[_Spinor, ___], loop, tree],
    Return[{}] ];

  Outer[WeylSq[fac], Sequence@@ abbr] //Flatten
]


(* colour matrix elements *)

SUNN = 3

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

sunTsum[i_, j_, k_, l_] :=
  (sunT[i, l] sunT[k, j] - sunT[i, j] sunT[k, l]/SUNN)/2

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

Conjugate[f_SUNF] ^:= f


sunOver[i:Index[Colour, _], n_:SUNN, ___] := sunSum[i, n]

sunOver[i:Index[Gluon, _], n_:(SUNN^2 - 1), ___] := sunSum[i, n]

sunOver[i__] := SumOver[i]


Attributes[sunExpand] = {Listable}

sunExpand[expr_] := expr /; FreeQ[expr, sunSum]

sunExpand[expr_. sunSum[i_, n_]] := sunExpand[n expr] /; FreeQ[expr, i]

sunExpand[expr_ sunSum[i_, _]] := sunExpand[expr /. i -> iname@@ i]

sunExpand[expr_] := sunExpand/@ expr


ColourSimplify[tree_:1, loop_] :=
Block[ {res},
  res = sunExpand[Conjugate[tree] loop /.
    {IndexDelta -> idelta, SumOver -> sunOver}];
  res = Expand[res /. 
    {SUNT -> sunText, SUNF -> sunF,
     SUNEps -> sunEps, SUNTSum -> sunTsum} /.
    sunTrace[a__]^n_. :> Times@@ Table[sunTr[a], {n}]];
  Simplify[res] /. {
    sunT[i_, j_] :> Sort[SUNT[i, j]],
    sunT[g__, i_, i_] :> SUNT[g, 0, 0],
    sunT -> SUNT, sunEps -> SUNEps}
]


ColourGrouping[tops_] := DiagramGrouping[ tops,
  Replace[
    ColourSimplify[Times@@
      FeynAmpCases[_[Index[Colour | Gluon, _], ___]][##]],
    _?NumberQ r_ :> r]& ]


ColourSq[f_ -> x_, fc_ -> xc_] := Mat[f, fc] -> ColourSimplify[x, xc]


Options[ColourME] = {
  TreeSquare :> $TreeSquare,
  LoopSquare :> $LoopSquare
}

ColourME[tree_, opt___?OptionQ] := ColourME[tree, tree, opt]

ColourME[tree_, loop_, opt___?OptionQ] :=
Block[ {abbr},
  {$TreeSquare, $LoopSquare} = ParseOpt[ColourME, opt];

  abbr = Check[
    SelectAbbr[ColourME][Abbr[], SUNObjs, loop, tree],
    Return[{}] ];

  Outer[ColourSq, Sequence@@ abbr] //Flatten
]


(* squaring the matrix element *)

UniqIndices[tree_, loop_] :=
Block[ {ind},
  ind = Intersection@@
    (Union[Cases[#, SumOver[i_, ___] :> i, Infinity]]&)/@ {tree, loop};
  tree /. Thread[ind -> (ToSymbol[#, "c"]&)/@ ind]
]


SquaredME[amp_] := SquaredME[amp, amp]

SquaredME[tree:Amp[_][___], loop:Amp[_][___]] := (
  ChkProc[{tree, loop}, SquaredME];
  squaredME[Plus@@ tree, Plus@@ loop]
)

squaredME[0, _] = squaredME[_, 0] = 0

squaredME[tree_, loop_] :=
Block[ {ffdef = {}, ff, ffc},
  {ff, ffc} = (#1/@ ToList[Collect[#2, _Mat, FF]])&@@@ {
    matFF[FF, Identity] -> loop,
    matFF[FFC, Conjugate] -> UniqIndices[tree, loop] };
  { TermCollect[Plus@@ (FF[#] Plus@@ matSq[#]/@ ffc &)/@ ff],
    Flatten[ffdef] }
]

matSq[1][mc_] := FFC[mc]

matSq[m_Times][mc_] := FFC[mc] Inner[Mat, m, mc, Times]

matSq[m_][mc_] := FFC[mc] Mat[m, mc]


matFF[h_, c_][FF[f_] Mat[m_]] := (ffdef = {ffdef, h[m] -> c[f]}; m)

matFF[h_, c_][FF[f_]] := (ffdef = {ffdef, h[1] -> c[f]}; 1)


Unprotect[Conjugate]

Format[Conjugate[x_]] := Superscript[x, "*"]

(*
Format[Conjugate[t_Times]] := Superscript[SequenceForm["(", t, ")"], "*"]
*)

Conjugate[D] = D

Conjugate[Dminus4] = Dminus4

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

Conjugate[(x_?RealQ)^n_.] := x^n

Conjugate[HoldCode[x_]] := HoldCode[Conjugate[x]]

Protect[Conjugate]


Unprotect[Re]

Re[D] = D

Re[Dminus4] ^= Dminus4

Re[Finite] ^= Finite

Re[Divergence] ^= Divergence

Re[x_. o_SumOver] := Re[x] o

Re[d_IndexDelta] := d

Re[p_Plus] := Re/@ p

Re[d_Den] := Re/@ d

Re[(x_?RealQ)^n_.] := x^n

Re[(x_?RealQ)^n_. y_] := x^n Re[y]

Protect[Re]


RealQ[_Integer] = True

RealQ[_Real] = True

RealQ[p_Plus] := VectorQ[List@@ p, RealQ]

RealQ[t_Times] := VectorQ[List@@ t, RealQ]


(* performing the polarization sum analytically *)

Options[PolarizationSum] = {
  SumLegs -> All,
  Dimension -> 4,
  GaugeTerms -> True,
  NoBracket -> NoBracket,
  EditCode -> False,
  RetainFile -> False }

PolarizationSum::noprocess = "No process defined so far.  \
PolarizationSum works only after DeclareProcess or CalcFeynAmp."

PolarizationSum::incomp = "PolarizationSum used on an amplitude \
other than the last one set up by DeclareProcess or CalcFeynAmp."

PolarizationSum::polchain =
"Warning: Input contains polarization vectors in fermion chains.  \
Result likely not correct, please insert the HelicityME first."

PolarizationSum[amp:Amp[_][___].., opt___?OptionQ] :=
Block[ {Hel},
  ChkProc[{amp, CurrentProc}, PolarizationSum, Abort[]];
  _Hel = 0;
  PolarizationSum[
    { $HelicityME = HelicityME[amp, FilterOpt[HelicityME, opt]],
      $ColourME = ColourME[amp, FilterOpt[ColourME, opt]],
      SquaredME[amp] },
    opt ]
]

PolarizationSum[expr_, opt___?OptionQ] :=
Block[ {slegs, dim, gauge, nobrk, edit, retain,
fexpr, lor, indices, legs, masses, e, ferm = {},
etasubst, vars, hh},

  If[ CurrentProc === {},
    Message[PolarizationSum::noprocess];
    Abort[] ];

  {slegs, dim, gauge, nobrk, edit, retain} =
    ParseOpt[PolarizationSum, opt] /. Options[CalcFeynAmp];

  fexpr = Flatten[{expr}];
  fexpr = fexpr /. Cases[fexpr, (lhs_ -> _) :>
    lhs -> "\\[" <> ToString[lhs, CForm] <> "\\]"];
  fexpr = Unabbr[fexpr, !_Mat, DiracChain | WeylChain] /.
    FinalFormRules;
  If[ !FreeQ[Unabbr[fexpr], (DiracChain | WeylChain)[__, _e | _ec, __]],
    Message[PolarizationSum::polchain] ];

  lor = Cases[fexpr, _Lor, Infinity] //Union;
  indices = FormIndices[[ Level[lor, {2}] ]];
  fexpr = fexpr /. Thread[lor -> indices];

  legs = Cases[fexpr, Alt[ExtWF][i_] :> i,
    Infinity, Heads -> True] //Union;
  If[ slegs =!= All, legs = Intersection[legs, Flatten[{slegs}]] ];
  masses = Masses[CurrentProc][[legs]];

  fexpr = fexpr /. s[i_] :> (ferm = {ferm, i}; e[i]) /.
    Reverse/@ FromFormRules /.
    {Eps -> "e_", MetricTensor -> "d_", Pair -> Dot} /.
    NoExpandRule /.
    FinalFormRules;

  etasubst = Block[{dv = DownValues[eta], eta},
    Cases[dv, _[_[lhs_], rhs_] :> (lhs -> rhs)] /. Reverse/@ FromFormRules];

  dim = If[dim === 0, D, 4];
  vars = FormVars[dim, {fexpr, masses}, indices];

  hh = OpenForm["fc-pol-"];
  WriteString[hh, "\
#define Dim \"", ToString[dim], "\"\n\
#define GaugeTerms \"" <> ToString[gauge] <> "\"\n\
#define FermionChains \"None\"\n\n" <>
    vars[[1]] <> "\n\
#procedure EtaSubst\n" <>
    FormId[etasubst] <> "\
#endprocedure\n" <>
    FormProcs <>
    FormConst[vars, nobrk] <>
    FormCode["Common.frm"] <>
    FormCode["PolarizationSum.frm"]];

  If[ MemberQ[fexpr, _Rule],
    FormWrite[hh, fexpr];
    WriteString[hh, "\
#call Prepare\n\
#define Prepared\n\
.sort\n\
drop;\n\n"] ];
  Write[hh, "L SquaredME = ", Plus@@ DeleteCases[fexpr, _Rule], ";"];

  WriteString[hh,
    "\n#call eiei\n" <>
    MapThread[{"\n#call PolSum(", ToString[#1], ", ", ToForm[#2], ", ",
        ToString[If[FreeQ[fexpr, (z | zc)[#1]], dim, Dminus4]], ")"}&,
      {legs, masses}] <>
    "\n\n#call Emit\n"];
  Close[hh];

  nobrk = Alt[nobrk];	(* for DotSimplify *)
  (e[#] = s[#])&/@ Union[Flatten[ferm]];
  FormPre[fexpr];

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

MkDir[dirs__] := MkDir[Flatten[{dirs}]]


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

  CopyFile[ToFileName[$FormCalcBin, #], ToFileName[path, #]]&/@
    {"util.a", "simd.h"};

  path
]


(* Fortran code generation *)

Dim[i_Integer] := i

Dim[i_] := DoDim[i]

LoopInd[] := #[[1,1,1]]&/@ Select[DownValues[DoDim], FreeQ[#, Pattern]&]


FunctionNames[base_, ind___] := {
  #1 <> abridge[0, #2, Infinity],
  #1 <> abridge[0, #2, $MaxFunctionName -
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


InvDef[h_[i_], h_[j_]] := Invariant[1, i, j] -> SInvariant[i, j]

InvDef[_[i_], _[j_]] := Invariant[-1, i, j] -> TInvariant[i, j]


InvList[n_, r__] := {InvDef[n, #]&/@ {r}, InvList[r]}

InvList[_] = {}


pdefs[proc_] :=
Block[ {qn, n = 0},
  qn = (##4&)@@@ Level[proc, {2}];
  qn = Union[Flatten[DeleteCases[qn, _?NumberQ, Infinity]]];
  { "\n\n\
#undef Compose\n\
#define Compose(f,c", #1, ") ",
    Fold[{"c(", #1, ",", #2, ")"}&, First[#2], Rest[#2]], 
    MapThread[ {"\n\n\
#undef ", #1, "\n\
#define ", #1, "(f,c) Compose(f,c", {",", #}&/@ #2, ")"}&,
      {Flatten[{"Generic", "Anti", "Mass", ToCode/@ qn}], {##3}} ]
  }&@@ Transpose[Level[MapIndexed[pspec, proc, {2}], {2}]]
]

pspec[{(s_Integer:1) p_, _, m_, q_:0}, {o_, _}] := Flatten[{
  "," <> #,
  "f(" <> # <> "," <> ToString[n] <> "," <> ToString[o] <> ")",
  ptype[p, m],
  ToString[s],
  ToCode[m],
  pqnum[Plus@@ Flatten[{q}]]/@ qn
}]&[ FromCharacterCode[++n + 64] ]

ptype[_V, 0] := "PHOTON";
ptype[_V, _] := "VECTOR";
ptype[_F, _] := "UFERMION";
ptype[_F, _] := "FERMION";
ptype[_S, _] := "SCALAR";
ptype[_U, _] := "GHOST"

pqnum[expr_][qn_] := ToCode[Coefficient[expr, qn]]


phel[hmax_, {h_, _, {_. _S | _U, ___}}] :=
  (Sow[h -> 0]; hmax);
phel[hmax_, {h_, h_, {_. _F, ___}}] := dfCode[
  (Sow[h -> "hval_df"@@ h]; hmax),
  (Sow[h -> "hval"[hmax, 2, 2, 1]]; 2 hmax) ];
phel[hmax_, {h_, h_, {_. _V, _, 0, ___}}] :=
  (Sow[h -> "hval"[hmax, 2, 2, 1]]; 2 hmax);
phel[hmax_, {h_, h_, {_. _V, ___}}] :=
  (Sow[h -> "hval"[hmax, 3, 1, 1]]; 3 hmax);
phel[hmax_, {h_, x_, _}] :=
  (Sow[h -> x]; hmax)


ProcCheck[p_] := (
  proc = p;
  name = ToString[Map[First, p, {2}], InputForm];
  legs = Plus@@ Length/@ p;
  invs = If[ legs < 4, {},
    InvList@@ MapIndexed[Apply, Drop[Signs[p], -1]] //Flatten ];
  $SymbolPrefix = StringJoin[ProcSubst[symprefix]];
)

ProcCheck[p_, p_] := 0

_ProcCheck := Message[WriteSquaredME::incomp]


Attributes[ProcSubst] = {Listable}

ProcSubst[s_String] := s

ProcSubst[f_] := f[proc]


DefModName[dir_] := (
  ModName[mod_, ext_:$CodeExt] := ModName[mod, ext] =
    ToFileName[dir, file = prefix <> mod <> ext];
  header = StringReplace[header, "%t" -> TimeStamp[]];
  Hdr[desc_] := StringReplace[header, {"%f" -> file, "%d" -> desc}];
)

MakefileName[file__] := {" $(DIR)", prefix, file}

Attributes[MakefileObjs] = {Listable}

MakefileObjs[mod_] := {" \\\n ", MakefileName[mod, ".o"]}


FFPut[Amp[p_][amp__], array_, file_] := (
  ProcCheck[p, proc];
  FFWrite[#, array, file]&/@ {amp}
)

FFPut[_[], __] = {}

FFPut[other_, __] := (Message[WriteSquaredME::noamp, other]; Abort[])


FFWrite[0, __] = {}

FFWrite[amp_, array_, file_] :=
Block[ {ind, ff, mods},
  ind = Union[Cases[amp,
    SumOver[i_, r_] :> (inds = {inds, i}; DoDim[i] = r; i), Infinity]];
  ff = amp /. unused[array] -> 0 /. xrules /. {
    _SumOver -> 1,
    int:LoopIntegral[__] :> abbint[int] };
  ff = FFList[ff, array];
  mods = FileSplit[ff, FunctionNames[file, ind], delay[FFMod]];
  (Indices[#] = ind)&/@ (#[[2,2]] &)/@ mods;
  mods
]


FFMod[ff_, {fmod_, mod_}] :=
Block[ {hh},
  FCPrint[3, "  writing ", mod];
  hh = OpenCode[ModName[fmod]];
  WriteString[hh,
    Hdr["form factors for " <> name] <>
    vPre <> fPost <>
    SubroutineDecl[mod] <>
    vPre];
  WriteExpr[hh, {sPost, ff, sEnd, SubroutineEnd[]},
    TmpType -> dupType,
    Optimize -> True,
    DebugLines -> $DebugFF, DebugLabel -> fmod];
  Close[hh];
  {fmod, mod}
]

$DebugFF = 1


FFList[0, _] = {}

FFList[amp_, array_] := (
  maxmat[array] = {Mat[1] -> 1};
  RuleAdd[array[1], amp] ) /; FreeQ[amp, Mat]

FFList[amp_Plus, array_] := FFList[#, array]&/@ List@@ amp

FFList[Mat[m_] x_., array_] := (
  maxmat[array] = {maxmat[array], #};
  RuleAdd[Level[#, {3}, array], x]
)&[ (ToArray[#] -> #)&/@ Level[m, {-1}] ]


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


setdep[Tag[___, r_]] := setdep[r]

setdep[v_ -> _] := v = TAG[++c]

markdep[Tag[___, r_]] := markdep[r]

markdep[v_ -> x_] := (v = TAG[++c]; {}) /; !FreeQ[x, TAG]

markdep[r_] := {} /; !FreeQ[r, TAG]


MoveDepsRight[li_List] := {li}

MoveDepsRight[f__List, li_List] :=
Block[ {pos, c = 0, cc = 0},
  Block[ #,
    setdep/@ li;
    pos = Map[markdep, {f}, {2}];
    While[c != cc, cc = c; pos = pos]
  ]&[ Union@@ Map[Kind, {f, li}, {2}] ];
  pos = Position[pos, {}, {2}, Heads -> False];
  { Sequence@@ MoveDepsRight@@ Delete[{f}, pos],
    JoinDeps[li, Extract[{f}, pos]] }
]


JoinDeps = Join


depcat[expr_][r:Tag[___, v_ -> _]] := (cpos = 1; {r, {}}) /; FreeQ[expr, v]

depcat[expr_][r:(v_ -> _)] := (cpos = 1; {r, {}}) /; FreeQ[expr, v]

depcat[_][r_] := (cnew = 1; {{}, r})


MoveDepsLeft[li_List] := {li}

MoveDepsLeft[li_List, f__List] :=
Block[ {pos = {f}, old = {}, new = li, cpos = 1, cnew = 1},
  While[ cpos + cnew > 1,
    old = {new, old};
    cpos = cnew = 0;
    {pos, new} = Transpose[ ToCat[2, depcat[new]/@ #1]&/@ pos ];
    new = Flatten[new] ];
  { JoinDeps[new, Flatten[old]],
    Sequence@@ MoveDepsLeft@@ pos }
]


toTicket[Tag[___, x_]] := toTicket[x]

toTicket[v_ -> x_] := (v = Dep[v]; Ticket[Dep[v], x]) /; FreeQ[v, Pattern]

toTicket[v_ -> x_] := {
  Ticket[Dep[v], x /. v -> 1],
  v = Dep[v] /. HoldPattern[Pattern][a_, _] :> a }[[1]]

toTicket[x_] := Ticket[x]

fromTicket[v_, x_] := v -> x

fromTicket[x_] := x

opo[li_] :=
Block[ {c = 0, l = Length[li], Dep, Ticket, posmap, prev},
  Attributes[Dep] = {HoldFirst};
  Block[{##, AddrOf},
    _AddrOf = 4711;
    posmap = Hold@@ Evaluate[toTicket/@ li]
  ]&@@ Union[Kind/@ li];
  Ticket[v_, x_] := (v = Random[]; ++c) /; FreeQ[x, Dep];
  Ticket[x_] := ++c /; FreeQ[x, Dep];
  While[ c < l,
    prev = c;
    posmap = Evaluate/@ posmap;
    If[ c === prev,
      $OnePassDebug = Join[$OnePassDebug,
        Cases[posmap, t_Ticket :> fromTicket@@ t] /. Dep -> Identity];
      Message[OnePassOrder::recurs, $OnePassDebug];
      Return[li] ]
  ];
  li[[ Level[Sort[MapIndexed[List, posmap]], {3}] ]]
]

$OnePassDebug = {}

OnePassOrder::recurs = "Recursive definition among `` \
(full list in $OnePassDebug).  Returning list unordered."

OnePassOrder[li_List] :=
Block[ {jf, jd},
  $OnePassDebug = {};
  Attributes[jf] = {HoldFirst};
  jf[v_ -> IndexIf[a__]] := jf[ v, {a}[[1;;;;2]], {a}[[2;;;;2]] ];
  jf[v_, cond_, x_] := (
    jf[w_, cond, y_] := (jd[cond] = {jd[cond], w -> y}; {});
    jd[cond] = v -> x;
    jx[cond] );
  jf[other_] := other;
  opo[ Flatten[jf/@ li] /. jx[cond_] :> Level[
    Delete[Transpose @
      {cond, opo/@ Transpose[Thread/@ Flatten[{jd[cond]}]]}, {-1, 1}],
    {2}, IndexIf] ]
]

OnePassOrder[i_IndexIf] :=
Block[ {defs, com, lcom, rcom},
  defs = (List@@ i)[[2;;;;2]];
  com = First/@ Intersection@@ defs;
  defs = fromIf/@ defs;
  lcom = OnePassOrder[Intersection@@ First/@ defs];
  rcom = OnePassOrder[Intersection@@ Last/@ defs];
  Flatten[{lcom, Level[
    MapThread[toIf[lcom, rcom], {(List@@ i)[[1;;;;2]], defs}],
    {2}, IndexIf ], rcom}]
]

comQ[def:_[(Tag[___, var_] | var_), _]] := {def, {}} /; MemberQ[com, var]

comQ[def_] := {{}, def}

fromIf[li_] :=
Block[ {l, m, r},
  {l, m} = ToCat[2, comQ/@ li];
  {l, m, r} = Block[{JoinDeps = Sequence}, MoveDepsRight[l, m]];
  {m, r} = MoveDepsLeft[m, r];
  {l, m, r}
]

toIf[lmin_, rmin_][cond_, {l_, m_, r_}] := {cond,
  OnePassOrder[Flatten[{Complement[l, lmin], m, Complement[r, rmin]}]]}


FindDeps[li_, patt_] :=
Block[ {rul, p = patt, old = {}},
  rul = depsow@@@ Cases[li, _Rule, Infinity];
  While[ Length[ new = Reap[rul = rul /. p -> Indeterminate][[2]] ] > 0,
    old = {old, new};
    p = Alt[new] ];
  Union[Flatten[old]]
]

depsow[var_, val_] := (Sow[var]; 0) /; !FreeQ[val, Indeterminate]


_AbbrMod[{}, _] := Sequence[]

AbbrMod[cat_][abbr_, mod_] :=
Block[ {hh},
  FCPrint[3, "  writing ", mod];
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["abbreviations for " <> name] <>
    vPre <> fPost <>
    SubroutineDecl[mod] <>
    vPre ];
  WriteExpr[hh, {sPost, abbr, sEnd},
    TmpType -> tmpType[cat],
    DebugLines -> $DebugAbbr[cat], DebugLabel -> mod];
  WriteString[hh, SubroutineEnd[]];
  Close[hh];
  mod
]

tmpType["h"] := helType

_tmpType = "ComplexType"

_$DebugAbbr = 0


NumMod[{}, _] := Sequence[]

NumMod[expr_, mod_] :=
Block[ {hh},
  FCPrint[3, "  writing ", mod];
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["numerators for " <> name] <>
    vPre <> fPost ];
  WriteNum[hh]/@ expr;
  Close[hh];
  mod
]


NumName[h_[___], ___] := _h -> $SymbolPrefix <> ToCode[h]

NumName[h_, ___] := h -> $SymbolPrefix <> ToCode[h]


NumExpr[name_][sub_, {{0} -> muexp_}] :=
  {"MuExp", name, sub, njcoeff[HelAll[1]] -> muexp}

NumExpr[name_][subexpr_, coefflist:{i_ -> _, ___}] :=
Block[ {type, coeff, res = {}, sub = subexpr},
  type = {"MuExp", "T3Exp", "T2Exp"}[[ Length[i] ]];
  coeff[j_ -> c_, {n_}] := (
    res = {res,
      type <> "Coeff(" <> ToSeq[j] <> ")\n",
      Select[sub, !FreeQ[ c, #[[1]] ]&],
      njcoeff[HelAll[n]] -> c};
    sub = Select[sub, FreeQ[ c, #[[1]] ]&] );
  MapIndexed[coeff, coefflist];
  {type, name, Flatten[res]}
]

NumExpr[name_][sub_, num_] :=
  {"Num", name, sub, "Result(" <> name <> ")" -> num}


WriteNum[hh_][var_ -> Num[expr__]] :=
Block[ {$AbbPrefix = Sequence["qx", 0]},
  SymbolNumber["qx"] = 0;
  writeNum[hh]@@ NumExpr[ NumName[var][[2]] ]@@
    Subexpr[{expr},
      MatchQ[#, _WeylChain | _Pair | _Eps]&,
      Deny -> {}, Fuse -> False, MinLeafCount -> 2] /.
    Pair[q1, q1] -> Pair[q1, q1] - MuTildeSq
]

writeNumF[hh_][type_, name_, expr__] := (
  WriteString[hh, "\
\t" <> type <> "Function(" <> name <> ")\n\
\timplicit none\n"];
  WriteExpr[hh,
    { vPre <>
      "#define C" <> type <> "\n\n" <>
      sPost <>
      "#undef C" <> type <> "\n\n\
!" <> type <> "Debug(\"" <> name <> "\")\n\n",
      expr,
      sEnd,
      "\tend\n\n\n" },
    Type -> helType,
    FinalTouch -> Simplify,
    FinalCollect -> True,
    DebugLines -> $DebugNum, DebugLabel -> name]
)

$DebugNum = 0


$SymbolPrefix = ""


VarDecl[vars__] := StringJoin[(varArgs@@ ({vars} /. (x_ -> _) :> x))[]]

varArgs[{}, _, r___] := varArgs[r]

varArgs[{l__}, r___] := varArgs[l, r] /; !FreeQ[{l}, List, {2, Infinity}]

varArgs[s_String, r___][d___] := varArgs[r][d, s]

varArgs[vars_List, type_, r___][d___] := varArgs[r][d, chkDecl[
  DeleteCases[
    Replace[unpatt[vars],
      i_Symbol :> Dim[i], {2, Infinity}] /. DoDim -> Identity,
    (*_[0] |*) _[] ],
  type ]]

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


chkDecl[{}, _] := Sequence[]

chkDecl[vars_, type_] := {vars, type}


Unprotect[Span];
Format[i_ ;; j_, FortranForm] := SequenceForm[i, ":", j];
Protect[Span]

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
  {"struct v", com, " {\n", varMapC[com]@@@ vars, "} v", com, "_;\n"}

varDeclC[vars_, type_String] := {
  "  ", type, " ", StringReplace[StringTake[ToString[
    Replace[vars, h_[i_, j__] :>
      Fold[#1[#2]&, h, RanDim/@ Reverse[{i, j}]], 1],
    InputForm, PageWidth -> 63 - StringLength[type] ], {2, -2}],
  ", \n" -> ";\n  " <> type], ";\n" }

varDeclC[vars_, decl_, com_String] := {
  "struct ", com, " {\n", decl, "} ", $SymbolPrefix, com, "_;\n",
  {"\n#define ", #1, " ", $SymbolPrefix, com, "_.", #2}&@@@
    varSubstC/@ vars,
  "\n" }


varMapF[com_][vars_, type_] :=
Block[ {off = 1, v = com <> StringTake[type, 1]},
  arr = {arr, v};
  { varMap[v <> "(", ")\n"]/@ vars,
    "#define ", n = v <> "Len", " ", ToDef[off - 1], "\n",
    varDeclF[{v[n]}, type] }
]

varMapC[com_][vars_, type_] :=
Block[ {off = 0, v = StringTake[type, 1], n},
  { varMap[com <> "." <> v <> "[", "]\n"]/@ vars,
    "#define ", n = com <> v <> "Len", " ", ToDef[off], "\n",
    varDeclC[{v[n]}, type] }
]

varMap[arrL_, arrR_][var_[i___]] :=
Block[ {ind = off, stride = 1, c = 104, lhs, HelDim = Identity},
  lhs = Reap[
    ( ind += stride Sow[FromCharacterCode[++c]] - RanOff[#] stride;
      stride *= RanDim[#] )&/@ {i} ];
  off += stride;
  {"#define ", ToDef[Level[lhs, {3}, var]], " ", arrL, ToDef[ind], arrR}
]

varMap[arrL_, arrR_][var_] :=
  {"#define ", ToCode[var], " ", arrL, ToCode[off++], arrR}


varSubstC[var_[i___]] :=
  {varSeq[#1, "(", #2, ",", ")"], varCarr[##]}&[
    ToCode[var],
    Array[FromCharacterCode, Length[{i}], 105], {i} ]

varSubstC[var_] := {#, #}& @ ToCode[var]


varSeq[___, {}, ___] = {}

varSeq[s1___, {i___, j_}, c_, s2___] := s1 <> ({#, c}&/@ {i}) <> j <> s2


varCarr[h_, i__] := h <> Reverse[MapThread[varCind, {i}]]


varCind[i_, n_Symbol] := {"[", i, "]"} /; Context[n] === "LoopTools`"

varCind[i_, n_ ;; _] := {"[", i, "-", ToString[n], "]"}

varCind[i_, n_] := {"[", i, "-1]"}


ToFunction[h:_[___]] := h

ToFunction[h_] := h[]


subroutineDeclF[name_, decl___String] :=
  "\tsubroutine " <> $SymbolPrefix <> ToCode[name] <>
    "\n\timplicit none\n" <> decl <> "\n"

subroutineDeclC[name_, decl___String] :=
  "void " <> $SymbolPrefix <> ToCode[ToFunction[name]] <> " {\n" <>
    decl <> "\n"


newline[""] = newline[{}] = ""

newline[s_] := s <> "\n"


IfDecl[var_String, if_, else_] := {
  "#ifndef " <> var <> "\n#define " <> var <> "\n\n" <> if <> "#else\n\n",
  else, "#endif\n\n" }

IfDecl[_, decl_, expr_] := {decl, expr}

ifDeclF[var_, com_, fdecl_, sdecl_, x___] := IfDecl[var, fdecl, {sdecl, com, x}]

ifDeclC[var_, com_, fdecl_, sdecl_, ___] := IfDecl[var, {fdecl, com}, sdecl]


CallDecl[li_List] := StringJoin[CallDecl/@ li]

CallDecl[DoLoop[name_, ind__]] :=
  ("\n" <> #1 <> CallDecl[name] <> #2 &)@@ DoDecl[ind]

CallDecl[name_] := callDecl[name]

callDeclC[name_] :=  "  " <> $SymbolPrefix <> name <> "();\n"

callDeclF[name_] := "\tcall " <> $SymbolPrefix <> name <> "\n"


DoDecl[{var_}] := DoDecl[{var, DoDim[var]}]

DoDecl[{_, _DoDim}] := {{}, {}}

DoDecl[{var_, Span[i__]}] := DoDecl[{var, i}]

DoDecl[{var_, from_:1, to_, step_:1}] := {
  "!LOOP(" <> ToCode[var] <> ", " <> ToCode[{from, to, step}] <> ")\n",
  "!ENDLOOP(" <> ToCode[var] <> ")\n" }

DoDecl[var_] := DoDecl[{var, DoDim[var]}]

DoDecl[vars__] := {StringJoin[#1], StringJoin[Reverse[#2]]}&@@
  Transpose[DoDecl/@ ReverseDo[{vars}]]

ReverseDo = (*Identity*) Reverse


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

MatType[h_[i_], h_[j_]] := MatSym[h][i, j]

MatType[h_[i_], 0] := MatSym[h][i, 0]

MatType[0, h_[j_]] := Conjugate[MatSym[h][j, 0]]

MatSym[h_] := MatSym[h] = ToSymbol["Mat", h]


defMat[mat:(m_[imin_ ;; imax_ | imax_, jmin_ ;; jmax_ | jmax_])] := (
  matsel[_[m[i_, j_], _]] := imin <= i <= imax && jmin <= j <= jmax;
  matsub[_m, _] := (matsub[_m, _] =.; mat) )

defTree[m_[imin_ ;; imax_ | imax_, jmin_ ;; jmax_ | jmax_]] :=
  defcat[r:_[m[i_, j_], _]] := {r, {}, {}} /;
    imin <= i <= imax && jmin <= j <= jmax

defLoop[m_[imin_ ;; imax_ | imax_, jmin_ ;; jmax_ | jmax_]] :=
  defcat[r:_[m[i_, j_], _]] := {{}, r, {}} /;
    imin <= i <= imax && jmin <= j <= jmax


matnan[m_[i__], {h_, ___}] :=
  {ToString[m], RanDim/@ {i}, {"HelNaN", "NaN"}[[h]]}


Assort[m_Mat, x_] := {m -> x, {}, {}}

Assort[v_, n_Num] := {{}, v -> n, {}}

Assort[v_, x_] := {{}, {}, v -> x} /; FreeQ[x, DiracChain | SUNT]

_Assort := Sequence[]


ChainHead[om_, n_] := ChainHead[om, n] =
  ToSymbol["Chain", {"V", "B"}[[om - 5]], n]

SplitChain[Spinor[_[i_], _, s1_, d1_, e1_], om_, g___,
    Spinor[_[j_], _, s2_, d2_, e2_]] :=
  s1^(om - e1) s2^(om + Length[{g}] + e2 + 1) *
  ChainHead[om, Length[{g}]][Spinor[i, s1, d1], e1, g, e2, Spinor[j, s2, d2]]


Attributes[Cond] = {HoldRest}

Cond[True, args__] := {args}

Cond[False, __] = {}


Attributes[helDim] =
Attributes[helDimF] =
Attributes[helDimC] = {HoldFirst}

helDimF[h_[i_, j___], ___] :=
  (Sow[h[a_, b___] :> h[HelAll[a], b]]; h[HelDim[i], j])

helDimF[h_, ___] := h


helDimC[h_Hel, ___] := (Sow[x_Hel :> hel[x]]; h)

helDimC[h_[i___], ___] := (Sow[x_h :> vec[x]]; h[i])

helDimC[h_, ___] := (Sow[h :> vec[h]]; h)


vcat[vec[v_]] := {{}, {}, {}, {v}}

vcat[hel[h_]] := {{}, {}, {h}, {}}

vcat[n_?NumberQ] := {{n}, {}, {}, {}}

vcat[other_] := {{}, {other}, {}, {}}


vtimes[a__] := Times[a] /; FreeQ[{a}, vec | hel]

vtimes[a__] := vt@@ Apply[Times, Transpose[vcat/@ {a}], {1, 2}]

vt[n_, 1, h_, 1] := hel[n h]

vt[n_, c_, h_, 1] := vec[SxI[n c, h]]

vt[n_, c_, 1, v_] := vec[SxH[n c, hxh[v]]]

vt[n_, c_, h_, v_] := vec[HxH[SxI[n c, h], hxh[v]]]

hxh[v_ w__] := Fold[HxH, v, {w}]

hxh[v_] := v

SxI[1, h_] := ItoH[h]

SxH[1, h_] := h

HxH[ItoH[h_], v_] := IxH[h, v]


vplus[a__] := Plus[a] /; FreeQ[{a}, vec | hel]

vplus[a__] := vp@@ Apply[Plus, Transpose[vcat/@ {a}], {1, 2}]

vp[n_, 0, h_, 0] := hel[n + h]

vp[n_, c_, 0, v_] := vec[StoH[n + c] + v]

vp[n_, 0, h_, v_] := vec[ItoH[n + h] + v]

vp[n_, c_, h_, v_] := vec[StoH[n + c] + ItoH[h] + v]

StoH[0] = 0

ItoH[0] = 0


vec/: vec[x_]^y_ := vec[x^y]

vec/: vec[x_] -> y_?NumberQ := vec[x] -> ToH[y]

vec/: RuleAdd[vec[x_], y_?NumberQ] := RuleAdd[vec[x], ToH[y]]

vec/: Conjugate[vec[x_]] := vec[ConjugateH[x]]

vec/: Re[vec[x_]] := vec[ReH[x]]


ToH[0] = "HelZero"

ToH[n_] := StoH[n]


comlim[defs_, end_] :=
  ToCode[{Kind[First[Flatten[{defs, end -> 0}]]], end}]


fDecl[{s___}] := fDecl[s]

fDecl[pre_String:"", post___String] :=
  {newline[pre], newline[StringJoin[post]] <> "\n"}

sDecl[{s___}] := sDecl[s]

sDecl[pre_String:"", post_String:"", end_String:""] :=
  {newline[pre], newline[post], newline[end]}


Attributes[WriteSquaredME] = {HoldAll}

Options[WriteSquaredME] = {
  TreeSquare :> $TreeSquare,
  LoopSquare :> $LoopSquare,
  Folder -> "squaredme" (* {"squaredme", ProcName} *),
  ExtraRules -> {},
  FilePrefix -> "",
  SymbolPrefix -> "",
  FileHeader -> "#if 0\n* %f\n* %d\n* generated by " <>
    $FormCalcVersion <> " on %t\n#endif\n\n",
  FileIncludes -> {"#include \"decl.h\"\n",
                   "#include \"inline.h\"\n",
                   "#include \"contains.h\"\n"},
  SubroutineIncludes -> FileIncludes
}

WriteSquaredME::noamp = "`` is not an amplitude."

WriteSquaredME::empty = "Warning: no amplitudes were specified."

WriteSquaredME::badmat = "Warning: Incompatible matrix elements `` and ``."

WriteSquaredME[tree_, tree_, r___] := WriteSquaredME[tree, {}, r]

WriteSquaredME[tree_, loop_, dir_, opt___?OptionQ] :=
  WriteSquaredME[tree, loop, Abbr[], dir, opt]

WriteSquaredME[tree_, loop_, abbr__, dir_, opt___?OptionQ] :=
Block[ {folder, xrules, prefix, symprefix, header, fincl, sincl,
vPre, fPre, fPost, sPre, sPost, sEnd, dfCode, ModName, Hdr, 
proc = Sequence[], name, legs, invs, $SymbolPrefix,
mat, nums, abrs, matsel, matsub, defcat, angledep,
abbint, intc, lint = {},
inds = {}, defs, Indices, pos, file, files, hh,
unused, maxmat, mats, mmat, nmat, mat1, ntree, nloop,
ffmods, nummods, abbrmods,
com, helrul, helvars, dupType, hmax, hfun},

  {$TreeSquare, $LoopSquare, folder, xrules, prefix, symprefix,
    header, fincl, sincl} =
    ParseOpt[WriteSquaredME, opt] //. Options[WriteRenConst];

  vPre = "#include \"" <> prefix <> "vars.h\"\n";
  {fPre, fPost} = fDecl[fincl];
  {sPre, sPost, sEnd} = sDecl[sincl];

  abbint[f_] := AbbrevInt[f];
  _intc = 0;

  Attributes[dfCode] = {HoldAll};
  If[ FreeQ[{abbr}, DiracChain],
    dfCode[_, x___] := x,
    dfCode[x_, ___] := x ];

  {mat, nums, abrs} = ToCat[3, Assort@@@ Flatten[{abbr}]];

  mats = First/@ DeleteCases[mat, _ -> 0];
  unused[Ctree] = Alt@@
    Select[Union[#[[1,2]]&/@ mat], FreeQ[mats, Mat[_, #]]&];
  unused[Cloop] = Alt@@
    Select[Union[#[[1,1]]&/@ mat], FreeQ[mats, Mat[#, _]]&];

(* Part 1: the form factors *)

  FCPrint[2, "formatting form factors"];

  _maxmat = {};
  ffmods = Flatten/@ {
    Block[{modnum = 0}, WriteFF[tree, Ctree]],
    Block[{modnum = 0}, WriteFF[loop, Cloop]] };
  If[ Plus@@ Length/@ ffmods === 0,
    Message[WriteSquaredME::empty]; Return[{}] ];

  abrs = abrs /. xrules /. int:LoopIntegral[__] :> abbint[int];
  lint = RhsMap[Prio[1], Flatten[lint]];

  mmat = Split[Union[Flatten[#]], #1[[1,0]] === #2[[1,0]]&]&/@
    {maxmat[Cloop], maxmat[Ctree]};
  {nloop, ntree} = nmat = Map[#[[1,1,0]][Level[#, {3}, Ran]]&, mmat, {2}];
  mat1 = {};
  If[ Length[ntree] === 0, ntree = nloop,
    If[ Length[nloop] === 0, nloop = ntree,
      nloop = MaxDims[nloop, ntree];
      If[ $LoopSquare, ntree = nloop ];
      If[ UnsameQ@@ Map[Head, nmat, {2}],
        Message[WriteSquaredME::badmat, nmat[[2]], nmat[[1]]];
        mat1 = Complement[##, SameTest -> (#1[[1,1,0]] === #2[[1,1,0]]&)]&@@@
          {mmat, Reverse[mmat]};
        mat1 = Map[
          Block[ {h = #[[1,1,0]]},
            ntree = MaxDims[ntree, h[0]];
            h = MatSym[h];
            (h[#1[[1]], 0] -> #2)&@@@ # ]&, mat1, {2} ];
      ];
  ] ];
  mats = Select[MapThread[MatType, {nloop, ntree}], Length[#] > 0 &];

(* Part 2: the numerators and abbreviations *)

  FCPrint[2, "formatting abbreviations"];

  Scan[defMat, mats];
  mat = (MatType@@ ToArray/@ #1 -> #2)&@@@ mat;
  mat = Select[mat, matsel];

  defs = Flatten[{abrs, mat, mat1, nums, lint}];

  (* split into tree/loop *)
  Scan[ defTree[MatType[#, #]]&, nmat[[2]] ];
  Scan[ defLoop[MatType[#, #]]&, Complement@@ nmat ];
  defcat[r:_[v_, _]] := {r, {}, {}} /; !FreeQ[ffmods[[1]], v];
  defcat[r:_[v_, _]] := {{}, r, {}} /; !FreeQ[ffmods[[2]], v];
  defcat[r_] := {{}, {}, r};
  defs = Join[#1, Tag/@ #2]&@@ MoveDepsLeft@@ ToCat[3, defcat/@ defs];

  (* split into s/angle/hel *)
  angledep = Alt[(Range[#2] + #1)&@@ Length/@ proc];
  angledep = Alt[{
    Cases[invs, _[x_, r_] :> x /; MemberQ[r, angledep, {1}]],
    (k | s)[angledep],
    MomEncoding[_, angledep] }];
  defs = OnePassOrder/@ VarSort/@ MoveDepsRight@@
    ToCat[3, Category/@ defs];
  defs = Append[
    {#1, #2} /. {Pair -> Pair0, Eps -> Eps0, k -> k0, s -> s0},
    #3 ]&@@ defs;
  pos = Take[#, 2]&/@ Position[defs, _Num];
  nums = Extract[defs, pos] /. Tag -> Identity;
  defs = ToCat[2, #]&/@ Replace[ Delete[defs, pos],
    {Tag[r_] :> {{}, r}, r_ :> {r, {}}}, {2} ];

  nummods = FileSplit[nums, "num", delay[NumMod]];
  nums = NumName@@@ nums;

  abbrmods = MapThread[
    FileSplit[ToDoLoops[#1], #2[[2]] <> #2[[1]], delay[AbbrMod[ #2[[1]] ]]]&,
    {defs, Outer[List, {"s", "a", "h"}, {"abbr0", "abbr1"}]},
    2 ] /. nums;

  _matsub = {};
  mats = Apply[matsub, Reverse[defs], {3}];
  mats = {Flatten[#1], Flatten[{##2}]}&@@ mats;

  ff = ToCat[2, MapThread[
    If[FreeQ[#1, Alt[First/@ #2]], {{}, #3}, {#3, {}}]&,
    {ffmods, defs[[3]], {Level[nmat[[2]], {2}, Ctree],
                         Level[nmat[[1]], {2}, Cloop]}} ]];

  abrs = First/@ Join[abrs, lint];
  defs = Map[Select[#, MemberQ[abrs, First[#]]&]&, defs, {2}];

  helvars = Alt[kind/@ First/@ Flatten[ defs[[3]] ]];
  _dupType = If[FreeQ[#, helvars], "ComplexType", helType]&;

  {com, helrul} = Reap[helDim/@ Array[Hel, legs]; VarDecl[
    { Common["varXs"][defs[[1,1]], "ComplexType",
                      defs[[1,2]], "ComplexType"],
      Common["varXa"][defs[[2,1]], "ComplexType",
                      defs[[2,2]], "ComplexType",
                      invs, "RealType"],
      Common["varXh"][helDim@@@ defs[[3,1]], helType,
                      helDim@@@ defs[[3,2]], helType] },
    Common["indices"][Union[Flatten[inds]], "integer"],
    Common["formfactors"][
      helDim/@ Join[ ff[[1]], mats[[1]] ], helType,
      Join[ ff[[2]], mats[[2]] ], "ComplexType" ],
    NotEmpty["#ifndef NUM_H\n",
      Last/@ nums, Extern,
      "#endif\n\n"] ]];
  helrul = Dispatch[Flatten[helrul] //Union];
  mats = Level[MapIndexed[matnan, mats, {2}], {2}];

  FCPrint[2, "writing code modules"];

  DefModName[MkDir[dir, ProcSubst[folder]]];

  {ffmods, abbrmods, nummods} = {ffmods, abbrmods, nummods} /.
    helrul /.
    vxrul /.
    vec | hel -> Identity /.
    delay -> Identity;

(* Part 3: the variable declarations *)

  FCPrint[2, "writing variable declarations"];

  hh = OpenCode[ModName["vars", ".h"]];
  WriteString[hh, Hdr["variable declarations"] <>
    ifDecl["vars_h", com, "\
#define SQUAREDME\n\
#define LEGS " <> ToString[legs] <> "\n\n" <> fPre, sPre]];
  Close[hh];

(* Part 4: the makefile *)

  FCPrint[2, "writing makefile"];

  hh = OpenWrite[ModName["makefile", ""]];

  WriteString[hh, "\
LIBS += $(LIB)\n\n\
OBJS :=" <>
    MakefileObjs[{nummods, abbrmods, Apply[#1&, ffmods, {2}], 
      "SquaredME"}] <> "\n\n\
$(LIB): $(LIB)($(OBJS))\n\n\
$(LIB)($(OBJS)): " <> $MakeDeps[[1]] <> MakefileName["vars.h"] <> "\n\n"];

  Close[hh];

  ffmods = Apply[#2&, ffmods, {2}];

(* Part 6: the process specs *)

  FCPrint[2, "writing specs.h"];

  hh = OpenFortran[ModName["specs", ".h"]];

  WriteString[hh, (Hdr["process specifications"] <> "\
#undef FCVERSION\n\
#define FCVERSION \"" <> $FormCalcVersion <> "\"\n\n\
#undef PROCNAME\n\
#define PROCNAME \"" <> ProcName[proc] <> "\"\n\n\
#undef SQUAREDME_FUNC\n\
#define SQUAREDME_FUNC " <> $SymbolPrefix <> "SquaredME\n\n\
#undef KIN\n\
#define KIN \"" <> #1 <> "to" <> #2 <> ".F\"\n\n\
#undef IDENTICALFACTOR\n\
#define IDENTICALFACTOR " <>
    ToString[Times@@ (Factorial[Length[#]]&)/@
      Split[DeleteCases[First/@ proc[[2]], _Symbol, {-1}]//Sort]] <> "\n\n\
#undef DIRACFERMIONS\n\
#define DIRACFERMIONS " <> dfCode["1", "0"] <>
    pdefs[proc] <> "\n\n")&@@ ToString/@ Length/@ proc];

  Close[hh];

(* Part 5: the master subroutine SquaredME *)

  FCPrint[2, "writing SquaredME routine"];

  {maxmat[Cloop], maxmat[Ctree]} = LoopReduce[nmat] /. helrul;

  proc = Level[proc, {2}];
  {hmax, hfun} = Reap[Fold[phel, 1, Reverse @ Transpose @
    {Array[Hel0, legs], Array[Hel, legs] /. Hel -> Hel0, proc}]];
  hfun = ToCode/@ OnePassOrder[Reverse@@ hfun];

  hh = OpenCode[ModName["SquaredME"]];
  writeSquaredME[hh, proc, hfun, hmax];
  Close[hh];

  Cases[DownValues[ModName], _[_, s_String] :> s]
]


	(* LoopElem[Ctree, j] gives back e.g.
1. Fortran: Ctree(HelAll(jF),jSUN),
   C: Ctree(jF,jSUN)
2. {F[jF], SUN[jSUN]}
3. "\n\tLOOP(jF, 1,3,1)
    \n\tLOOP(jSUN, 1,5,1)"
4. "\n\tENDLOOP(jSUN)
    \n\tENDLOOP(jF)" *)

LoopVar[_][h_[-1]] = {h[1], "", ""}

LoopVar[j_][h_[n_]] := LoopVar[j][h[1 ;; n]]

LoopVar[j_][h_[m_ ;; n_]] := ( ffind = {ffind, #};
  { h[#],
    {"\n", $CodeIndent, "LOOP(", #, ", ", ToString[m], ",", ToString[n], ",1)"},
    {"\n", $CodeIndent, "ENDLOOP(", #, ")"} }
)&[ j <> ToString[h] ]


LoopElem[ff_, j_] := LoopElem[maxmat[ff], ff, j]

LoopElem[{}, __] = ""

LoopElem[mm_, ff_, j_] := {Level[#1, {2}, ff] /. helrul, ##}&@@
  Transpose[LoopVar[j]/@ mm]


mpset[h_[n_], {i_, ___}] := (mp[h, i] = h[n]; h)


sumup[l1_, l2_] := sumup[l2, l1] /; Length[ l1[[2]] ] > Length[ l2[[2]] ]

sumup[{ff1_, v1_, l1_, e1_}, {ff2_, v2_, l2_, e2_}] :=
Block[ {mp, v2v1, a = vx["amp"]},
  _mp = 0;
  v2v1 = Times@@ (MatType[mp[#, 1], mp[#, 2]]&)/@
    Union[MapIndexed[mpset, {v2, v1}, {2}] //Flatten] /. helrul;
  { l1,
    $CodeBoln, ToCode[a -> 0], $CodeEoln,
    l2,
    $CodeBoln, ToCode[RuleAdd[a, ff2 v2v1 /. vxrul]], $CodeEoln,
    e2,
    $CodeBoln, ToCode[RuleAdd["ampsq", Re[Conjugate[ff1] a /. vxrul]]], $CodeEoln,
    e1 }
]

_sumup = {}


Attributes[vxCode] = {HoldAll}


writeSquaredMEF[hh_, part_, hfun_, hmax_] :=
Block[ {jtree, jloop, ffcode, ffind = {}},
  jtree = LoopElem[Ctree, "j"];
  jloop = LoopElem[Cloop, "j"];
  ffcode = "\n\
\tampsq = 0" <> ({"\n\
* ", prefix, "BEGIN FF_TREE",
    #3, "\n\t", ToCode[#1 -> 0], #4, "\n\n",
    CallDecl[ToDoLoops[ffmods[[1]], Indices]], "\
* ", prefix, "END FF_TREE\n\n\
* ", prefix, "BEGIN M2_TREE",
    Cond[ TrueQ[$TreeSquare],
      sumup[jtree, LoopElem[Ctree, "i"]] ], "\n\
* ", prefix, "END M2_TREE"
  }&)@@ jtree <> "\n\
\tres(" <> vxCode["HelAll(1)", "1"] <> ") = ampsq\n\n\
\tampsq = 0" <>
  ({"\n\
\tTEST(flags, BIT_LOOP)\n\
* ", prefix, "BEGIN FF_LOOP",
    #3, "\n\t", ToCode[#1 -> 0], #4, "\n\n",
    CallDecl[ToDoLoops[ffmods[[2]], Indices]], "\
* ", prefix, "END FF_LOOP\n\n\
* ", prefix, "BEGIN M2_LOOP",
    Cond[ jtree =!= "",
      sumup[jtree, LoopElem[Cloop, "i"]],
      "\n\tampsq = ampsq + ampsq" ],
    Cond[ jtree === "" || TrueQ[$LoopSquare],
      sumup[jloop, LoopElem[Cloop, "i"]] ], "\n\
* ", prefix, "END M2_LOOP\n\
\tENDTEST(flags, BIT_LOOP)"
  }&)@@ jloop <> "\n\
\tres(" <> vxCode["HelAll(2)", "2"] <> ") = ampsq\n\n";

  WriteString[hh, "\
*#define CHECK\n\n" <>
    Hdr["assembly of squared matrix element"] <> "\
#include \"" <> prefix <> "specs.h\"\n" <>
    vPre <> fPost <> "\
************************************************************************\n\n\
\tsubroutine " <> $SymbolPrefix <> "SquaredMEHel(" <>
      vxCode["HelInd(vmax,res)", "res"] <> ", flags)\n\
\timplicit none\n" <>
    vxCode["\
\tSIMD_ONLY(integer vmax)\n"] <> "\
\t" <> resType <> vxCode[" res(HelDim(*))", " res(*)"] <> "\n\
\tinteger flags\n\n" <>
    vPre <> "\n\
\t" <> helType <> " amp\n\
\t" <> resType <> " ampsq\n" <>
    VarDecl[Union[Flatten[ffind]], "integer"] <>
    vxCode[
      $DebugPre[1, 3] <> "\
\tSIMD_ONLY(integer v)\n" <>
      $DebugPost[1]
    ] <> "\n\
* " <> prefix <> "BEGIN ABBR_HEL\n" <>
    Invoke@@ abbrmods[[3]] <> "\
* " <> prefix <> "END ABBR_HEL\n" <>
    ffcode <>
    $DebugPre[1, 3] <>
    vxCode["\
\tprint 1, HelLoop(Hel(HelInd(v,1:LEGS)), res(HelInd(v,1:2)), v,vmax)\n",
    (* vxElse *) "\
\tprint 1, Hel(1:LEGS), res(1:2)\n"
    ] <> "\
1\tformat(' sqme(', LEGS I3, ') =', 2F25.12)\n",
    $DebugPost[1] <> "\
\tend\n\n\
************************************************************************\n\n\
\tsubroutine " <> $SymbolPrefix <> "SquaredME(result, helicities, flags)\n\
\timplicit none\n\
\tRealType result(*)\n\
\tinteger*8 helicities\n\
\tinteger flags\n\n" <>
    vPre <> "\n\
* " <> prefix <> "BEGIN VAR_DECL\n\
\tinteger*8 hbits, hlast\n\
\tinteger seq(2)\n\
\tsave hlast, seq\n\
\tinteger i, h, hmax\n\
\tparameter (hmax = " <> ToString[hmax] <> ")\n" <>
    vxCode["\
\tSIMD_ONLY(integer v)\n\
\tinteger hsimd\n\
\tparameter (hsimd = SIMD_CEIL(hmax))\n\
\tResType res(HelDim(2),hsimd)\n",
    (* vxElse *) "\
\tRealType res(2,hmax)\n"
    ] <> "\
\tRealType rtree, rloop\n\
\texternal " <> $SymbolPrefix <> "SquaredMEHel\n\
* " <> prefix <> "END VAR_DECL\n\n\
* " <> prefix <> "BEGIN HSEL_DECL\n\
\tRealType norm\n\
\tRealType hseltest(0:hmax-1)\n" <>
    vxCode["\
\tResType hseltest_v(HelDim(hsimd))\n\
\tequivalence (hseltest_v, hseltest)\n"
    ] <> "\
\tRealType hselmin\n\
\tinteger hseli\n\
\tsave hseltest, hselmin, hseli\n\
* " <> prefix <> "END HSEL_DECL\n\n" <>
    ({"\tdata ", #1, " /", ToCode[#3[Times@@ #2]], "/\n"}&@@@ mats) <> "\n" <>
    sPost <> "\
* " <> prefix <> "BEGIN SETMASS\n\
\tTEST(flags, BIT_SETMASS)" <>
    MapIndexed[{"\n\tresult(", ToString[ #2[[1]] ], ") = ",
      ToCode[ #1[[3]] ]}&, part] <> "\n\
\treturn\n\
\tENDTEST(flags, BIT_SETMASS)\n\
* " <> prefix <> "END SETMASS\n\n"];

  WriteExpr[hh, {"\
* " <> prefix <> "BEGIN INVARIANTS\n",
    invs, "\
* " <> prefix <> "END INVARIANTS\n\n"}, Newline -> ""];

  WriteString[hh, "\
\tCHK_INI(seq)\n\n\
\tTEST(flags, BIT_RESET)\n\
\thlast = 0\n\
\tseq(1) = seq(1) + 1\n\
\tINI_S()\n\
* " <> prefix <> "BEGIN ABBR_S\n" <>
    Invoke@@ abbrmods[[1]] <> "\
* " <> prefix <> "END ABBR_S\n\
\tENDTEST(flags, BIT_RESET)\n\n\
\tseq(2) = seq(2) + 1\n\
\tINI_A()\n\
* " <> prefix <> "BEGIN ABBR_ANGLE\n" <>
    Invoke@@ abbrmods[[2]] <> "\
* " <> prefix <> "END ABBR_ANGLE\n\n\
* " <> prefix <> "BEGIN HEL_LOOP\n\
\thelicities = iand(helicities, Generic(ARG_Phel,JOIN_HEL))\n\n" <>
    vxCode["\
\tSIMD_ONLY(v = 1)\n"
    ] <> "\
\th = 0\n\
\tdo i = 0, hmax - 1\n\
* " <> prefix <> "BEGIN HSEL_IF\n\
\t  if( hseltest(i) .lt. hselmin ) cycle\n\
* " <> prefix <> "END HSEL_IF\n\n\
#define hval(s,m,a,b) a*mod(i/s,m)-b\n" <>
    dfCode["\
#define hval_df(i) HelBit(helicities,i,1) - HelBit(helicities,i,-1)\n"
    ] <>
    Thread[{"\n\t  ", hfun}] <>
    If[ hmax === 1, "\n", "\n\
\t  hbits = Generic(ARG_Bhel,JOIN_HEL)\n\
\t  if( iand(hbits, helicities) .ne. hbits ) cycle\n"
    ] <>
    $DebugPre[1, 1] <> "\
\t  if( helicities .ne. hlast ) print 1, ' helicities: ', Hel0\n" <>
    $DebugPost[1] <>
    vxCode["\
\t  SIMD_ONLY(call VecCopy(v, LEGS))\n\
\t  SIMD_ONLY(v = mod(v, SIMD) + 1)\n\
\t  SIMD_ONLY(if( v .eq. 1 ) then)\n\
\t  h = h + 1\n\
\t  call " <> $SymbolPrefix <> "SquaredMEHel(HelInd(SIMD,res(HelAll(1),h)), flags)\n\
\t  SIMD_ONLY(endif)\n\
\tenddo\n\n\
\tSIMD_ONLY(if( v .ne. 1 ) then)\n\
\tSIMD_ONLY(h = h + 1)\n\
\tSIMD_ONLY(call " <> $SymbolPrefix <> "SquaredMEHel(HelInd(v,res(HelAll(1),h)), flags))\n\
\tSIMD_ONLY(endif)\n\
\tDEINI()\n\
\tSIMD_ONLY(if( v .ne. 1 ) res(v:SIMD,:,h) = 0)\n",
    (* vxElse *) "\
\t  h = h + 1\n\
\t  call " <> $SymbolPrefix <> "SquaredMEHel(res(1,h), flags)\n\
\tenddo\n\n\
\tDEINI()\n"
    ] <> "\
* " <> prefix <> "END HEL_LOOP\n\n\
* " <> prefix <> "BEGIN RESULT\n\
\trtree = 0\n\
\trloop = 0\n\
\tdo i = 1, h\n" <>
    vxCode["\
\t  rtree = rtree + HelSum(res(HelAll(1),i))\n\
\t  rloop = rloop + HelSum(res(HelAll(2),i))\n",
    (* vxElse *) "\
\t  rtree = rtree + res(1,i)\n\
\t  rloop = rloop + res(2,i)\n"
    ] <> "\
\tenddo\n" <>
  dfCode["\
\thbits = iand(helicities, Generic(ARG_Ferm,JOIN_HEL))\n\
\ti = ibset(0, BitCount(iand(hbits, ishft(hbits, -2))))\n\
\trtree = rtree*i\n\
\trloop = rloop*i\n"
  ] <> "\
\tresult(1) = rtree\n\
\tTEST(flags, BIT_LOOP)\n\
\tresult(2) = rloop\n\
\tENDTEST(flags, BIT_LOOP)\n\n" <>
    $DebugPre[1, 2] <> "\
\tprint *, PROCNAME, ' =', rtree, rloop\n" <>
    $DebugPost[1] <> "\
* " <> prefix <> "END RESULT\n\n\
* " <> prefix <> "BEGIN HSEL_SET\n\
\tif( helicities .ne. hlast ) then\n\
\t  hseltest = 0\n\
\t  hselmin = 0\n\
\t  hseli = 0\n\
\tendif\n\
\tif( hseli .lt. hseln .and. rtree + rloop .ne. 0 ) then\n\
\t  norm = 1/(rtree + rloop)\n" <>
    vxCode["\
\t  do i = 1, hsimd\n\
\t    hseltest_v(HelAll(i)) = hseltest_v(HelAll(i)) +\n\
     &        abs(norm*(res(HelAll(1),i) + res(HelAll(2),i)))\n",
    (* vxElse *) "\
\t  do i = 1, hmax\n\
\t    hseltest(i) = hseltest(i) +\n\
     &        abs(norm*(res(1,i) + res(2,i)))\n"
    ] <> "\
\t  enddo\n\
\t  hseli = hseli + 1\n\
\t  if( hseli .eq. hseln ) then\n\
\t    hselmin = 0\n\
\t    do i = 0, hmax - 1\n\
\t      hselmin = max(hselmin, hseltest(i))\n\
\t    enddo\n\
\t    hselmin = hselmin*hseleps\n\
\t    do i = 0, hmax - 1\n\
\t      if( hseltest(i) .ge. hselmin ) cycle" <>
  Thread[{"\n\t      ", hfun}] <> "\n\
\t      print 1, ' neglecting ', Hel0\n\
1\t      format(A, LEGS I3)\n\
\t    enddo\n\
\t  endif\n\
\tendif\n\
* " <> prefix <> "END HSEL_SET\n\n\
#ifdef CHECK" <>
    ({"\n\tprint *, '", #, " =', ", #}&)/@
      (ToString[#1]&)@@@ invs <> "\n\
\tprint *, 'tree =', rtree\n\
\tprint *, 'loop =', rloop\n\
\tstop\n\
#endif\n\n\
\thlast = helicities\n\n" <>
    sEnd <> "\
\tend\n\n"];
]


writeSquaredMEC[hh_, part_, hfun_, hmax_] :=
Block[ {jtree, jloop, ffcode, ffind = {}, Global`h},
  jtree = LoopElem[Ctree, "j"];
  jloop = LoopElem[Cloop, "j"];
  ffcode = "\n\
  Zero(ampsq);" <> ({ "\n\
// ", prefix, "BEGIN FF_TREE",
    #3, "\n  ", ToCode[#1 -> 0], ";", #4, "\n\n",
    CallDecl[ToDoLoops[ffmods[[1]], Indices]], "\
// ", prefix, "END FF_TREE\n\n\
// ", prefix, "BEGIN M2_TREE",
    Cond[ TrueQ[$TreeSquare],
      sumup[jtree, LoopElem[Ctree, "i"]] ], "\n\
// ", prefix, "END M2_TREE"
  }&)@@ jtree <> "\n\
  res[0] = ampsq;\n\n\
  Zero(ampsq);" <> ({"\n\
  TEST(flags, BIT_LOOP)\n\
// ", prefix, "BEGIN FF_LOOP",
    #3, "\n  ", ToCode[#1 -> 0], ";", #4, "\n\n",
    CallDecl[ToDoLoops[ffmods[[2]], Indices]], "\
// ", prefix, "END FF_LOOP\n\n\
// ", prefix, "BEGIN M2_LOOP",
    Cond[ jtree =!= "",
      sumup[jtree, LoopElem[Cloop, "i"]],
      "\n  ampsq += ampsq;" ],
    Cond[ jtree === "" || TrueQ[$LoopSquare],
      sumup[jloop, LoopElem[Cloop, "i"]] ], "\n\
// ",  prefix, "END M2_LOOP\n\
  ENDTEST(flags, BIT_LOOP)"
  }&)@@ jloop <> "\n\
  res[1] = ampsq;\n\n";

  WriteString[hh, "\
//#define CHECK\n\n" <>
    Hdr["assembly of squared matrix element"] <> "\
#include <math.h>\n\
#include \"" <> prefix <> "specs.h\"\n" <>
    vPre <> fPost <> "\
#if NOUNDERSCORE\n\
#define " <> $SymbolPrefix <> "SquaredME " <> $SymbolPrefix <> "squaredme\n\
#else\n\
#define " <> $SymbolPrefix <> "SquaredME " <> $SymbolPrefix <> "squaredme_\n\
#endif\n\n" <>
    ({"void ", $SymbolPrefix, #, "(void);\n"}&)/@
      Flatten[{abbrmods, ffmods}] <>
    "\n" <>
    varSeq["struct formfactors ", $SymbolPrefix, "formfactors = {\n",
      { "  .", ToString[#1],
        {"[0 ... ", ToString[# - 1], "]"}&/@ Reverse[#2],
        " = ", #3 }&@@@ mats,
      ",\n", " };\n\n"] <> "\
#define RESFMT \"%25.12f\"\n\
#define HELFMT \"" <> Table["%3d", {legs}] <> "\"\n\
#define HELOUT(h,v) " <>
      StringDrop[StringJoin[Array[{"h(", ToString[#], ")v,"}&, legs]], -1] <> "\n\n\
/**********************************************************************/\n\n\
static void " <> $SymbolPrefix <> "SquaredMEHel(" <>
      vxCode["HelArg(cinteger vmax, ResType *res)",
             "RealType *res"] <> ", cinteger *pflags) {\n\n" <>
    vPre <> "\n\
  cinteger flags = *pflags;\n\
  " <> helType <> " amp;\n\
  " <> resType <> " ampsq;\n" <>
    vxCode["  SIMD_ONLY(integer v;)\n"] <>
    VarDecl[Union[Flatten[ffind]], "integer"] <> "\n\
// " <> prefix <> "BEGIN ABBR_HEL\n" <>
    Invoke@@ abbrmods[[3]] <> "\
// " <> prefix <> "END ABBR_HEL\n" <>
    ffcode <>
    $DebugPre[1, 3] <>
    vxCode["\
  SIMD_ONLY(for( v = 0; v < vmax; ++v ))\n\
!printf(\" sqme(\" HELFMT \") = \" RESFMT RESFMT \"\\n\",\n\
    HELOUT(Hel,HelInd(v)), res[0] HelInd(v), res[1] HelInd(v));\n",
    (* vxElse *) "\
!printf(\" sqme(\" HELFMT \") = \" RESFMT RESFMT \"\\n\",\n\
    HELOUT(Hel,), res[0], res[1]);\n"
    ] <>
    $DebugPost[1] <> "\
}\n\n\
/**********************************************************************/\n\n\
void " <> $SymbolPrefix <>
    "SquaredME(RealType *result, integer8 *phelicities, cinteger *pflags) {\n\n" <>
    vPre <> "\n\
// " <> prefix <> "BEGIN VAR_DECL\n\
  cinteger flags = *pflags;\n\
  integer8 helicities, hbits;\n\
  static integer8 hlast;\n\
  static integer seq[2];
  integer i, h;\n\
  enum { hmax = " <> ToString[hmax] <> " };\n" <>
    vxCode["\
  SIMD_ONLY(integer v;)\n\
  enum { hsimd = SIMD_CEIL(hmax) };\n\
  ResType res[hsimd][2];\n",
    (* vxElse *) "\
  RealType res[hmax][2];\n"
    ] <> "\
  RealType rtree, rloop;\n\
// " <> prefix <> "END VAR_DECL\n\n\
// " <> prefix <> "BEGIN HSEL_DECL\n\
  RealType norm;\n" <>
    vxCode["\
  static union {\n\
    RealType s[hmax];\n\
    ResType v[hsimd];\n\
  } hseltest;\n",
    (* vxElse *) "\
  RealType hseltest[hmax];\n"
    ] <> "\
  static RealType hselmin;\n\
  static integer hseli;\n\
// " <> prefix <> "END HSEL_DECL\n\n" <>
    sPost <> "\
// " <> prefix <> "BEGIN SETMASS\n\
  TEST(flags, BIT_SETMASS)" <>
    MapIndexed[{"\n  result[", ToString[#2[[1]] - 1], "] = ",
      ToCode[ #1[[3]] ], ";"}&, part] <> "\n\
  return;\n\
  ENDTEST(flags, BIT_SETMASS)\n\
// " <> prefix <> "END SETMASS\n\n"];

  WriteExpr[hh, {"\
// " <> prefix <> "BEGIN INVARIANTS\n",
    invs, "\
// " <> prefix <> "END INVARIANTS\n\n"}, Newline -> ""];

  WriteString[hh, "\
  CHK_INI(seq);\n\n\
  TEST(flags, BIT_RESET)\n\
  hlast = 0;\n\
  ++seq[0];\n\
  INI_S();\n\
// " <> prefix <> "BEGIN ABBR_S\n" <>
    Invoke@@ abbrmods[[1]] <> "\
// " <> prefix <> "END ABBR_S\n\
  ENDTEST(flags, BIT_RESET)\n\n\
  ++seq[1];\n\
  INI_A();\n\
// " <> prefix <> "BEGIN ABBR_ANGLE\n" <>
    Invoke@@ abbrmods[[2]] <> "\
// " <> prefix <> "END ABBR_ANGLE\n\n\
// " <> prefix <> "BEGIN HEL_LOOP\n\
  helicities = *phelicities &= Generic(ARG_Phel,JOIN_HEL);\n\n" <>
    vxCode["\
  SIMD_ONLY(v = 1;)\n"
    ] <> "\
  h = 0;\n\
  for( i = 0; i < hmax; ++i ) {\n\
// " <> prefix <> "BEGIN HSEL_IF\n\
    if( hseltest" <> vxCode[".s"] <> "[i] < hselmin ) continue;\n\
// " <> prefix <> "END HSEL_IF\n\n\
#define hval(s,m,a,b) a*(i/s % m)-b\n" <>
    dfCode["\
#define hval_df(i) HelBit(helicities,i,1) - HelBit(helicities,i,-1)\n"
    ] <>
    Thread[{"\n    ", hfun, ";"}] <>
    If[ hmax === 1, "\n", "\n\
    hbits = Generic(ARG_Bhel,JOIN_HEL);\n\
    if( (hbits & helicities) != hbits ) continue;\n"
    ] <>
    $DebugPre[1, 1] <> "\
!  if( helicities != hlast ) printf(\" helicities: \" HELFMT \"\\n\", HELOUT(Hel0,));\n" <>
    $DebugPost[1] <>
    vxCode["\
    SIMD_ONLY(VecCopy(v, LEGS);)\n\
    SIMD_ONLY(v = v % SIMD + 1;)\n\
    SIMD_ONLY(if( v == 1 ))\n\
    " <> $SymbolPrefix <> "SquaredMEHel(HelArg(SIMD, res[h++]), pflags);\n\
  }\n\n\
  SIMD_ONLY(if( v != 1 ) " <> $SymbolPrefix <> "SquaredMEHel(HelArg(v, res[h++]), pflags);)\n\
  DEINI();\n\
  SIMD_ONLY(if( v != 1 ) res[h-1][0][1] = res[h-1][1][1] = 0;)\n",
    (* vxElse *) "\
    " <> $SymbolPrefix <> "SquaredMEHel(res[h++], pflags);\n\
  }\n\n\
  DEINI();\n"
    ] <> "\
// " <> prefix <> "END HEL_LOOP\n\n\
// " <> prefix <> "BEGIN RESULT\n\
  rtree = rloop = 0;\n\
  for( i = 0; i < h; ++i ) {\n" <>
    vxCode["\
    rtree += HelSum(res[i][0]);\n\
    rloop += HelSum(res[i][1]);\n",
    (* vxElse *) "\
    rtree += res[i][0];\n\
    rloop += res[i][1];\n"
    ] <> "\
  }\n" <>
  dfCode["\
  hbits = helicities & Generic(ARG_Ferm,JOIN_HEL);\n\
  i = 1 << BitCount(hbits & (hbits >> 2));\n\
  rtree = rtree*i;\n\
  rloop = rloop*i;\n"
  ] <> "\
  result[0] = rtree;\n\
  TEST(flags, BIT_LOOP)\n\
  result[1] = rloop;\n\
  ENDTEST(flags, BIT_LOOP)\n\n" <>
    $DebugPre[1, 2] <> "\
!printf(\" \" PROCNAME \" =\" RESFMT RESFMT \"\\n\", rtree, rloop);\n" <>
    $DebugPost[1] <> "\
// " <> prefix <> "END RESULT\n\n\
// " <> prefix <> "BEGIN HSEL_SET\n\
  if( helicities != hlast ) {\n\
    Zero(hseltest" <> vxCode[".v"] <> ");\n\
    hselmin = 0;\n\
    hseli = 0;\n\
  }\n\
  if( hseli < hseln && rtree + rloop != 0 ) {\n\
    norm = 1/(rtree + rloop);\n" <>
    vxCode["\
    for( i = 0; i < hsimd; ++i )\n\
      hseltest.v[i] += AbsRxRes(norm, res[i][0] + res[i][1]);\n",
    (* vxElse *) "\
    for( i = 0; i < hmax; ++i )\n\
      hseltest[i] += fabs(norm*(res[i][0] + res[i][1]));\n"
    ] <> "\
    if( ++hseli == hseln ) {\n\
      hselmin = 0;\n\
      for( i = 0; i < hmax; ++i )\n\
        hselmin = fmax(hselmin, hseltest" <> vxCode[".s"] <> "[i]);\n\
      hselmin *= hseleps;\n\
      for( i = 0; i < hmax; ++i ) {\n\
        if( hseltest" <> vxCode[".s"] <> "[i] >= hselmin ) continue;" <>
    Thread[{"\n        ", hfun, ";"}] <> "\n\
!      printf(\" neglecting \" HELFMT \"\\n\", HELOUT(Hel0,));\n\
      }\n\
    }\n\
  }\n\
// " <> prefix <> "END HSEL_SET\n\n\
#ifdef CHECK\n" <>
    ({"!printf(\"", #, " = %g\\n\", ", #, ");\n"}&)/@
      (ToString[#1]&)@@@ invs <> "\
!printf(\"tree = %g\\n\", rtree);\n\
!printf(\"loop = %g\\n\", rloop);\n\
!exit(1);\n\
#endif\n\n\
  hlast = helicities;\n\n" <>
    sEnd <> "\
}\n"];
]


(* renormalization constants *)

rcpatt[_[_[_[h_Symbol[___]]], _]] := HoldPattern[_h]

rcpatt[_[_[_[h_Symbol]], _]] := HoldPattern[h]

RCPattern[h_, other___] :=
  Level[{{other}, Union[rcpatt/@ DownValues[h]]}, {2}, Alt]


Attributes[WithOpt] = {HoldFirst}

WithOpt[foo_] := (Needs["FeynArts`"]; foo)

WithOpt[foo_, opt__] := (
  Needs["FeynArts`"];
  (Options[CreateTopologies] = #1;
   Options[InsertFields] = #2;
   Options[CreateFeynAmp] = #3;
   #4)&[
    Options[CreateTopologies],
    Options[InsertFields],
    Options[CreateFeynAmp],
    SetOptions[CreateTopologies,
      ExcludeTopologies -> Internal,
      FilterOpt[CreateTopologies, opt]];
    SetOptions[InsertFields,
      FilterOpt[InsertFields, opt]];
    SetOptions[CreateFeynAmp,
      FilterOpt[CreateFeynAmp, opt]];
    foo
  ] )


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
  a /. Cases[{cond}, i_ == j_ :> (IndexDelta[i, j] -> 1), Infinity],
  b /. Cases[{cond}, i_ == j_ :> (IndexDelta[i, j] -> 0)] ]


Attributes[FermionRC] = {HoldRest}
	(* HoldRest postpones TheMass until model is initialized *)

FermionRC[se_, m_, a_, b_] :=
  m (a SEPart[LVectorCoeff, se] + b SEPart[RVectorCoeff, se]) +
     b SEPart[LScalarCoeff, se] + a SEPart[RScalarCoeff, se]

BosonRC[se_] := SEPart[Identity, se]


RCArg[(f:P$Field)[k_]] := f :> k

RCArg[f:P$Field] := f :> TheMass[f]

RCArg[other_] := other


MassRC[f_, f_, opt___Rule] := MassRC[f, opt]

MassRC[f__, opt___Rule] := (RCArg/@ massRC[f])[opt]

massRC[f_F :> m_][opt___] :=
  FermionRC[ReTilde[SelfEnergy[f, m, opt]], m, 1/2, 1/2]

massRC[f_ :> m_][opt___] := BosonRC[ReTilde[SelfEnergy[f, m, opt]]]

massRC[f1_ :> m1_, f2_ :> m2_][opt___] :=
  1/2 (BosonRC[ReTilde[SelfEnergy[f2 -> f1, m1, opt]]] +
       BosonRC[ReTilde[SelfEnergy[f2 -> f1, m2, opt]]])


FieldRC[f_, f_, opt___Rule] := FieldRC[f, opt]

FieldRC[f__, opt___Rule] := (RCArg/@ fieldRC[f])[opt]

fieldRC[f_F :> m_][opt___] := (
  -{SEPart[LVectorCoeff, #1], SEPart[RVectorCoeff, #1]} - FermionRC[#2, m, m, m]
)&[ ReTilde[SelfEnergy[f, m, opt]], ReTilde[DSelfEnergy[f, m, opt]] ]

fieldRC[f_ :> m_][opt___] := -BosonRC[ReTilde[DSelfEnergy[f, m, opt]]]

fieldRC[f1_, f2_, c_:0][opt___] := OnlyIf[
  And@@ MapThread[Equal, Flatten[{##}]&@@@ {First[f1], First[f2]}],
  fieldRC[f1][opt],
  fieldRC2[f1, f2, c][opt]
]

fieldRC2[f1_F :> m1_, f2_F :> m2_, c_][opt___] := (
  2/(m1^2 - m2^2) (FermionRC[#, m2, {m2, m1}, {m1, m2}] - c)
)&[ ReTilde[SelfEnergy[f2 -> f1, m2, opt]] ]

fieldRC2[f1_ :> m1_, f2_ :> m2_, c_][opt___] := (
  2/(m1^2 - m2^2) (BosonRC[#] - c)
)&[ ReTilde[SelfEnergy[f2 -> f1, m2, opt]] ]


TadpoleRC[f_, opt___Rule] :=
  -BosonRC[ReTilde[SelfEnergy[f -> {}, Indeterminate, opt]]]


WidthRC[f_, opt___Rule] := widthRC[RCArg[f]][opt]

widthRC[f_F :> m_][opt___] :=
  FermionRC[ImTilde[SelfEnergy[f, m, opt]], m, 1, 1]

widthRC[f_ :> m_][opt___] :=
  BosonRC[ImTilde[SelfEnergy[f, m, opt]]]/m


TreeCoupling[proc_Rule, opt___Rule] :=
  WithOpt[CalcProcess[0, proc], opt]


VertexFunc[proc_Rule, opt___Rule] :=
  WithOpt[CalcProcess[1, proc], opt]


Attributes[SEHook] = {HoldAll}

SEHook[se_, amp_, rul_] := (
  FCPrint[2, HoldForm[se]];
  amp /. rul
)


Attributes[SelfEnergy] = {HoldRest}

se:SelfEnergy[proc_Rule, m_, opt___Rule] := SEHook[se,
  WithOpt[CalcSelfEnergy[proc], opt], K2 -> m^2]

SelfEnergy[f:Except[_Rule], opt___Rule] :=
  SelfEnergy[f -> f, TheMass[f], opt]

SelfEnergy[f:Except[_Rule], m__] := SelfEnergy[f -> f, m]


Attributes[DSelfEnergy] = {HoldRest}

se:DSelfEnergy[proc_Rule, m_, opt___Rule] := SEHook[se,
  WithOpt[D[CalcSelfEnergy[proc], K2], opt], K2 -> m^2]

DSelfEnergy[f:Except[_Rule], opt___Rule] :=
  DSelfEnergy[f -> f, TheMass[f], opt]

DSelfEnergy[f:Except[_Rule], m__] := DSelfEnergy[f -> f, m]


RCSub = Simplify

RCInt = Simplify


ClearSE[] := Clear[RCCache]

CalcSelfEnergy[proc_] := CalcSelfEnergy[RCCache @
  Hash[{proc, Options[CreateTopologies], Options[InsertFields]}], proc]

CalcSelfEnergy[c_RCCache, proc_] := c = CalcProcess[1, proc] /. {
  Pair[_k, _k] -> K2,
  Pair[_e | _ec, _k] -> If[ MatchQ[proc, _V -> _V],
	(* default: take only the transverse part of vector-boson SEs *)
    If[TrueQ[$LongitudinalSE], I Sqrt[K2], 0],
    1 ],
  Pair[_e, _ec] -> -1,
  SUNT[_, _] -> 1,
  SUNT[_, _, 0, 0] -> 1/2 }

CalcSelfEnergy[expr_, _] := expr


CalcProcess[loop_, proc_] :=
Block[ {Neglect, FormSub = RCSub},
  FCPrint[2, "calculating ", loop, "-loop ", proc];
  ClearProcess[];
  $RCTop = CreateTopologiesHook[loop, Length[Flatten[{#}]]&/@ proc];
  $RCIns = InsertFieldsHook[$RCTop, proc];
  PaintSE[$RCIns];
  $RCAmp = CreateFeynAmpHook[$RCIns, Truncated -> !FreeQ[proc, F]];
  $RCRes = Unabbr @ CalcFeynAmp[$RCAmp,
    OnShell -> False, Transverse -> False,
    FermionChains -> Chiral, FermionOrder -> None,
    OPP -> False, FileTag -> "rc"];
  PutSE[{$RCTop, $RCIns, $RCAmp, $RCRes}];
  Plus@@ $RCRes /. Mat -> Identity
]


CreateTopologiesHook[args__] := CreateTopologies[args]

InsertFieldsHook[args__] := InsertFields[args]

CreateFeynAmpHook[args__] := CreateFeynAmp[args]


PaintSE[ins_] := PaintSE[ins, Sequence@@ $PaintSE]

PaintSE[ins_, True, ___] := (
  Paint[ins];
  ins
)

PaintSE[ins_, prefix_String, suffix___String] :=
Block[ {file = ChkExist[prefix <> ProcName[ins] <> suffix <> ".ps"]},
  FCPrint[2, "diagrams in ", file];
  Paint[ins, DisplayFunction -> (Export[file, #]&)];
  ins
]

PaintSE[ins_, __] := ins


PutSE[expr_] := PutSE[expr, Sequence@@ $PutSE]

PutSE[{top_, ins_, amp_, res_}, prefix_String, suffix___String] :=
Block[ {base = ChkExist[prefix <> ProcName[ins] <> suffix]},
  FCPrint[2, "results in ", base];
  Put[top, base <> ".top"];
  Put[ins, base <> ".ins"];
  Put[amp, base <> ".amp"];
  Put[res, base <> ".res"];
  res
]


ProcName[(TopologyList | FeynAmpList)[info__][___]] :=
  ProcName@@ ({ Process,
    Model, GenericModel,
    ExcludeParticles, ExcludeFieldPoints, LastSelections
  } /. {info})

ProcName[Amp[proc_][___]] := ProcName[proc]

ProcName[proc_, opt__] :=
  ProcName[proc] <> "_" <>
    FromCharacterCode[IntegerDigits[Hash[{opt}], 26] + 97]

ProcName[proc_] := StringJoin@@ Map[pname, proc, {2}]


pname[{f_, ___}] := pname[f]

pname[-f_[a___]] := pname[f["b"][a]]

pname[f_] := StringJoin[ToString/@ DeleteCases[
  Level[f /. _Index -> "i" /.
      i_List :> Replace[i, _Symbol -> "i", {1, Infinity}],
    {-1}, Heads -> True],
  s_Symbol /; Context[s] === "System`" ]]


$RenConst = {RenConst (*, MassShift*)}


Attributes[FindRC] = {Listable}

FindRC[expr_, h_] :=
Block[ {test = expr, rcp = RCPattern[h], rcs = {}, new},
  While[
    Apply[(orbit[#1] = Range[##2])&,
      { Cases[test, SumOver[i_, r_, ___] :> {i, r}, Infinity],
        Cases[test, IndexSum[_, r___] :> r, Infinity] }, {2}];
    Length[new = Complement[
      Flatten[Cases[test, rc:rcp :> Distribute[orbit/@ rc, List], Infinity]],
      rcs, SameTest -> (isym/@ #1 === isym/@ #2 &) ]] =!= 0,
    rcs = Flatten[{new, rcs}];
    test = h/@ new ];
  RenConstList[h]@@ Cases[rcs, rcp]
]


General::nodefs = "Warning: no definitions for `` found."

Off[MassShift::nodefs]

FindRenConst::nodef =
"Warning: `` might be renormalization constants but have no definition."

FindRenConst[expr_] := FindRenConst[expr, $RenConst]

FindRenConst[expr_, h_] :=
Block[ {test, orbit, isym, rcs, x, SelfEnergy, DSelfEnergy},
  test = x[expr] /. Rule -> (#2 &) //. Dispatch[Subexpr[]];

  Needs["FeynArts`"];
  If[ $Model === "",
    InitializeModel[
      Model /. Cases[test, HoldPattern[Model -> _], Infinity] /.
        Options[InsertFields],
      GenericModel -> (GenericModel /. Options[InsertFields]),
      Reinitialize -> False ] ];

  Scan[(isym[#] = x)&, LoopInd[]];
  isym[other_] := other;

  orbit[other_] := other;

  rcs = FindRC[test, h];
  Message[#::nodefs, #]&/@ Cases[rcs, RenConstList[hx_][] :> hx, {-2}];

  x = Select[ Names["Global`d*"],
    (FreeQ[rcs, #] && !FreeQ[test, #]&)[
      ToExpression[#, InputForm, HoldPattern] ]& ];
  If[ Length[x] =!= 0, Message[FindRenConst::nodef, x] ];

  rcs
]


Attributes[RenConstHook] = {HoldRest}

RenConstHook[rc_, expr_] := (FCPrint[1, rc]; rc -> expr)


Attributes[CalcRC] = {Listable}

CalcRC[rcs:RenConstList[h_][___], opt___] := RenConstHook[#, WithOpt[
  Expand[h[#], IndexSum] /.
    IndexSum[x_, {i_, r__}]^n_. :> Product[
      ((x /. i -> #) SumOver[#, r])&[ NewSymbol[ToString[i]] ], {n} ],
  Options[kind[#]], opt
]]&/@ rcs


CalcRenConst[expr_, h___Symbol, opt___Rule] :=
  CalcRenConst[FindRenConst[expr, h], opt] /; FreeQ[expr, RenConstList]

CalcRenConst[expr_, h___Symbol, opt___Rule] :=
  CalcRC[expr, opt] /. Plus -> IntCollect


IntCollect[p__] := Plus[p] /; FreeQ[{p}, Re]

(*IntCollect[p__] := Collect[Plus[p], _Re, RCInt]*)

IntCollect[p__] :=
Block[ {RCInt},
  Replace[
    Collect[ Plus[p],
      First/@ DeleteCases[Split @ Sort @ Cases[Plus[p], _Re, Infinity], {_}],
      RCInt ],
    RCInt[x_] :> x, {1} ]
]


RCMod[rcs_, mod_] :=
Block[ {hh},
  FCPrint[3, "  writing ", mod];
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["renormalization constants"] <>
    fPre <> fPost <>
    SubroutineDecl[mod] <>
    sPre <> "\n" <>
    VarDecl[Union[Cases[rcs, SumOver[i_, _] :> i, Infinity]], "integer"]];
  WriteExpr[hh, {sPost, rcs, sEnd},
    Optimize -> True,
    DebugLines -> $DebugRC, DebugLabel -> mod];
  WriteString[hh, SubroutineEnd[]];
  Close[hh];
  mod
]

$DebugRC = 1


RCAll[mod_, mods_] :=
Block[ {hh},
  hh = OpenCode[ModName[mod]];
  WriteString[hh,
    Hdr["RC invocations"] <>
    fPre <> "\n\n" <>
    SubroutineDecl[mod] <>
    ({"\n\tcall ", $SymbolPrefix, #}&/@ mods) <>
    "\n" <> SubroutineEnd[]];
  Close[hh];
  {mod, mods}
]


Attributes[WriteRC] = {Listable}

WriteRC[RenConstList[h_][rcs___]] :=
Block[ {fname = fname /. Automatic :> ToString[h], rcmods},
  Block[ {varDecl, hh},
    varDecl = varDeclF;
    hh = OpenFortran[ModName[fname, ".h.F"]];
    WriteString[hh, #1 <> IfDecl[#2, {}, VarDecl[#3]]];
    Close[hh];
    varDecl = varDeclC;
    hh = OpenC[ModName[fname, ".h.c"]];
    WriteString[hh, #1 <> IfDecl[#2, VarDecl[#3], {}]];
    Close[hh]
  ]&[ Hdr["RC declarations"],
    fname <> "_h",
    Common["var" <> fname][MaxDims[Map[Dim, First/@ {rcs}, {2}]],
      "ComplexType"]
  ];

  rcmods = FileSplit[ToDoLoops[OnePassOrder[{rcs} /. xrules]],
    fname, RCMod, RCAll];

  {"\
OBJS :=", MakefileObjs[rcmods], "\n\n\
OBJS_RC += $(OBJS)\n\n\
$(LIB)($(OBJS)): ", $MakeDeps[[2]], MakefileName[fname, ".h" <> $CodeExt], "\n\n\n"}
]


Options[WriteRenConst] = {
  Folder -> "renconst",
  ExtraRules -> ExtraRules (* i.e. from WriteSquaredME *),
  FilePrefix -> FilePrefix,
  SymbolPrefix -> SymbolPrefix,
  FileName -> Automatic,
  FileHeader -> FileHeader,
  FileIncludes -> FileIncludes,
  SubroutineIncludes -> SubroutineIncludes }

WriteRenConst[expr_, h___Symbol, dir_, opt___Rule] :=
  WriteRenConst[CalcRenConst[expr, h], dir, opt] /;
    FreeQ[expr, _RenConstList[___Rule]]

WriteRenConst[rcs_, dir_, opt___Rule] :=
Block[ {folder, xrules, prefix, $SymbolPrefix, fname, fincl, sincl,
fPre, fPost, sPre, sPost, sEnd, header, ModName, Hdr, file, hh, mkcmds},

  {folder, xrules, prefix, $SymbolPrefix, fname, header, fincl, sincl} =
    ParseOpt[WriteRenConst, opt] //. Options[WriteSquaredME];

  {fPre, fPost} = fDecl[fincl];
  fPre = "#define RENCONST\n\n" <> fPre;
  {sPre, sPost, sEnd} = sDecl[sincl];

  FCPrint[2, "writing renconst modules"];

  DefModName[MkDir[dir, folder]];
  mkcmds = WriteRC[rcs];

  FCPrint[2, "writing makefile"];

  hh = OpenWrite[ModName["makefile", ""]];
  WriteString[hh, "\
LIBS += $(LIB)\n\
VPATH := $(DIR):$(VPATH)\n\
OBJS_RC :=\n\n" <>
    mkcmds <> "\
$(LIB): $(LIB)($(OBJS_RC))\n\n"];
  Close[hh];

  Cases[DownValues[ModName], _[_, s_String] :> s]
]


(* low-level output functions *)

SetLanguage::badlang = "Unknown language ``."

SetLanguage["C"] := (
  toCode[x_] := Block[{vec = Identity}, ToString[x, CForm]];
  varDecl = varDeclC;
  SubroutineDecl = subroutineDeclC;
  SubroutineEnd[] = "}\n\n";
  callDecl = callDeclC;
  ifDecl = ifDeclC;
  helDim = helDimC;
  helType = "HelType";
  resType = "ResType";
  vxCode[s_, ___] := s;
  vx = vec;
  vxrul = {Times -> vtimes, Plus -> vplus};
  writeSquaredME = writeSquaredMEC;
  writeNum = writeNumC;
  OpenCode = OpenC;
  $MakeDeps = {"C-squaredme.d", "C-renconst.d"};
  $CodeExt = ".c";
  $CodeIndent = "  ";
  $CodeBoln = "\n  ";
  $CodeEoln = ";";
  $CodeIf = "  if( ";
  $CodeThen = " ) {\n";
  $CodeElseif = "  }\n  else if( ";
  $CodeElse = "  }\n  else {\n";
  $CodeEndIf = "  }\n";
  $CodeCall = "  ";
  $Code = "C";
)

SetLanguage["Fortran"] = (
  toCode[x_] := ToString[x, FortranForm];
  varDecl = varDeclF;
  SubroutineDecl = subroutineDeclF;
  SubroutineEnd[] = "\tend\n\n";
  callDecl = callDeclF;
  ifDecl = ifDeclF;
  helDim = helDimF;
  helType = "HelType";
  resType = "ResType";
  ampvar = "amp";
  vxCode[s_, ___] := s;
  vx = Identity;
  vxrul = {};
  writeSquaredME = writeSquaredMEF;
  writeNum = writeNumF;
  OpenCode = OpenFortran;
  $MakeDeps = {"F-squaredme.d", "F-renconst.d"};
  $CodeExt = ".F";
  $CodeIndent = "";
  $CodeBoln = "\n\t";
  $CodeEoln = Sequence[];
  $CodeIf = "if( ";
  $CodeThen = " ) then\n";
  $CodeElseIf = "else if( ";
  $CodeElse = "else\n";
  $CodeEndIf = "endif\n";
  $CodeCall = "call ";
  $Code = "Fortran" )

SetLanguage[lang_, "novec"] := (
  SetLanguage[lang];
  helDim = #1 &;
  helType = "ComplexType";
  resType = "RealType";
  vxCode[_, s___] := s;
  vx = Identity )

SetLanguage[lang_] := (
  Message[CodeLanguage::badlang, lang];
  OpenCode = Abort;
  $Failed )

OpenC[file_String, pipe___String, opt___Rule] :=
  OpenWrite[toc <> pipe <> " > " <> Escape[file],
    FormatType -> CForm, opt, PageWidth -> 67]

OpenFortran[file_String, pipe___String, opt___Rule] :=
  OpenWrite[tofortran <> pipe <> " > " <> Escape[file],
    FormatType -> FortranForm, opt, PageWidth -> 67]

tofortran = "!" <> Escape[ToFileName[$FormCalcBin, "ToFortran"]]

toc = "!" <> Escape[ToFileName[$FormCalcBin, "ToC"]]

Format[AddrOf[x_], CForm] := SequenceForm["&", x]

Format[AddrOf[x_], FortranForm] := x


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

Coalesce[(ru:Rule | RuleAdd)[v_, x_. p_Plus], r___] :=
  Level[
    { ReplacePart[
        RuleAdd[v, x Plus@@ #]&/@ Flatten[Coalesce@@ p],
        ru, {1, 0} ],
      {r} }, {2}, Coalesce ] /; AtomQ[x] && LeafCount[p] > size

Coalesce[DoLoop[rul_, ind__], r___] := Level[
  {DoLoop[{##}, ind]&@@@ Flatten[Coalesce@@ rul], {r}},
  {2}, Coalesce ] /; LeafCount[rul] > size

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
  {writemod[expr, mod /. Dot -> (#1&)]} /; LeafCount[expr] < $FileSize

FileSplit[expr_List, mod_, writemod_, writeall_:(#2&)] :=
Block[ {size = $FileSize, m = mod /. Dot -> StringJoin},
  writeall[mod, MapIndexed[
    writemod[List@@ #1, NumberMod[m, ToString@@ #2]]&,
    Flatten[Coalesce@@ expr] ]]
]

FileSplit[CodeExpr[vars__, expr_List], mod_,
  writemod_, writeall_:(#2&)] :=
Block[ {size = $FileSize, m = mod /. Dot -> StringJoin},
  writeall[mod, MapIndexed[
    writemod[CodeExpr[vars, List@@ #1], NumberMod[m, ToString@@ #2]]&,
    Flatten[Coalesce@@ expr] ]]
]

FileSplit[other_, r__] := FileSplit[{other}, r]


isdup[expr_, {i_, j__}] := {expr -> dup[expr], Min[First/@ {i, j}]}

_isdup = {}

RemoveDups[expr_CodeExpr] := expr

RemoveDups[expr_] :=
  Fold[RemoveDups, #, -Range[3, Depth[#] - 1]]& @ Flatten[{expr}]

RemoveDups[expr_, lev_] :=
Block[ {tmps, new},
  tmps = isdup[#, Position[expr, #]]&/@
    Union[Cases[expr /. {_DoLoop -> 1, _IndexIf -> 1},
      p_Plus /; LeafCount[N[p]] > minleaf, {lev}]];
  new = Block[{Plus}, Apply[Set, tmps, {2}]; expr];
  new = #1 /. Flatten[#2] & @@
    Reap[new /. r:(_dup -> _dup) :> (Sow[r]; {})];
  Fold[insdef, new, Reverse[tmps]] //Flatten
]

insdef[expr_, {var_ -> tmp_, pos_}] :=
  MapAt[{tmp -> var, #}&, expr, pos]

insdef[expr_, _] := expr


Attributes[TmpList] = {HoldFirst}

TmpList[expr_] := Reverse[Reap[expr]]

ToTmp[expr_] := (Sow[# -> expr]; #)& @ tmp[expr]

$TmpPrefix = "tmp"

$DupPrefix = "dup"

tmp[expr_] := NewSymbol[$TmpPrefix[Head[expr]], 0]

dup[expr_] := NewSymbol[$DupPrefix[Head[expr]], 0]


Attributes[SplitExpr] = {Listable}

SplitExpr[(ru:Rule | RuleAdd)[var_, expr_]] :=
  BlockSplit @ TmpList[ ru[var, psplit[expr]] ]

SplitExpr[other_] := other


psplit[x_?AtomQ] := x

psplit[x_?AtomQ p_Plus] := x psplit[p]

psplit[p_[t__]] := p@@ ({t} /.
  { h_HoldForm :> h,
    Plus -> (If[LeafCount[#] > $BlockSize, ToTmp[#], #]&[Plus[##]]&) })


Attributes[RhsMap] = {Listable}

RhsMap[foo_, (ru:Rule | RuleAdd)[var_, expr_]] := ru[var, foo[expr]]

RhsMap[_, other_] := other


Options[PrepareExpr] = {
  Optimize -> False,
  Expensive -> {},
  MinLeafCount -> 10,
  DebugLines -> 0,
  DebugLabel -> True,
  MakeTmp -> Identity,
  Declarations -> (Rule | RuleAdd)[var_, _] :> var,
  FinalTouch -> Identity,
  ResetNumbering -> True }

PrepareExpr[expr_, opt___Rule] :=
Block[ {optim, expen, minleaf, debug, debtag, mktmp, decl, final, reset,
process, doloop, new, vars, tmps},
  {optim, expen, minleaf, debug, debtag, mktmp, decl, final, reset} =
    ParseOpt[PrepareExpr, opt];
  process = RhsMap[final, Flatten[SplitExpr @ Prep[#]]] &;
  If[ TrueQ[optim], process = process /. p_Prep :> RemoveDups[p] ];
  If[ reset, DownValues[SymbolNumber] = Select[DownValues[SymbolNumber],
    FreeQ[#, $TmpPrefix | $DupPrefix]&] ];
  doloop = Hoist[Alt[expen]];
  new = unpatt[expr];
  vars = Cases[new, decl, Infinity];
  If[ debug > 0, new = addDebug[new] ];
  new = process[mktmp[new]];
  new = new /. CodeExpr[_, t_, x_] :> (vars = Complement[vars, t]; x);
  tmps = Cases[new, decl, Infinity];
  If[ debug < 0, new = addDebug[new] ];
  CodeExpr[MaxDims[vars], MaxDims[Complement[tmps, vars]], Flatten[new]]
]


Attributes[addDebug] = {Listable}

addDebug[DoLoop[expr_, ind__]] := DoLoop[
  { DebLine[1, debug, DebTag[debtag][(#1&@@@ {ind}) -> DoLoop]],
    addDebug[expr] },
  ind ]

addDebug[i_IndexIf] := MapIf[addDebug, i]

addDebug[NoDebug[ru_]] := ru

addDebug[ru:_Rule | _RuleAdd] := {
  DebLine[-2, ##], ru, DebLine[1, ##], DebLine[2, ##]
}&[ debug, DebTag[debtag] @ DebPatt[ru] ]

addDebug[other_] := other


DebLine[__, {}] = {}

DebLine[i_, j_, _] := {} /; BitAnd[i, Abs[j]] === 0

DebLine[i_, _, {var_, tag___}] := DebugLine[i, var[[1]], tag]


DebPatt[ru_[_, expr_Call]] :=
  ru[Cases[expr, _AddrOf, Infinity, 1][[1,1]], expr]

DebPatt[other_] := other


Attributes[DebTag] = {HoldRest}

DebTag[False, ___][_] = {}

DebTag[True, ___][ru_] := {ru}

DebTag[tag_String, ___][ru_] := {ru, tag}

DebTag[tag_][ru_] := DebTag[tag[ru], tag[ru], tag][ru]

DebTag[t_, t_, tag_][ru_] := {ru, tag}

DebTag[tag_, __][ru_] := {ru, tag}


Attributes[Prep] = {Listable}

Prep[DoLoop[expr_, ind__]] := doloop[process[expr], ind]

Prep[i_IndexIf] := MapIf[process, i]

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

Hoist[patt_][expr_, i__] :=
Block[ {veto, abb, got = {}, c},
  veto = Alt@@ Union[
    DoIndex/@ Cases[expr, DoLoop[_, j__] :> j, Infinity],
    Cases[expr, SumOver[j_, ___] :> j, Infinity] ];
  abb = Cases[expr, x:patt /; FreeQ[x, veto], Infinity] //Union;
  abb = (got = Join[got, c = Complement[#2, got]]; #1 -> c)&@@@
    SortBy[hsel[abb]/@ {i}, Length[ #[[2]] ]&];
  Fold[hdo, expr, ReplacePart[abb, {1,2} -> {}]]
]


hsel[abb_][it:{s_, ___}] := it -> Select[abb, !FreeQ[#, s]&]

hsel[abb_][it_] := it -> Select[abb, !FreeQ[#, it]&]


hdo[DoLoop[expr_, j__], i_ -> {}] := DoLoop[expr, i, j]

hdo[expr_, i_ -> {}] := DoLoop[expr, i]

hdo[expr_, i_ -> abb_] :=
  DoLoop[RotateLeft[Flatten[Fold[htmp, expr, abb]], 1], i]


htmp[expr_, abb_] := htmp[expr, abb,
  Position[expr, _ -> abb, Infinity, 1, Heads -> False]]

htmp[expr_, abb_, {}] := {expr /. abb -> #, # -> abb}& @ tmp[abb]

htmp[expr_, abb_, pos_] := {Delete[expr, pos], Extract[expr, pos]}


ToVars[patt_, stub_String] := ToVars[patt, stub &]

ToVars[patt_List, f_] := ToVars[ToAlt[patt], f]

ToVars[patt_, namefun_][expr_] :=
Block[ {var, mkvar, process, doloop = DoLoop, new, varc = 0},
  var[val_] := var[val] =
    (Sow[#1[++varc, #2 -> val]]; #2)&[#, NewSymbol[#, 0]]& @ namefun[val];

  Attributes[mkvar] = {Listable};
  mkvar[(ru:Rule | RuleAdd)[lhs_, rhs:Except[patt]]] :=
    {Last/@ Sort[Flatten[#2]], #1}&@@
      Reap[ru[lhs, rhs /. p:patt :> var[p]]];
  mkvar[other_] := other;

  process = First[{mkvar[Prep[#1]], DownValues[var] = #2}]&[
    #, DownValues[var] ]&;

  Flatten[{process[expr]}]
]


SplitSums[li_List, wrap___] := SplitSums[Plus@@ li, wrap]

SplitSums[x_, wrap_:Identity] := {wrap[x]} /; FreeQ[x, SumOver]

SplitSums[x_, wrap_:Identity] :=
Block[ {term},
  term[_] = 0;
  assign[Expand[x, SumOver]];
  term[_] =.;
  #[[1,1,1]] wrap[Plus@@ Flatten[ #[[2]] ]]&/@ DownValues[term]
]

assign[p_Plus] := assign/@ p

assign[t_Times] := (term[#1] = {term[#1], #2})&@@ cull/@ t

assign[other_] := term[1] = {term[1], other}

cull[o_SumOver] := {o, 1}

cull[other_] := {1, other}


Attributes[ivalid] = {HoldRest}

ivalid[i_, i_] := Sequence[]

ivalid[_, DoDim[s_]] := s


FindIndices[var_ -> _] :=
  Union[Cases[var, s_Symbol :> ivalid[DoDim[s], DoDim[s]]]]

FindIndices[t_Times] := Cases[t, SumOver[i__] :> {i}]

FindIndices[i_IndexIf] := Union@@ FindIndices/@ Flatten[List@@ i]

_FindIndices = {}


ToDoLoops[li:{__}, indices_:FindIndices] :=
Block[ {do, si},
  _do = {};
  Scan[(do[#1] = {do[#1], Tag[##]})&[indices[unpatt[#]], #]&, li];
  si = Flatten[Cases[DownValues[do], _[_[_[{___}]], a_] :> a]];
  DoLoop[ Last/@ #, Sequence@@ #[[1,1]] ]&/@ 
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
  DeclIf -> False,
  RealArgs -> Level[{PaVeIntegral, CutIntegral, CutMasters,
    Bput, Cput, Dput, Eput, Fput, Log, Sqrt}, {-1}],
  Newline -> "\n" }

WriteExpr[_, _[], ___] = {}

WriteExpr[hh_, CodeExpr[vars_, tmpvars_, expr_], opt___Rule] :=
Block[ {horner, fcoll, ffun, type, tmptype, indtype, declif, rargs, $Newline,
var, decl, $DoLoop, $IndexIf},
  {horner, fcoll, ffun, type, tmptype, indtype, declif, rargs, $Newline} =
    ParseOpt[WriteExpr, opt];
  rargs = Alt[rargs];
  horner = If[ horner =!= True, {}, {h_HoldForm :> h,
    p_Plus :> (HornerForm[p] /. HornerForm -> Identity) /; Depth[p] < 6} ];
  fcoll = If[ TrueQ[fcoll], TermCollect, Identity ];
  _var = {};
  $DoLoop = $IndexIf = 0;
  varType[vars, type, expr];
  varType[tmpvars, tmptype /. Type -> type, expr];
  varType[Union[DoIndex/@
    Cases[expr, DoLoop[_, i__] :> i, Infinity]], indtype, expr];
  _var =.;
  decl = newline[VarDecl[ Flatten[#2], #1[[1,1]] ]&@@@ DownValues[var]];
  Flatten[{WriteBlock[hh, IfDecl[declif, decl, expr]]}]
]

WriteExpr[hh_, expr_, opt___Rule] := WriteExpr[hh,
  $WriteExprDebug = PrepareExpr[expr, FilterOpt[PrepareExpr, opt]],
  FilterOpt[WriteExpr, opt]]


varType[_, False, _] = 0

varType[v:{__}, s:_String | {__String}, _] := var[s] = {var[s], v}

varType[v_List, r:_Rule | {__Rule}, _] :=
  MapThread[varSet[], {v, Replace[v, r, {1}]}]

varType[v_List, f_, expr_] := MapThread[varSet[expr], {v, f/@ v}]


_varSet[v_, v_] = 0

_varSet[v_, s:_String | {__String}] := var[s] = {var[s], v}

varSet[expr_][v_, f_] := varSet[][v, f[Cases[expr, _[v, _]]]]


$DebugPre[1, level_:4] := "#if DEBUG >= " <> ToString[level] <> "\n"

$DebugPost[1] = "#endif\n"

_$DebugPre = _$DebugPost = ""

$DebugCmd[1] = "DEB(\"`1``2`\", `3`)\n"

$DebugCmd[-2] = "CHK_PRE(`3`)\n"

$DebugCmd[2] = "CHK_POST(\"`1``2`\", `3`)\n"


Attributes[DebStatement] = {Listable}

DebStatement[var_, cmd_, tag_] := "!" <>
  ToString[StringForm[cmd, tag, ToDef[var /. HelAll -> Identity], 
    ToDef[var]]]


Attributes[WriteBlock] = {Listable}

WriteBlock[_, Hold[expr_]] := (expr; {})

WriteBlock[hh_, s_String] := (WriteString[hh, s]; s)

WriteBlock[hh_, Newline[s___]] := WriteBlock[hh, s <> $Newline]

WriteBlock[_, DebugLine[_]] = {}

WriteBlock[hh_, DebugLine[i_, var_, tag___]] := WriteBlock[hh, Newline[
  $DebugPre[i],
  DebStatement[var, StringJoin[$DebugCmd[i]],
    StringJoin[({ToString[#], ":"}&)/@ {tag}]],
  $DebugPost[i] ]]

WriteBlock[hh_, DoLoop[expr_, ind__]] := WriteBlock[hh, {
  Hold[++$DoLoop],
  Newline[#1],
  expr,
  Hold[--$DoLoop],
  Newline[#2]
}]&@@ DoDecl[ind]

WriteBlock[hh_, i_IndexIf] := WriteBlock[hh, {
  Hold[++$IndexIf],
  Newline[$CodeIf, ToCode[#1], $CodeThen],
  #2, ElseIf[##3],
  Hold[--$IndexIf],
  Newline[$CodeEndIf]
}]&@@ i

ElseIf[] = ElseIf[True, RuleAdd[_, 0]] = {}

ElseIf[True, a_] := {Newline[$CodeElse], a}

ElseIf[cond_, a_, r___] := {
  Newline[$CodeElseIf, ToCode[cond], $CodeThen],
  a, ElseIf[r] }

WriteBlock[hh_, ru_[var_, {sub__, expr_}]] :=
  WriteBlock[hh, {sub, ru[var, expr]}]

WriteBlock[_, RuleAdd[_, 0]] := Sequence[]

WriteBlock[hh_, RuleAdd[var_, expr_]] := (
  Write[hh, $CodeIndent,
    (var -> HoldForm[var + #])& @ ffun[FExpr[expr]], $CodeEoln];
  WriteString[hh, $Newline];
  var -> var + expr
)

WriteBlock[hh_, var_ -> expr_Call] := (
  Write[hh, $CodeCall, ffun[FExpr@@ expr], $CodeEoln];
  WriteString[hh, $Newline];
  var -> expr
)

WriteBlock[hh_, var_ -> expr_] := (
  Write[hh, $CodeIndent, var -> ffun[FExpr[expr]], $CodeEoln];
  WriteString[hh, $Newline];
  var -> expr
)

WriteBlock[hh_, other_] := (
  Write[hh, FExpr[other]];
  other
)


FExpr[expr_] := fcoll[ expr /.
    Complex[a_, b_] :> a + cI b /.
    Dminus4 -> -2/Divergence /.
    E^x_ :> exp[x] /.
    f:rargs[__] :> RealArgs[f] /.
    Den[p_, m_, d___] :> (p - m)^-d /.
    IGram[x_] :> 1/x /.
    horner ] /.
  {h_HoldForm :> h, Times -> OptTimes} /.
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
  Format[p_Integer^r_Rational, #] := (HoldForm[#]^r &)[ N[p] ];
)&/@ {FortranForm, CForm}

Format[RuleAdd[a_, b_], FortranForm] := SequenceForm[a, " = ", a + b];
Format[RuleAdd[a_, b_], CForm] := SequenceForm[a, " += ", b];

Format[s_Symbol^2, CForm] := HoldForm[s s];
Format[s_Symbol^-2, CForm] := 1/HoldForm[s s]

Protect[Rule, Rational, Power]

Format[HoldCode[x_], FortranForm] := HoldForm[x];
Format[HoldCode[x_], CForm] := HoldForm[x]


OptTimes[t__] :=
Block[ {p = Position[N[{t}], _Real, 1, Heads -> False]},
  OptNum[Times@@ Extract[{t}, p], Times@@ Delete[{t}, p]]
]

OptNum[const_, 1] := const

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
format 75;\n\
auto s sM;\n\
auto i iM;\n\n"

$BlockSize = 700

$FileSize = 30 $BlockSize

$MaxFunctionName = 30

$RecursionLimit = 1024

$PaintSE = False

$PutSE = False

EndPackage[]


If[ $NoModelSpecific =!= True,
  Get[ToFileName[$FormCalcSrc, "ModelSpecific.m"]] ]


(* make Needs["FeynArts`"] work after loading FormCalc *)

If[ !NameQ["FeynArts`$FeynArts"],
  Unprotect[$Packages];
  $Packages = DeleteCases[$Packages, "FeynArts`"];
  Protect[$Packages] ]

Null

