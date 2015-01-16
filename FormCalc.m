(*

This is FormCalc, Version 1.5
Copyright by Thomas Hahn 1999
last modified 28 Oct 99 by Thomas Hahn

Release notes:

FormCalc is free software, but is not in the public domain.
Instead it is covered by the GNU library general public license.
In plain English this means:

1. We don't promise that this software works.
   (But if you find any bugs, please let us know!)

2. You can use this software for whatever you want.
   You don't have to pay us.

3. You may not pretend that you wrote this software.
   If you use it in a program, you must acknowledge
   somewhere in your documentation that you've used
   our code.

If you're a lawyer, you can find the legal stuff at
http://www.fsf.org/copyleft/lgpl.html.

The user guide for this program can be found at
http://www-itp.physik.uni-karlsruhe.de/formcalc

If you find any bugs, or want to make suggestions, or
just write fan mail, address it to:
	Thomas Hahn
	Institut fuer Theoretische Physik
	Universitaet Karlsruhe
	e-mail: hahn@particle.uni-karlsruhe.de

To join the FormCalc mailing list, send a mail (any text) to
	hahn-formcalc-subscribe@particle.uni-karlsruhe.de

Have fun!

*)

Print[""];
Print["FormCalc 1.5"];
Print["by Thomas Hahn"];
Print["last revision: 28 Oct 99"];


BeginPackage["FormCalc`"]

(* FeynArts symbols and their FormCalc counterparts *)

FeynAmpList::usage = "FeynAmpList[info][amps] is the head of a list of
FeynAmp objects. info contains additional information handed down by
FeynArts about the process, the momenta etc."

FeynAmp::usage = "FeynAmp[gname, mom, amp] is the FeynArts way of writing
a Feynman amplitude amp with name gname and integration momentum mom."

GraphName::usage = "GraphName[identifiers] gives the name of a Feynman
diagram."

Insertions::usage = "Insertions[lev][ins] gives a list of insertion rules
for the generic amplitude at level lev."

G::usage = "G[sym][cto][fi][kin] is the FeynArts notation of a generic
coupling constant. It is converted to a symbol like GM1 in FormCalc."

Process::usage = "Process -> {proc} contains the process specification
handed down by FeynArts in the information field of a FeynAmpList."

PropagatorDenominator::usage = "PropagatorDenominator[p, m] is the
FeynArts expression for 1/(p^2 - m^2)."

FeynAmpDenominator::usage = "FeynAmpDenominator[prden..] is the head
wrapped around the PropagatorDenominators inside a loop."

DEN::usage = "DEN[p2, m2] stands for 1/(p2 - m2). Note that in contrast
to PropagatorDenominator p2 and m2 are the momentum and mass _squared_."


Spinor::usage = "Spinor[p, m] is the spinor with momentum p and mass m.
Spinor corresponds to the more conventional way of writing spinors by\n
   Spinor[p, m, 1] ** ...  -> \\bar u\n
   Spinor[p, m, -1] ** ... -> \\bar v\n
   ... ** Spinor[p, m, 1]  -> u\n
   ... ** Spinor[p, m, -1] -> v."

LeptonSpinor = Spinor

QuarkSpinor = Spinor

DiracSlash::usage = "DiracSlash[p] represents p_mu gamma^mu."

DiracMatrix::usage = "DiracMatrix[mu] represents the Dirac matrix with
Lorentz index mu."

ga::usage = "ga[li] represents the gamma matrix with Lorentz index li."

ga5::usage = "ga5 represents gamma_5."

ChiralityProjector::usage = "ChiralityProjector[+-1] represents the
chirality projectors omega_{+-} = (1 +- gamma_5)/2."

omp::usage = "omp represents the right handed chirality projector."

omm::usage = "omm represents the left handed chirality projector."

DiracTrace::usage = "DiracTrace represents the trace over Dirac matrices
in a fermion loop."

SpinorChain::usage = "SpinorChain[Spinor[...], ..., Spinor[...]]
represents an open fermion chain."


Index::usage = "Index[t, n] represents an index of type t with number n."

IndexDelta::usage = "IndexDelta[i1, i2] represents the Kronecker delta."

SumOver::usage = "SumOver[i, r, ext] specifies that the amplitude which
it is multiplied with is to be summed in the index i over the range r.
For an index belonging to an external particle there is a third argument,
External."

MetricTensor::usage = "MetricTensor[mu, nu] represents the metric tensor
with Lorentz indices mu and nu."

FourVector::usage = "FourVector[p, mu] represents the four-vector p with
Lorentz index mu."

PolarizationVector::usage = "PolarizationVector[p, mu] represents the
polarization vector belonging to the momentum p with Lorentz index mu."

e::usage = "e[n] is the nth polarization vector."

Momentum::usage = "Momentum[p] is an obsolete form to represent the
four-momentum p."

k::usage = "k[n] is the nth momentum."

q1::usage = "q1 is the integration momentum in a one-loop amplitude."

ScalarProduct::usage = "ScalarProduct[p1, p2] is the scalar product of two
four-vectors p1 and p2. It is converted to Pair[p1, p2] in FormCalc."

GaugeXi::usage = "GaugeXi[v] is the way FeynArts denotes the gauge
parameter for the gauge boson v. It is converted to xi[v] in FormCalc."

xi::usage = "xi[v] is the gauge parameter of gauge boson v."


SUNT::usage = "SUNT[g, c1, c2] are the generators of SU(N)."

SUNF::usage = "SUNF[g1, g2, g3] are the structure constants of SU(N)."

SUNFSum::usage = "SUNFSum[g1, g2, g3, g4] is short for
Sum[SUNF[g1, g2, i] SUNF[i, g3, g4], i]."

SUND::usage = "SUND[g1, g2, g3] are the totally symmetric d-tensors
of SU(N)."


S::usage = "S is the Mandelstam variable s. If p1 and p2 denote the
incoming momenta, S = (p1 + p2)^2."

T::usage = "T is the Mandelstam variable t. If p1 denotes the first
incoming and k1 the first outgoing momentum, T = (p1 - k1)^2."

U::usage = "U is the Mandelstam variable u. If p1 denotes the first
incoming and k2 the second outgoing momentum, U = (p1 - k2)^2."

Sf::usage = "Sf is an extended Mandelstam variable for a 2 -> 3 reaction.
If k1 and k2 denote the first two outgoing momenta, Sf = (k1 + k2)^2."

Tf::usage = "Tf is an extended Mandelstam variable for a 2 -> 3 reaction.
If p2 denotes the second incoming and k2 the second outgoing momentum,
Tf = (p2 - k2)^2."

Uf::usage = "Uf is an extended Mandelstam variable for a 2 -> 3 reaction.
If p2 denotes the second incoming and k1 the first outgoing momentum,
Uf = (p2 - k1)^2."



(* the main functions: OneLoop, ProcessFile, and their options *)

OneLoop::usage = "OneLoop[amps] is the basic function for evaluating
Feynman diagrams. It takes as an argument either a single diagram (with
head FeynAmp) or a list of diagrams and returns the results in a collected
and abbreviated form."

AmplitudeLevel::usage = "AmplitudeLevel is an option of OneLoop. It is
used in amplitudes generated with FeynArts 2 only, and specifies the level
(Classes or Particles) at which to calculate the amplitudes. The default
setting Automatic takes the deepest level available."

ChiralME::usage = "ChiralME is an option of OneLoop. It specifies whether
spinor chains are returned in terms of the chirality projectors omega_+
and omega_-, or in the vector/axial vector decomposition 1 and gamma_5.
The former is the default since it yields much more compact helicity
matrix elements."

DiracSimplify::usage = "DiracSimplify is an option of OneLoop. With
DiracSimplify -> True (the default), OneLoop simplifies fermionic matrix
elements by replacing in turn all momenta via momentum conservation and
then using the Dirac equation. However, this procedure can be rather slow,
so DiracSimplify -> False can be used to turn it off."

DotSimplify::usage = "DotSimplify is an option of OneLoop. When set to
True (the default), FormCalc tries to reduce the number of independent
matrix elements by eliminating one momentum via momentum conservation in
dot products and epsilon tensors."

NoExpand::usage = "NoExpand is an option of OneLoop. NoExpand -> {sym1,
sym2, ...} specifies that sums containing any of sym1, sym2, ... are
expanded in FORM."

EditCode::usage = "EditCode is an option of OneLoop and HelicityME. It
edits the temporary file passed to FORM using $Editor and is of course
used only for debugging."

RetainFile::usage = "RetainFile is an option of OneLoop and HelicityME.
When set to True, it prevents removing the temporary file which contains
the FORM input after running FORM."

ProcessFile::usage = "ProcessFile[amps, outtag] evaluates amps using
OneLoop. The results are split into blocks according to the Classification
option, but at least such that each block contains only diagrams with the
same index summations. These blocks are then written to files whose
names are outtag <> an identifier for the block. If a file name is given
for amps, the corresponding file is loaded.\n
ProcessFile[amps] does the same as ProcessFile[amps, outtag] but returns
the blocks in a list instead of writing them to files."

Classification::usage = "Classification is an option of ProcessFile.
It can take the values IndexSumsOnly, Standard (default), and Tough which
control how ProcessFile divides the diagrams into different classes.
In all cases, diagrams with different index summations are separated,
i.e. each class contains only diagrams which have the same index
summations. For Standard, the amplitude is split into bosonic and
fermionic parts, too, and for Tough also into different topologies and (if
applicable) different fermion families. The latter is mainly used for
large amplitudes."

IndexSumsOnly::usage = "IndexSumsOnly is a value the Classification
option of ProcessFile can take. The criterium for two diagrams to be put
into different classes is in this case that they have different index
summations."

Standard::usage = "Standard is a value the Classification option of
ProcessFile can take. The criterium for two diagrams to be put into
different classes is in this case that they either have different index
summations or are bosonic and fermionic."

Tough::usage = "Tough is a value the Classification option of ProcessFile
can take. The criterium for two diagrams to be put into different classes
is in this case that they either have different index summations, or are
bosonic and fermionic, or belong to different topologies or (if they are
fermionic) fermion families. This is not so much a meaningful as a
practical criterium to break down large amplitudes."

FermionFamily::usage = "FermionFamily[patt] = \"string\" returns the
identifier appended to the filename by ProcessFile if Classification ->
Tough is set, the amplitude is fermionic, and it contains patt. Typically
used like FermionFamily[ME | MM | ML] = \"e\" which gives files containing
leptonic amplitudes an extra e in the name."


(* internal stuff (some of it needs to be visible because of ReadForm) *)

ReadForm::usage = "ReadForm[\"file\"] reads file which it assumes to
contain FORM output into Mathematica.\n
ReadForm[\"!cmd\"] executes cmd and pipes its output back into
Mathematica."

PowerCountingFor::usage = "PowerCountingFor[momlist] tells ReadForm to
multiply each abbreviation in the FORM result by $Scale^n, where n is the
number of momlist members appearing in the term. Usually, momlist is the
list of momenta, e.g. {p1, p2, k1, k2}. This function is intended for
internal use only."

FromForm::usage = "FromForm is used internally to reverse substitutions
that were necessary for FORM."

a0::usage = "a0[m] is used internally to represent a pre-form of the
one-point function."

b0::usage = "b0[p, m1, m2] is used internally to represent a pre-form of
the two-point function."

c0::usage = "c0[p1, p2, m1, m2, m3] is used internally to represent a
pre-form of the three-point function."

d0::usage = "d0[p1, p2, p3, m1, m2, m3, m4] is used internally to
represent a pre-form of the four-point function."

e0::usage = "e0[p1, p2, p3, p4, m1, m2, m3, m4, m5] is used internally to
represent a pre-form of the five-point function."

(* more internal symbols: *)

{ abb, abbsum, mat, fme, sun, r2, amp3, Unhide, d$, e$, gi$, ncm,
  pave4, pave3, pave2, pave1, B0m, B1m, B00m, B11m,
  IND1, IND2, IND3, IND4, IND1c, IND2c, IND3c, IND4c }


(* abbreviationing-related functions *)

Abbreviations::usage = "Abbreviations[] returns a list of all
abbreviations introduced so far. It is typically used at the end of a
calculation to save the abbreviations with a command like
Abbreviations[] >> abbrfile."

UseAbbreviations::usage = "UseAbbreviations[abbr] registers the
abbreviations abbr with FormCalc. abbr may be either the abbreviations
themselves (a list of rules) or a filename containing them. This is useful
when continuing a former FormCalc session. UseAbbreviations requires that
the value of $Scale that was used in the former session is also set in
the current session."

OptimizeAbbreviations::usage = "OptimizeAbbreviations[abbr] optimizes the
set of abbreviations returned by Abbreviations[] by eliminating common
subexpressions."

Scale::usage = "Scale is a scale introduced via $Scale to scale matrix
elements, e.g. to make them dimensionless. The variable Scale itself
appears in the Abbreviations[]."

Pair::usage = "Pair[a, b] represents the Minkovskian scalar product of
the four-vectors a and b."

Eps::usage = "Eps[a, b, c, d] represents -I times the antisymmetric
Levi-Civita tensor: Eps[a, b, c, d] = -I a[mu] b[nu] c[rho] d[sigma]
epsilon[mu, nu, rho, sigma], where a, b, c, d are four-vectors and the
sign convention is epsilon[0, 1, 2, 3] = +1."

o1::usage = "o1 is a head wrapped around the prefactor of an abbreviation
in an amplitude. Its default value is o1 = Identity (doing nothing). To
get the shortest possible amplitude, set o1 = Simplify or similar."

o2::usage = "o2 is a head wrapped around a linear combination of
abbreviations. Its default value is o2 = Identity (doing nothing). To get
the shortest possible amplitude, set o2 = Simplify or similar."


(* miscellaneous functions *)

DiagramType::usage = "DiagramType[diag] returns 2 for a self-energy
diagram, 1 for a vertex diagram, and 0 for a box. In other words, it
returns the number of denominators not containing the integration
momentum."

FermionicQ::usage = "FermionicQ[diag] gives True for a diagram containing
fermions and False otherwise."

Small::usage = "Small[sym] = 0 tells FormCalc to put sym = 0 except when
it appears in negative powers or in loop integrals."

PickLevel::usage = "PickLevel[lev][amps] selects amplitude level lev from
amps."

Pick::usage = "Pick[amp, graphs] picks diagrams from amp where graphs is
of the form {5} which picks diagram # 5, {{5, 9}} which picks diagrams # 5
through 9, or combinations of these. For example, Pick[amp, {6, {10, 12},
27}] picks diagrams # 6, 10, 11, 12, and 27 from amp."


(* FeynCalc compatibility functions *)

FeynCalcGet::usage = "FeynCalcGet[mask] reads files produced with
FeynCalc. mask is taken as input to the Mathematica function FileNames, so
it might be FeynCalcGet[\"file.m\"] or FeynCalcGet[\"*.m\"] or
FeynCalcGet[\"*.m\", \"~/feyncalcfiles\"]."

FeynCalcPut::usage = "FeynCalcPut[expr, file] writes expr to file in
FeynCalc format."

C0::usage = "C0[p1, p2, p1p2, m1, m2, m3] is the scalar three-point
Passarino-Veltman function as used by FeynCalc. It is converted to C0i in
FormCalc."

D0::usage = "D0[p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the scalar
four-point Passarino-Veltman function as used by FeynCalc. It is converted
to D0i in FormCalc."

PaVe::usage = "PaVe[ind, {pi}, {mi}] is the generalized Passarino-Veltman
function as used by FeynCalc. It is converted to C0i or D0i in FormCalc."


(* one-loop integrals *)

A0::usage = "A0[m] is the one-point scalar Passarino-Veltman function
where m is the mass squared."

B0::usage = "B0[p, m1, m2] is the scalar two-point Passarino-Veltman
function where p is the external momentum squared and m1 and m2 are the
masses squared."

B1::usage = "B1[p, m1, m2] is the tensor two-point Passarino-Veltman
function B_1 where p is the external momentum squared and m1 and m2 are
the masses squared."

B00::usage = "B00[p, m1, m2] is the tensor two-point Passarino-Veltman
function B_00 where p is the external momentum squared and m1 and m2 are
the masses squared."

B11::usage = "B11[p, m1, m2] is the tensor two-point Passarino-Veltman
function B_11 where p is the external momentum squared and m1 and m2 are
the masses squared."

DB0::usage = "DB0[p, m1, m2] is the derivative of B0[p, m1, m2] with
respect to p."

DB1::usage = "DB1[p, m1, m2] is the derivative of B1[p, m1, m2] with
respect to p."

DB00::usage = "DB00[p, m1, m2] is the derivative of B00[p, m1, m2] with
respect to p."

DB11::usage = "DB11[p, m1, m2] is the derivative of B11[p, m1, m2] with
respect to p."

C0i::usage = "C0i[id, p1, p2, p1p2, m1, m2, m3] is the generic three-point
Passarino-Veltman function which includes both scalar and tensor integrals
specified by id. For example, C0i[cc0, ...] is the scalar function C0,
C0i[cc112, ...] the tensor function C_112 etc. Call the external momenta
k1...k3, then the arguments are given as p1 = k1^2, p2 = k2^2,
p1p2 = (k1 + k2)^2, and m1...m3 are the masses squared."

D0i::usage = "D0i[id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the
generic Passarino-Veltman four-point function which includes both scalar
and tensor integrals specified by id. For example, D0i[dd0, ...] is the
scalar function D0, D0i[dd1233, ...] the tensor function D_1233 etc.
Call the external momenta k1...k4, then the arguments are given as
p1 = k1^2, p2 = k2^2, p3 = k3^2, p4 = k4^2, p1p2 = (k1 + k2)^2,
p2p3 = (k2 + k3)^2, and m1...m4 are the masses squared."


(* helicity matrix elements *)

HelicityME::usage = "HelicityME[plain, conj] calculates the helicity
matrix elements for all combinations of spinor chains that appear in the
expression (plain conj^*). Terms of this kind arise in the calculation of
the squared matrix element, where typically plain is the one-loop result
and conj the Born expression. The arguments don't even have to be
amplitudes since they are only used to determine which spinor chains to
select from the abbreviations. The symbol All can be used to select all
spinor chains currently defined in the abbreviations."

All::usage = "All is a possible input value for HelicityME, indicating
that all spinor chains currently defined in the abbreviations should be
used instead of just those appearing in a particular expression."

AbbreviationsToUse::usage = "AbbreviationsToUse is an option of
HelicityME. It specifies which abbreviations are used to calculate the
helicity matrix elements."

Hel::usage = "Hel[i] denotes the helicity of the ith external particle.
It can take the values +-1."

s::usage = "s[n] is the nth helicity reference vector."

Mat::usage = "Mat[Fi SUNi] denotes a matrix element in an amplitude. This
is an auxillary symbol and it used to aid FormCalc when calculating the
squared matrix elements with HelicityME.\n
Mat[Fi, Fj] denotes a squared matrix element composed from the simple
matrix elements Fi and Fj."

SquaredME::usage = "SquaredME[plain, conj] returns the matrix element
plain Conjugate[conj]. This performs a nontrivial task only for fermionic
amplitudes: the product of two fermionic amplitudes\n
    M1 = a1 F1 + a2 F2 + ... and\n
    M2 = b1 F1 + b2 F2 + ... is returned as\n
    M1 M2^* = a1 b1^* Mat[F1, F1] + a2 b1^* Mat[F2, F1] + ...\n
The special case of plain === conj can be written as SquaredME[plain]
which is of course equivalent to SquaredME[plain, plain]."

RealQ::usage = "RealQ[sym] is True if sym represents a real quantity
which means in particular that Conjugate[sym] = sym."


(* system variables *)

$Editor::usage = "$Editor specifies the editor used in debugging FORM
code."

$TempFile::usage = "$TempFile is the name for the temporary file used to
write out FORM code in OneLoop. Unless an error occurs in FORM or
RetainFile -> True is used, it is deleted after usage."

$OnShell::usage = "$OnShell specifies whether the external particles may
be assumed on shell. Caution: This variable must be set *before* loading
amplitudes."

$FormCalcDir::usage = "$FormCalcDir is the directory in which FormCalc
lives."

$Platform::usage = "$Platform is a string that identifies the platform
FormCalc is running on. It is the value of the environment variable
HOSTTYPE in the Bourne shell (sh) and is used to distinguish the ReadForm
binaries of different platforms."

$FormCmd::usage = "$FormCmd gives the name of the actual FORM executable.
It may contain a path."

$Transversality::usage = "$Transversality specifies whether FormCalc may
assume transversality for polarization vectors, i.e. e[i].k[i] = 0."

$Scale::usage = "The abbreviations introduced by FormCalc are scaled
with $Scale for every momentum they include, e.g. e1.k2 is abbreviated
by Pair1/$Scale, e1.k2 e2.k3 by Pair1 Pair2/$Scale^2, etc. The idea is
to make the abbreviations dimensionless by cleverly choosing a scale
factor, e.g. $Scale = Sqrt[S] for a 2 -> 2 process."

$Dimension::usage = "$Dimension specifies the dimension FormCalc works in.
It can be D for dimensional regularization and 4 for dimensional reduction
and constrained differential renormalization."


Begin["`Private`"]

$FormCalcDir =
  If[ FileType[$Input] === File, $Input,
	(* if FormCalc was loaded from a directory in $Path: *)
    Block[ {full},
      Scan[
        If[ FileType[full = # <> "/" <> $Input] === File, Return[full] ]&,
        $Path ] ]
  ]

Block[ {pos = StringPosition[$FormCalcDir, "/"]},
  If[Length[pos] === 0, $FormCalcDir = "",
    $FormCalcDir =
      SetDirectory[StringTake[ $FormCalcDir, pos[[-1, -1]] ]] <> "/";
    ResetDirectory[] ]
]


$Platform = Environment["HOSTTYPE"]

If[ Head[$Platform] =!= String, $Platform = "" ]


If[ $LinkSupported,
  Install[$FormCalcDir <> "ReadForm_" <> $Platform],
(* else *)
  Print["WARNING: Your Mathematica kernel does not support MathLink."];
  Print["This means that the main functions OneLoop, ProcessFile, and"];
  Print["HelicityME will not work."];
  ReadForm[_] = $Failed ]


Off[General::spell1, General::spell, Unset::norep]


(* generic functions *)

ToSymbol[x__] := ToExpression[ StringJoin[ToString/@ Flatten[{x}]] ]


Pick[amp_, graphs_] :=
  amp[[ Flatten[ If[Head[#] === List, Range@@ #, #]&/@ graphs /.
    a_Integer (Repeated | RepeatedNull)[b_] :>
      Range@@ Sort[Floor[{b, a}]] ]//Union ]]


Attributes[DiagramType] = {Listable}

DiagramType[a_FeynAmp] := Exponent[a[[3]] /. DEN[__] -> DEN, DEN]


FermionicQ[a_] :=
  !FreeQ[a, Spinor | ga | ga5 | omp | omm | ChiralityProjector |
    DiracMatrix | DiracSlash | DiracTrace | SpinorChain]


TakeGraph[gr_ -> _] = gr

TakeGraph[gr_] = gr

PickLevel[lev_][ FeynAmp[n__, a_, coup_ -> ins_] ] :=
Block[ {sel},
  sel = TakeGraph/@
    Cases[{ins}, Insertions[lev][rulz__] :> rulz, Infinity];
  If[ Length[coup] === 0, FeynAmp[n, Length[sel] a],
    FeynAmp[n, a, coup -> Insertions[lev]@@ sel] ]
]

PickLevel[_][a_FeynAmp] = a

PickLevel[lev_][a_] := PickLevel[lev]/@ a


PropagatorDenominator[p_, m_] :=
  PropagatorDenominator[-p, m] /; !FreeQ[p, -q1]

PropagatorDenominator[0, m_] = -1/m^2

FeynAmpDenominator[ a___, p:PropagatorDenominator[q1, _], b___ ] :=
  den[p, b, a]

den[ a___, PropagatorDenominator[q_, m1_], b___,
  PropagatorDenominator[q_, m2_], c___ ] :=
  1/Factor[m1^2 - m2^2] (den[a, PropagatorDenominator[q, m1], b, c] -
    den[a, PropagatorDenominator[q, m2], b, c])

den[ PropagatorDenominator[q1, m1_] ] :=
  I Pi^2 a0[m1^2]

den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_] ] :=
  I Pi^2 b0[p2 - p1, m1^2, m2^2]

den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_] ] :=
  I Pi^2 c0[p2 - p1, p3 - p1, m1^2, m2^2, m3^2]

den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_],
     PropagatorDenominator[p4_, m4_] ] :=
  I Pi^2 d0[p2 - p1, p3 - p1, p4 - p1, m1^2, m2^2, m3^2, m4^2]

den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_],
     PropagatorDenominator[p4_, m4_],
     PropagatorDenominator[p5_, m5_] ] :=
  I Pi^2 e0[p2 - p1, p3 - p1, p4 - p1, p5 - p1,
            m1^2, m2^2, m3^2, m4^2, m5^2]


Attributes[IndexDelta] = {Orderless}

IndexDelta[n_, n_] = 1

IndexDelta[_Integer, n_Integer] = 0


(* Fermion stuff: Dirac algebra etc *)

(* DiracSlash[p_] := -DiracSlash[-p] /; !FreeQ[p, -q1] *)


Spinor[(s_Integer:1) p_, m_] := Spinor[p, m, s]

Spinor[p_, m_, s_] := Spinor[p, Small[m], s] /; Head[Small[m]] =!= Small

SpinorType[-1] = "v";
SpinorType[1] = "u"


Unprotect[NonCommutativeMultiply]

Format[ Spinor[p1_, m1_, s1_] ** c___ ** Spinor[p2_, m2_, s2_] ] :=
  SequenceForm[ "(",
    ColumnForm[{"_", SpinorType[s1][p1, m1]}, Left, Above] **
      c ** SpinorType[s2][p2, m2], ")" ]

___ ** 0 ** ___ = 0

Protect[NonCommutativeMultiply]


Unprotect[Dot]

Format[ sp_Spinor . c__ ] := sp ** c

Protect[Dot]

ga[0] = 0


(* Boson stuff: 4-vectors etc *)

GaugeXi = xi

ScalarProduct = Pair


FourVector[p_Plus, li_] := FourVector[#, li]&/@ p

FourVector[n_?NumberQ k_, li_] := n FourVector[k, li]

FourVector[k_ x_xi, li_] := x FourVector[k, li]

FourVector[k_, li_] := k[li]


Attributes[MetricTensor] = {Orderless}

MetricTensor/: MetricTensor[li__]^2 :=
  MetricTensor[li] ** MetricTensor[li]


Momentum[p_] = p


Conjugate[PolarizationVector] ^= PolarizationVector

PolarizationVector[n_?NumberQ p_, li_] := n PolarizationVector[p, li]

PolarizationVector[p_, li_] :=
  PolarizationVector[p, li] = ToSymbol["e", p][li]



(* Reading and analysing FeynArts files *)

FormMandelstam = FormDotSimplify = FormDiracSimplify = ""

FromFormRules = FormVectors = FormIndices = IndexRanges = {}


SetPair22[p1_, p2_, sign_, kin_] :=
  p1/: Pair[p1, p2] = sign kin - sign MomSquare[p1] - sign MomSquare[p2]

SetPair12[mom_, i_] :=
Block[ {lhs = Pair@@ Delete[mom, i], rhs},
  rhs = Plus@@
    Array[ If[# == i, 1, -1]/2 MomSquare[ mom[[#]] ] &, Length[mom]];
  If[MatchQ[lhs, -_], lhs = -lhs; rhs = -rhs];
  TagSet@@ {lhs[[1]], lhs, rhs}
]


id[0, _] = {}

id[lhs_, rhs_] :=
  { "id ", ToString[lhs//InputForm], " = ",
           ToString[rhs//InputForm], ";\n" }


FeynAmpList::noFAfile = "Hold it! This is no FeynArts amplitude."

FeynAmpList[h__][a___] :=
Block[ {amps, proc, mom, ps, eps, sf, tf, uf, mandel, dotprods = {},
dot, iddot, Index, PropagatorDenominator, MetricTensor = "d_"},
  proc = Process /. {h};
  If[ Head[proc] =!= Rule,
    Message[FeynAmpList::noFAfile];
    Return[$Failed] ];

  mom = #[[2]]&/@ Join[ proc[[1]], -proc[[2]] ];
  proc = Select[Join@@ proc, Head[ #[[2]] ] === Symbol &];

  Scan[Clear, ps = #[[2]]&/@ proc];
  Clear[U];

  eps = ToSymbol["e", #]&/@ ps;
  FormVectors = Flatten[{q1, eps, ps}];
  PowerCountingFor[ps];
  FromFormRules = Join[
    MapIndexed[#1 -> k@@ #2 &, ps], MapIndexed[#1 -> e@@ #2 &, eps] ];
  Block[ {Set},
    (FromForm[x_] := Block[#, ToExpression[x]])& @
      Append[Apply[Set, FromFormRules, 1], Dot = PairAbbr]
  ];

  Index[t_, n_] := Index[t, n] =
    ToSymbol[StringTake[ToString[t], 3], n];

  amps = {a} /. Conjugate[ep:(Alternatives@@ eps)[__]] -> ep;

  FormIndices = Cases[DownValues[Index], _[_, s_Symbol] -> s];
  IndexRanges =
  Block[ {dim},
    amps /. SumOver[x1_, x2_, ___] :> (dim[x1] = x1 -> x1 == x2);
    dim[x_] = x -> x == 0;
    dim/@ FormIndices ];

	(* this is p^2 = m^2: *)
  If[ $OnShell,
    mandel = Plus@@ Apply[UpSet@@ {MomSquare[#2], #3^2} &, proc, 1] ];

	(* introduce kinematical invariants: *)
  Switch[ Length[ps],
    2,
      amps = amps /. Rule@@ Reverse[ps],
    3,
      SetPair12[mom, 1];
      SetPair12[mom, 2];
      SetPair12[mom, 3],
    4 | 5,
      {sf, tf, uf} =
        If[ Length[ps] === 4,
          If[ $OnShell,
            U/: S + T + U = mandel /.
              s_Symbol :> Small[s] /. Small -> Identity ];
          {S, T, U},
        (* else *)
          {Sf, Tf, Uf} ];
      SetPair22[ ps[[1]], ps[[2]], 1/2, S ];
      SetPair22[ ps[[3]], ps[[1]], -1/2, T ];
      SetPair22[ ps[[4]], ps[[1]], -1/2, U ];
      SetPair22[ ps[[3]], ps[[4]], 1/2, sf ];
      SetPair22[ ps[[4]], ps[[2]], -1/2, tf ];
      SetPair22[ ps[[3]], ps[[2]], -1/2, uf ]
  ];

  PropagatorDenominator[p_, m_] :=
    PropagatorDenominator[p, m] = DEN[MomSquare[p], m^2];

  FormDiracSimplify = FormDotSimplify = "";

  mandel = Cases[ UpValues[#] /. Pair -> Dot,
    (_[p_] :> sp_) :> id[p, sp] ]&/@ ps;

  If[ $Transversality,
    dotprods = MapThread[(dot[#1][#2] = 0; id[#1 . #2, 0])&, {eps, ps}] ];

  If[ Length[ps] > 2,
    mom = Array[{-Plus@@ Delete[mom, #], mom[[#]]}&, Length[mom]];
    mom = Join[mom, -mom];
    proc = Cases[mom, {pb_, pa_Symbol} -> {pa, pb}];
    sf = Alternatives@@ ps;

    If[ Length[ps] > 4,
	(* in this case elim last mom completely *)
      AppendTo[dotprods, {mandel[[-1]], id@@ proc[[-1]]}];
      mandel = Drop[mandel, -1];
	(* use System` so that x is not declared as FORM symbol later: *)
      Apply[ (Pair[#1, System`x_] ^= Pair[#2, System`x])&, proc[[-1]] ],
    (* else *)
      If[ $Transversality,
        dot[x1_][x2_] := x1 . x2;
        iddot[pp_] := id@@ (pp /. p:sf :> dot[#][p])&/@ eps,
      (* else *)
        iddot[_] = {}
      ];
      FormDotSimplify = Apply[
        { iddot[{##}],
          id["e_"["J1?,J2?,J3?", #1], #2 /. p:sf :> "e_"["J1,J2,J3", p]],
          ".sort\n" }&,
        proc, 1 ]//StringJoin;
    ];

    FormDiracSimplify = Apply[
      { id[ga["C?", #1], #2 /. p:sf :> ga["C", p]],
        "#call DiracEquation{}\n" }&,
      proc, 1 ]//StringJoin;

	(* set the momentum-conservation rules: *)
    Apply[
      If[(sf = Select[#1, AtomQ, 1]) =!= 0, TagSet@@ {sf, ##}]&,
      mom, 1 ];
  ];
  FormMandelstam = dotprods <> mandel;

  amps /. (1/xi[f_])^r_Rational :> xi[f]^-r
]


GList = Array[ToSymbol["Global`GM", #]&, 50]

DeleteIns[ins_List] := Delete[ins, p]

DeleteIns[ins_] := DeleteIns/@ ins

FeynAmp[g_GraphName, ___, amp_, coup_ -> ins_] :=
Block[ {t, u, p, p2, r, c = 0, coupvar = Take[GList, Length[coup]]},
  t = Transpose[{ins} /. Rule | _Insertions -> Sequence];
  u = Union/@ t;
  p = Position[u, {_}, {1}];
  Apply[(coupvar[[#]] = u[[#, 1]])&, p, 1];
  t = ReplacePart[t, 0, p];
  While[ Length[t] > 1,
    u = t[[1]];
    t = Rest[t];
    ++c;
    If[ u =!= 0,
      If[Length[r = Position[t, -u, {1}]] =!= 0,
        t = ReplacePart[t, 0, r];
        coupvar = ReplacePart[coupvar, -coupvar[[c]], r += c];
        p = Join[p, r] ];
      If[Length[r = Position[t, u, {1}]] =!= 0,
        t = ReplacePart[t, 0, r];
        coupvar = ReplacePart[coupvar, coupvar[[c]], r += c];
        p = Join[p, r] ]
    ]
  ];
  FeynAmp[ g[[-1]], g[[2]],
    amp /. Thread[coup -> coupvar],
    DeleteIns[coupvar -> ins] ]
]

FeynAmp[g_GraphName, ___, amp_] :=
  FeynAmp[ g[[-1]], g[[2]], amp /.
    G[_][cto_][fi__][k_] :>
      (G[cto, fi, kinid[k]] /. List -> Sequence) ]

kinid[k_] := kinid[k] =
Block[ {DiracMatrix = "ga", DiracSlash = "ga", ChiralityProjector,
MetricTensor = "g", FourVector = List, ScalarProduct = List},
  ChiralityProjector[1] = "omp";
  ChiralityProjector[-1] = "omm";
  ToExpression[
    StringJoin[ToString/@
      DeleteCases[Flatten[k //. a_[b__] :> {a, b} /; a =!= List],
        sym_Symbol /; Context[sym] === "System`"]] ]
]


(* FORM interface and main function *)

FormWrite[ FeynAmp[gname_, _, amp_, rulz___] ] :=
Block[ {Pair = Dot, na = gname, fline = 0, ftrace = {}, diags, noins,
DiracTrace = FormDiracTrace, SpinorChain = FormSpinorChain},
  WriteString["stdout", na, " "];
  If[ MemberQ[amplist, na], na = Unique[ToString[na] <> "n"] ];
  AppendTo[amplist, na];
  noins = Length[{rulz}] === 0;
  WriteString[hh, "g ", If[ noins, na, na@@ rulz[[1]] ], " =\n  "];
  Write[hh, amp /. d_Dot :> FormSpinorChain@@ d /; FermionicQ[d] /.
    NonCommutativeMultiply[a_] -> a, ";"];
  AppendTo[alltrace, ftrace];
  If[ noins, WriteString[hh, ".store\n\n"],
    WriteString[hh, ".sort\n\n"];
    diags = Plus@@ Apply[na, Map[Replace[#, ins]&, rulz[[2]], {2}], 1];
    na = ToSymbol["Ins", na];
    inssum += na;
    AppendTo[names, na -> diags] ];
]


chp[1] = "g6_";
chp[-1] = "g7_"

FormDiracTrace[expr__] :=
Block[ {DiracMatrix, Dot = NonCommutativeMultiply},
  AppendTo[ftrace, ++fline];
  DiracMatrix[5] = "g5_"[fline];
  DiracMatrix[li_] = "g_"[fline, li];
  NonCommutativeMultiply[expr] /.
    f1_. DiracSlash[k_] + m_. :>
      f1 Map[If[Head[#] === Symbol, "g_"[fline, #], #]&, k, {-1}] +
        m "gi_"[fline] /.
    f1_. ChiralityProjector[1] + f1_. ChiralityProjector[-1] + rr_. :>
      (f1 + rr) "gi_"[fline] /.
    f1_. ChiralityProjector[1] + f2_. ChiralityProjector[-1] + rr_. :>
      f1/2 "g6_"[fline] + f2/2 "g7_"[fline] + rr "gi_"[fline] /.
    f1_. ChiralityProjector[pm_] + m_. :>
      f1/2 chp[pm][fline] + m "gi_"[fline]
]

FormSpinorChain[expr__] :=
Block[ {DiracMatrix, DiracSlash, ChiralityProjector,
Dot = NonCommutativeMultiply},
  extferm = "1";
  ++fline;
  DiracMatrix[5] = omp[fline]/2 - omm[fline]/2;
  DiracMatrix[li_] = ga[fline, li];
  DiracSlash[k_] :=
    Map[If[Head[#] === Symbol, ga[fline, #], #]&, k, {-1}];
  ChiralityProjector[1] = omp[fline];
  ChiralityProjector[-1] = omm[fline];
  NonCommutativeMultiply[expr]
]


FormDecl[type_, _[f_, v___]] :=
Block[ {l, ll, s},
  ll = StringLength[s = type <> ToString[f]];
  WriteString[hh,
    Fold[
      ( l = StringLength[s = ToString[#2]] + 2;
        {#1, If[(ll += l) > 70, ll = l; ",\n  ", ", "], s} )&,
      s, {v} ] <> ";\n"]
]


FormPattern[ _[_[ sym_Symbol ]], val_ ] :=
  (vars = {vars, sym}; id[sym, val])

FormPattern[ _[_[ (h_Symbol)[args__] ]], val_ ] :=
Block[ {patt, lhs, c1 = 0, c2 = 0},
  lhs = h[args] /. Pattern -> Arg2 /.
    _Blank :> "Patt" <> ToString[++c1] <> "?" /.
    _BlankSequence | _BlankNullSequence :>
      StringJoin[Table["?", {++c2}]];
  vars = {vars, ToSymbol["Patt", #]&/@ Range[c1], Cases[lhs, _Symbol]};
  AppendTo[func, h];
  id[lhs, val]
]

Arg2[_, x_] = x


DeclareVars[expr_, lorentzind_] :=
Block[ {theexpr, vars, indx, func = {}},
  theexpr = {expr, Last/@ Flatten[UpValues/@ FormVectors]};
  vars = Cases[theexpr, x_Symbol /; Context[x] =!= "System`",
    Infinity, Heads -> True];
  smalls = StringJoin[ Apply[FormPattern, DownValues[Small], 1] ];
  vars = Complement[Flatten[{vars, Pi}], FormVectors];
  func = Union[func,
    Select[Flatten[{vars, Conjugate, Sqrt}], !FreeQ[theexpr, #[___]]&] ];
  indx = Intersection[vars,
    Join[FormIndices, ToExpression/@ Names[lorentzind]]];
  vars = Complement[vars, func, indx];
  FormDecl["i ", indx /. IndexRanges];
  FormDecl["v ", FormVectors];
  FormDecl["s ", vars];
  FormDecl["cf ", DeleteCases[func,
    FeynAmp | Spinor | DiracMatrix | DiracSlash |
    ChiralityProjector | DiracTrace | tr]];
]


DoTrivialSums[ FeynAmp[n_, top_, a_, coup_ -> ins_] ] :=
Block[ {amp},
  amp = DoTrivialSumsAmp[a, ins];
  FeynAmp[n, top, amp, coup -> DoTrivialSumsIns/@ ins]
] /; !FreeQ[{a, ins}, SumOver[_, _]]

DoTrivialSums[ FeynAmp[n_, top_, a_] ] :=
  FeynAmp[n, top, DoTrivialSumsAmp[a]] /; !FreeQ[a, SumOver[_, _]]

DoTrivialSums[ a_ ] = a

DoTrivialSumsIns[{ins___, relcf_}] :=
  {ins, relcf /. SumOver[i_, v_] :> v /; FreeQ[{amp, ins}, i]}

DoTrivialSumsIns[a_] = a

DoTrivialSumsAmp[a_, ins___] :=
  Fold[
    If[FreeQ[ {#1, ins}, #2[[1]] ], #1 #2[[2]], #1 #2]&,
    DeleteCases[a, SumOver[_, _]],
    Cases[a, SumOver[_, _]] ]


TakeIns[FeynAmp[__, _ -> ins_]] := List@@ ins

TakeIns[_] = {}


kinobjs = Spinor | ChiralityProjector | DiracMatrix | DiracSlash |
  DiracTrace | SpinorChain | q1 | Global`p1 | Global`p2 | Global`p3 |
  Global`k1 | Global`k2 | Global`k3 | Global`k4 | _String

IList = Array[ToSymbol["Global`Ins", #]&, 500]

SUNobjs = SUNT | SUNF | SUNFSum | SUND


Options[OneLoop] = {
  AmplitudeLevel -> Automatic,
  DiracSimplify -> True,
  DotSimplify -> True,
  ChiralME -> True,
  NoExpand -> {},
  EditCode -> False,
  RetainFile -> False }

OneLoop[___Rule] = OneLoop[_[], ___Rule] = 0

OneLoop[amps_, opt___Rule] :=
Block[ {lev, diracsimp, dotsimp, chime, noexp, edcode, retain,
hh, amplist, res, res3, res4, smalls,
Hide, Unhide, hidec = 0, extferm = "0", inssum = 0, ins = {},
allins = {}, iabbr = {}, names = {}, alltrace = {}},

  {lev, diracsimp, dotsimp, chime, noexp, edcode, retain} =
    {AmplitudeLevel, DiracSimplify, DotSimplify,
      ChiralME, NoExpand, EditCode, RetainFile} /.
      {opt} /. Options[OneLoop];

  res = Flatten[{amps}] /. {
    2^(1/2) -> r2, 2^(-1/2) -> 1/r2,
    Complex[a_, b_] -> a + "i_" b,
    Eps -> "e_" };

  noexp = Alternatives@@ Flatten[{noexp}];
  If[ Length[noexp] =!= 0,
    Hide[p__] := Plus[p] /; FreeQ[{p}, noexp] || !FreeQ[{p}, kinobjs];
    Hide[a:-_, b__] := -Hide@@ (-#&)/@ {a, b};
    Hide[p__] := Hide[p] =
      ( Unhide[++hidec] = Plus[p] /. "Unhide" -> Unhide;
        "Unhide"[hidec] );
    res = res /. Plus -> Hide;
  ];

  If[ FreeQ[res, Insertions], res = DoTrivialSums/@ res,
    If[lev === Automatic,
      lev = Union[ Cases[res,
        Insertions[t_][__] -> t, Infinity] ][[-1]] ];
    res = DoTrivialSums/@ PickLevel[lev][res];
    allins = Union[Flatten[TakeIns/@ res]];
    ins = Select[allins, Head[#] =!= Symbol &];
    iabbr = Take[IList, Length[ins]];
    ins = Thread[ins -> iabbr];
  ];
  res3 = Select[res, Length[#] === 3 &];
  res4 = Select[res, Length[#] === 4 &];

  OpenFormTemp;

  If[ Head[$Dimension] === Symbol,
    WriteString[hh, "s ", $Dimension, ";\n"] ];
  WriteString[hh, "d ", $Dimension, ";\n"];

  DeclareVars[{#[[3]]&/@ res, allins, SumOver[], Unhide[]}, "li*"];
  WriteString[hh, "\
f Spinor;\n\
#if 'VERSION_' > 1\n\
NT ga, omp, omm, ga5;\n\
#else\n\
f ga, omp, omm, ga5;\n\
#endif\n\n\
.global\n\n\
#procedure Mandelstam()\n" <>
    FormMandelstam <>
    "#endprocedure\n\n#procedure DotSimplify()\n" <>
    If[TrueQ[dotsimp], FormDotSimplify, ""] <>
    "#endprocedure\n\n#procedure DiracSimplify()\n" <>
    If[TrueQ[diracsimp], FormDiracSimplify, ""] <>
    "#endprocedure\n\n"];

  amplist = {};
  If[ Length[res3] =!= 0,
    Scan[FormWrite, res3];
    WriteString[hh, "g amp3 = "];
    Write[hh, Plus@@ amplist, ";"];
    WriteString[hh, ".sort\n\n"];
    inssum = amp3 ];
  Scan[FormWrite, res4];

  Scan[ WriteString[hh, "trace4,", #, ";\n"]&,
    Union[Flatten[alltrace]] ];

  WriteString[hh,
    "\n#define VADecomp \"" <> If[chime =!= False, "0", "1"] <>
    "\"\n#define DIM \"" <> ToString[$Dimension /. D -> 0] <>
    "\"\n#define ExtFermions \"" <> extferm <>
    "\"\n#define SUNStuff \"" <> If[FreeQ[res, SUNobjs], "0", "1"] <>
    "\"\n\n#if 'VERSION_' > 1\n" <>
    "#include " <> $FormCalcDir <> "OneLoop.h2\n#else\n" <>
    "#include " <> $FormCalcDir <> "OneLoop.h1\n" <>
    "#endif\n\n.store\n\n"];
  If[ Length[res4] === 0, WriteString[hh, "l amp = amp3;\n"],
    FormDecl["s ", iabbr];
    WriteString[hh, "\n.global\n"];
    Apply[
      ( WriteString[hh, "\ng ", #1, " = "];
        Write[hh, #2, ";"];
        WriteString[hh, ".store\n"] )&, names, 1 ];
    WriteString[hh, "\nl amp = "];
    Write[hh, inssum, ";"];
    Apply[
      ( WriteString[hh, "id ", #2, " = "];
        Write[hh, #1, ";"] )&, ins, 1 ]
  ];

  WriteString[hh, "\n" <> smalls <> "\n\
#call TrivialSubst{}\n\n\
b SumOver, Mat, DEN, i_, Unhide,\n\
  A0, B0m, B1m, B00m, B11m,\n\
  pave1, pave2, pave3, pave4, pave5;\n\n\
print +s amp;\n\n\
.end\n"];
  Close[hh];

  Block[ {Set},
    (res := Block[#, RunForm])&[Apply[(#2 = #1)&, ins, 1]] ];
  res
]


OpenFormTemp := (
  Print["preparing FORM code in ", $TempFile];
  hh = OpenWrite[
    "!sed -e '/^[^#]/s/\"//g' -e 's/\\[/(/g' -e 's/]/)/g' " <>
    "-e 's/\\*\\*/\\*/g' -e 's/\\*\\\\//g' -e 's/==/=/g' > " <> $TempFile,
    FormatType -> InputForm, PageWidth -> 73];
  WriteString[hh, "#-\n"];
)

RunForm := (
  If[edcode, Pause[1]; Run[StringForm[$Editor, $TempFile]]; Pause[3]];
  WriteString["stdout", "\nrunning FORM... "];
  res = Block[ {Dot = Pair, r2 = Sqrt[2]},
    ReadForm["!" <> $FormCmd <>
      " -s " <> $FormCalcDir <> "form.set " <> $TempFile]
  ] /. FromFormRules /. Pair -> PairAbbr;
  Print["ok"];
  If[!retain, DeleteFile[$TempFile]];
  res
)


(* Things to do when amplitude comes back from FORM *)

	(* for maximum simplification set e.g. o1 = o2 = Simplify *)
o1 = o2 = Identity

FormCalc`d$ = MetricTensor

e$ = EpsAbbr

gi$[_] = Sequence[]

ncm[x_] := x;
ncm[x__] := NonCommutativeMultiply[x]


pave4[i__Integer, p1_, p2_, p3_, m1_, m2_, m3_, m4_] := D0i[
  ToSymbol["dd", Sort[{i}]],
  MomSquare[p1], MomSquare[p1 - p2], MomSquare[p2 - p3],
  MomSquare[p3], MomSquare[p2], MomSquare[p1 - p3],
  m1, m2, m3, m4 ]

pave3[i__Integer, p1_, p2_, m1_, m2_, m3_] := C0i[
  ToSymbol["cc", Sort[{i}]],
  MomSquare[p1], MomSquare[p1 - p2], MomSquare[p2],
  m1, m2, m3 ]

pave2[i__Integer, p_, m1_, m2_] := B0i[
  ToSymbol["bb", Sort[{i}]], MomSquare[p], m1, m2 ]

pave1[i__Integer, m_] := A0i[
  ToSymbol["aa", Sort[{i}]], m ]

B0m[p_, m1_, m2_] := B0[MomSquare[p], m1, m2];
B1m[p_, m1_, m2_] := B1[MomSquare[p], m1, m2];
B00m[p_, m1_, m2_] := B00[MomSquare[p], m1, m2];
B11m[p_, m1_, m2_] := B11[MomSquare[p], m1, m2]

B0[p_, m1_, m2_] := B0[p, m2, m1] /; !OrderedQ[{m1, m2}]

Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B11] = DB11;
Derivative[1, 0, 0][B00] = DB00


a0[0] = A0[0] = 0


(* Abbreviationing business *)

FromForm[x_] := Block[{Dot = PairAbbr}, ToExpression[x]]
(* preliminary def, will be overwritten by FeynAmpList *)


abb[x_?AtomQ] = x

abb[x_] := abb[x] = Unique["Abb"]

abbsum[x_?AtomQ] = x

abbsum[x_] := (abbsum[x] =
Block[ {ft = FactorTerms[x]},
  If[ Head[ft] === Times, MapAt[abbsum, ft, 2], ft ]
]) /; Union[Head/@ List@@ x] === {Times}

abbsum[x_] := abbsum[x] = Unique["AbbSum"]

fme[x_] := fme[x] = Unique["F"]

sun[x_] := sun[x] = Unique["SUN"]

(* mat = Identity *)


Attributes[PairAbbr] = {Orderless}

	(* due to different operator priority in FORM: *)
PairAbbr[a_, x_^n_] := PairAbbr[a, x]^n

PairAbbr[x__] := PairAbbr[x] = Unique["Pair"]

EpsAbbr[x__] := EpsAbbr[x] = Unique["Eps"]


Attributes[Pair] = {Orderless}

Pair[a_, x_Plus] := Block[ {Pair}, Distribute[Pair[a, x]] ]

Pair[a_, n_?NumberQ x_] := n Pair[a, x]

Pair[0, _] = 0


MomSquare[0] = 0

MomSquare[p_] := Pair[p, p]


abdim[x_] := x /; $Scale === 1

abdim[x_] :=
Block[ {c = 0},
  x /. Spinor[__] -> 1 /. k[_] :> ++c;
  x Scale^c
]

dv[sym_, h_:Identity] :=
  Cases[ DownValues[sym], _[ _[_[ab__]], p_Symbol ] :>
    (p -> abdim[h[ab]]) /; FreeQ[{ab}, Pattern] ]

Abbreviations[] :=
Block[ {ab},
  ab = Join[dv[PairAbbr, Pair], dv[EpsAbbr, Eps],
    dv[fme], dv[sun], dv[abb], dv[abbsum]];
  If[ $Scale === 1, ab, Prepend[ab, Scale -> $Scale^-1] ]
]


UseAbbreviations::mismatch =
"Cannot use these abbreviations because `1` were already defined
differently."

UseAbbreviations::scaleinfo =
"$Scale has been set to `1` to use these abbreviations."

UseAbbreviations::badscale =
"These abbreviations were introduced with $Scale = `1` and cannot be
used with the current $Scale = `2`."

UseAbbreviations[file_String] := UseAbbreviations[Get[file]]

UseAbbreviations[abbr_] :=
Block[ {sc, cur = Abbreviations[]},
  sc = 1/Scale /. abbr /. Scale -> 1;
  If[ sc =!= $Scale,
    If[ FreeQ[cur /. Spinor[__] -> 1, k[_]],
      Message[UseAbbreviations::scaleinfo, sc];
      $Scale = sc,
    (* else *)
      Message[UseAbbreviations::badscale, sc, $Scale];
      Return[$Failed] ]
  ];
  sc = Intersection[cur, abbr, SameTest -> (#1[[1]] === #2[[1]] &)];
  If[ Length[sc] =!= 0,
    Message[UseAbbreviations::mismatch, First/@ sc];
    Return[$Failed]
  ];
  newab/@ (abbr /. Scale -> 1)
]

newab[ s_Symbol -> Pair[a__] ] := PairAbbr[a] = s;
newab[ s_Symbol -> Eps[a__] ] := EpsAbbr[a] = s;
newab[ s_Symbol -> p_Plus ] := abbsum[p] = s;
newab[ s_Symbol -> f_?FermionicQ ] := fme[f] = s;
newab[ s_Symbol -> m_ ] := abb[m] = s;
newab[ _ ] = Sequence[]


IntSec[a__Plus] :=
Block[ {i = Intersection[a]},
  If[Head[i] === Plus, i, 0]
]

IntSec[___] = 0

TempInsert[li_, t_] := Insert[li, t,
  Max[0, Position[First/@ li, #, {1}, 1]&/@
    DeleteCases[Level[t[[2]], {-1}], -1]] + 1]

OptLevel12[rul_] :=
Block[ {l = 0, pl, rl = Length[rul], i, is, iss, isl, pr,
repl, nurul = {}},
  Print["level 1: redundancy removal"];
  Attributes[pl] = {Flat, Orderless};
  Apply[
    ( WriteString["stdout", "\r", --rl, " "];
      Set@@ {pl@@ -#2, -#1};
      Set@@ {pl@@ #2, #1} )&,
    Sort[rul, Length[ #1[[2]] ] < Length[ #2[[2]] ] &], 1 ];
  repl = Cases[DownValues[pl],
    (_[_[p__]] :> o_Symbol) :> (o -> Plus[p])];

  Print["\rlevel 2: common subexpression elimination"];
  rl = Length[repl] - 1;
  Do[
    pr = False;
    While[
      is = IntSec[ repl[[i, 2]], repl[[i + 1, 2]] ];
      If[ Length[is] < Length[ repl[[i, 2]] ]/2,
        iss = IntSec[ repl[[i, 2]], -repl[[i + 1, 2]] ];
        If[Length[is] < Length[iss], is = iss]
      ];
      Length[is] > 3,
    (* while body: *)
      isl = Ceiling[Length[is]/2];
      iss = IntSec@@ Select[
        (IntSec[#[[2]], is] + IntSec[-#[[2]], is])&/@
          Drop[repl, {i, i + 1}],
        Length[#] > isl &];
      If[Length[iss] < 4, iss = is];
      is = ToSymbol["tmp", ++l];
      AppendTo[nurul, is -> iss];
      repl = repl /. iss -> is /. -iss -> -is;
      pr = True;
    ];
    If[pr, WriteString["stdout", "\r", rl - i, " "]],
  {i, rl}];
  WriteString["stdout", "\r     \r"];
  Fold[TempInsert, repl, nurul]
]

OptimizeAbbreviations[file_String] :=
Block[ {abbr = Get[file], new, orig = file <> ".orig"},
  new = OptimizeAbbreviations[abbr];
  If[ FileType[orig] === File, DeleteFile[orig] ];
  RenameFile[file, orig];
  Put[new, file];
]

OptimizeAbbreviations[rul_] :=
  Join[ Select[rul, FreeQ[#, Plus]&],
    OptLevel12[Select[rul, !FreeQ[#, Plus]&]] ]


(* This is for backward compatibility with FeynCalc *)

C0[p1_, p2_, p1p2_, m1_, m2_, m3_] :=
  C0i[cc0, p1, p2, p1p2, m1, m2, m3]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_}] :=
  C0i[ToSymbol["cc", Sort[{i}]], p, m1, m2, m3]

D0[p1_, p2_, p3_, p4_, p1p2_, p2p3_, m1_, m2_, m3_, m4_] :=
  D0i[dd0, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_, m4_}] :=
  D0i[ToSymbol["dd", Sort[{i}]], p, m1, m2, m3, m4]


FeynCalcGet[mask___] :=
Block[ {Global`OneLoopResult, GraphName},
  GraphName[__] = 0;
  Plus@@ ((Print[#]; Get[#]; Global`OneLoopResult[0])&)/@
    FileNames[mask] /.
    (p:S | T | U - m_)^(n_?Negative) :> DEN[p, m]^(-n) /.
    (p:S | T | U)^(n_?Negative) :> DEN[p, 0]^(-n) /.
    (m_ - p:S | T | U)^(n_?Negative) :> -DEN[p, m]^(-n) /.
    ep_Eps :> I ep
]


FeynCalcPut[expr_, file_] :=
Block[ {PaVe, C0i, D0i, C0, D0},
  C0i[Global`cc0, args__] := C0[args];
  D0i[Global`dd0, args__] := D0[args];
  C0i[i_, p1_, p2_, p1p2_, m1_, m2_, m3_] := PaVe[
    Sequence@@ (ToExpression/@ Characters[StringDrop[ToString[i], 2]]),
    {p1, p2, p1p2}, {m1, m2, m3}];
  D0i[i_, p1_, p2_, p3_, p4_, p1p2_, p2p3_, m1_, m2_, m3_, m4_] := PaVe[
    Sequence@@ (ToExpression/@ Characters[StringDrop[ToString[i], 2]]),
    {p1, p2, p3, p4, p1p2, p2p3}, {m1, m2, m3, m4}];
  Put[expr, file]
]


(* Classifying amplitudes *)

FermionFamily[-x_] := FermionFamiliy[x]

FermionFamily[_] = ""


ClassIndex[ FeynAmp[_, top_, amp_, ___] ] :=
Block[ {cname = ""},
  If[howcls =!= IndexSumsOnly,
    If[FermionicQ[amp], cname = "F"];
    If[howcls === Tough && !FreeQ[amp, d0],
      If[cname === "F", cname = cname <> Union[Cases[
        Cases[amp, (DiracTrace | SpinorChain | Dot)[t__] :> t, Infinity],
        _. _DiracSlash + m_. :> FermionFamily[m] ]] ];
      cname = cname <> ToString[top]
    ]
  ];
  cname
]


SeparateSums[amp_Plus, cname__] :=
  { SeparateTerm[Select[amp, FreeQ[#, SumOver]&], cname],
    (SeparateTerm[#, cname]&)/@
      Select[List@@ amp, !FreeQ[#, SumOver]&] } /; !FreeQ[amp, SumOver]

SeparateSums[amp___] := SeparateTerm[amp]

SeparateTerm[amp_, outtag_, cname_] :=
Block[ {na},
  na = outtag <> cname <>
    Cases[amp, SumOver[i_, ___] :> "_" <> ToString[i]] <> ".m";
  Print["writing ", na];
  Put[amp, na];
  na
]

SeparateTerm[amp_, cname_] :=
Block[ {s},
  s = Cases[amp, SumOver[i_, ___] -> i];
  Print["part ", ++modnr,
    If[Length[s] === 0, " (no index sums)",
      " (sum over " <> StringTake[ToString[s], {2, -2}] <> ")"]
  ];
  amp
]


Options[ProcessFile] = {Classification -> Standard}

ProcessFile[filename_String, outtag___String, opt___Rule] :=
  ProcessFile[Get[filename], outtag, opt]

ProcessFile[amps_List, outtag___String, opt___Rule] :=
Block[ {howcls, class, theamps, res, modnr = 0},
  howcls = Classification /. {opt} /. Options[ProcessFile];
  theamps = If[ FreeQ[amps, Insertions], amps,
    lev = AmplitudeLevel /. {opt} /. Options[OneLoop];
    If[lev === Automatic,
      lev = Union[ Cases[amps,
        Insertions[t_][__] -> t, Infinity] ][[-1]] ];
    PickLevel[lev][amps] ];
  class[_] = {};
  Scan[AppendTo[class[ClassIndex[#]], #]&, theamps];
  Cases[ DownValues[class], (_[_[cname_String]] :> amp_) :>
    ( Print[""];
      Print["processing class [", cname, "]"];
      res = OneLoop[amp, opt];
      Share[];
      SeparateSums[res, outtag, cname] ) ]//Flatten
]


(* Helicity matrix elements *)

rho[fline_, k[n_], 0, sign_] :=
  "g_"[fline, k[n]] ** ("gi_"[fline] + sign Hel[n] "g5_"[fline])/2

rho[fline_, k[n_], m_, sign_] :=
  ("g_"[fline, k[n]] + sign m "gi_"[fline]) **
    ("gi_"[fline] + Hel[n] "g5_"[fline] ** "g_"[fline, e[n]])/2

rhoc[fline_, k[n_], 0, sign_] :=
  -("gi_"[fline] + sign Hel[n] "g5_"[fline]) ** "g_"[fline, k[n]]/2

rhoc[fline_, k[n_], m_, sign_] :=
  ("gi_"[fline] - Hel[n] "g_"[fline, e[n]] ** "g5_"[fline]) **
    (-"g_"[fline, k[n]] + sign m "gi_"[fline])/2


ToTrace[fi_ -> plain_, fj_ -> conj_] :=
Block[ {me, fline = 0},
  me = plain conj //. {
    (a___ ** Spinor[km__]) (Spinor[km__] ** b___) :> a ** rho[km] ** b,
	(* If the spinors at the ends don't match directly, we
	   have to reverse one chain. This is a charge conjugation,
	   not a hermitian conjugation like the one HelicityME
	   does for the "conj" part. The rules:
	     a) reverse the chain and exchange u <-> v
	     b) gamma_mu -> -gamma_mu
	     c) add a global minus sign to compensate for the change
	        in the permutation of the external fermions.
	   For details see the Denner/Eck/Hahn/Kueblbeck paper. *)
    (Spinor[k1_, m1_, _] ** a___ ** Spinor[k2_, m2_, s2_]) *
      (Spinor[k1_, m1_, s1_] ** b___) :>
      -(1 - 2 Mod[Count[{a}, _ga], 2]) *
        (Reverse[a ** Spinor[k2, m2, -s2]] /.
          {rho -> rhoc, rhoc -> rho}) ** rho[k1, m1, s1] ** b
  } /.
    Spinor[km__] ** a___ ** Spinor[km__] :> rho[km] ** a /.
    n_NonCommutativeMultiply :>
      ( ++fline;
        n /. { ga[x_] :> "g_"[fline, x],
               ga5 :> "g5_"[fline],
               omp :> "g6_"[fline]/2,
               omm :> "g7_"[fline]/2,
               (r:rho | rhoc)[x__] :> r[fline, x] } );
  Mat[fi, fj] -> tr[fline, me]
]


ToHel[k[n_], __] := {{}, {}} /; Head[Hel[n]] =!= Hel

ToHel[k[n_], 0, s_] :=
Block[ {h = ToSymbol["Hel", n]},
  {Hel[n] -> h - s, h -> Hel[n] + s}
]

ToHel[k_, __] := ToHel[k, 0, 0]


SelectArg[All] := fabbr

SelectArg[expr_] := Union[Select[fabbr, !FreeQ[ expr, #[[1]] ]&]]


Options[HelicityME] = {
  AbbreviationsToUse :> Abbreviations[],
  EditCode -> False,
  RetainFile -> False }

HelicityME[plain_, conj_, opt___] :=
Block[ {abbr, fabbr, edcode, retain, part, tohel, fromhel, hels,
hh, res, e, traces, smalls},
  Update[Spinor];
  {abbr, edcode, retain} = {AbbreviationsToUse, EditCode, RetainFile} /.
    {opt} /. Options[HelicityME];

  fabbr = Select[abbr, !FreeQ[#, Spinor]&];
  abbr = SelectArg/@ ({plain, conj} //. Select[abbr, FreeQ[#, Spinor]&]);

  part = Cases[abbr, _Spinor, Infinity]//Union;
  {tohel, fromhel} = Flatten/@ Transpose[Apply[ToHel, part, 1]];
  part = #[[1, 1]]&/@ part;
  hels = First/@ fromhel;

  traces = Outer[ ToTrace, abbr[[1]],
    abbr[[2]] /.
      {omp -> omm, omm -> omp, ga5 -> -ga5, ep_Eps -> -ep} /.
      n_NonCommutativeMultiply :> Reverse[n] /.
      { IND1 -> IND1c, IND2 -> IND2c,
        IND3 -> IND3c, IND4 -> IND4c } ]//Flatten;
  traces = traces /. tohel /. Reverse/@ FromFormRules /. Eps -> "e_";

  OpenFormTemp;
  DeclareVars[{Last/@ traces, hels}, "IND*"];

  WriteString[hh, "\n\
i J1, J2, J3, J4;\n\
v P1, P2;\n\
s X1, X2;\n\
cf abb;\n\n\
.global\n\n\
#procedure Simplify()\ncontract,0;\n" <>
    FormMandelstam <>
    smalls <>
    FormDotSimplify <> "\
id P1? . P2? = abb(P1 . P2);\n\
id e_(J1?, J2?, J3?, J4?) = abb(e_(J1, J2, J3, J4));\n\
repeat;\n\
  id abb(X1?) * abb(X2?) = abb(X1 * X2);\n\
endrepeat;\n" ];
  FormDecl["b ", hels];
  WriteString[hh, "print +s;\n.store\n#endprocedure\n"];    

  Apply[
    ( WriteString["stdout", x = "Mat" <> ToString/@ List@@ #1 <> " "];
      WriteString[hh, "\nl " <> x <> "= "];
      Write[hh, #2[[2]], ";"];
      Array[ WriteString[hh, "trace4,", #, ";\n"]&, #2[[1]] ];
      WriteString[hh, "#call Simplify{}\n"] )&,
    traces, 1 ];

  WriteString[hh, "\n.end\n"];
  Close[hh];

  (e[#] = s[#])&/@ part;

  Thread[First/@ traces -> RunForm] /. fromhel
]


ToMat[m1_Symbol, m2_Symbol] := Mat[m1, m2]

ToMat[m1_, m2_] := Inner[Mat, m1, m2, Times]

ToSquared[Mat[m1_] x1_., Mat[m2_] x2_.] :=
  x1 Conjugate[x2] ToMat[m1, m2]


SquaredME[plain_, conj_] :=
  plain Conjugate[conj] /; FreeQ[{plain, conj}, Mat]

SquaredME[plain_, conj_] :=
  Plus@@ Flatten[Outer[ToSquared,
    Flatten[{plain} //. m_Mat x_. + r_ :> {m x, r}],
    Flatten[{conj} //. m_Mat x_. + r_ :> {m x, r}]
  ]]
    
SquaredME[amp_] := SquaredME[amp, amp]


Unprotect[Conjugate];
Conjugate[p:_Plus | _DEN | _Pair] := Conjugate/@ p;
Conjugate[sym_?RealQ] := sym;
Protect[Conjugate]

(RealQ[#] = True)&/@ {S, T, U, Sf, Tf, Uf}


End[]


Format[ Continuation[_] ] = "    "
  (* eliminate those `>' in front of continuation lines so one can cut
     and paste more easily *)

$OnShell = True
  (* whether the external particles are on-shell *)

$Transversality = True
  (* whether to enforce transversality on the external polarization
     vectors, i.e. e[i].k[i] = 0 *)

$Dimension = D

$FormCmd = "form"
  (* the filename of the actual FORM executable; may contain a path *)

$Editor = "xterm -geometry 80x30 -e pico `1` &"
  (* which editor to use when debugging FORM code *)

$Scale = 1
  (* factor which scales momenta in abbreviations, see usage *)

$TempFile = "t" <> ToString[$ProcessID] <> ".frm"

EndPackage[]


(* global definitions for specific models *)

(* definitions for the Standard Model *)

EL/: EL^4 = 16 Pi^2 a2

GS/: GS^4 = 16 Pi^2 as2

MW^(n_?EvenQ) ^= MW2^(n/2);
MZ^(n_?EvenQ) ^= MZ2^(n/2);
MH^(n_?EvenQ) ^= MH2^(n/2)

ME^(n_?EvenQ) ^= ME2^(n/2);
MM^(n_?EvenQ) ^= MM2^(n/2);
ML^(n_?EvenQ) ^= ML2^(n/2);
MLE[a__]^(n_?EvenQ) ^= MLE2[a]^(n/2)

MU^(n_?EvenQ) ^= MU2^(n/2);
MC^(n_?EvenQ) ^= MC2^(n/2);
MT^(n_?EvenQ) ^= MT2^(n/2);
MQU[a__]^(n_?EvenQ) ^= MQU2[a]^(n/2)

MD^(n_?EvenQ) ^= MD2^(n/2);
MS^(n_?EvenQ) ^= MS2^(n/2);
MB^(n_?EvenQ) ^= MB2^(n/2);
MQD[a__]^(n_?EvenQ) ^= MQD2[a]^(n/2)

(* these are the identifiers for fermion classes used by ProcessFile *)

FermionFamily[0] = "N";				(* neutrinos *)
FermionFamily[ME | MM | ML | MLE[_]] = "E";	(* massive leptons *)
FermionFamily[MU | MC | MT | MQU[_]] = "U";	(* up-type quarks *)
FermionFamily[MD | MS | MB | MQD[_]] = "D"	(* down-type quarks *)

(* these symbols represent real quantities, i.e. Conjugate[sym] = sym
   for any of these. Thinking e.g. of complex masses this looks
   dangerous but then again it's easy to remove any such definition.
   The function that really needs this is SquaredME. *)

(RealQ[#] = True)&/@
  { EL, a2, GS, as2, SW, CW, S2, C2, S4, C4,
    MW, MW2, MZ, MZ2, MH, MH2,
    ME, ME2, MM, MM2, ML, ML2, _MLE, _MLE2,
    MU, MU2, MC, MC2, MT, MT2, _MQU, _MQU2,
    MD, MD2, MS, MS2, MB, MB2, _MQD, _MQD2 }


(* definitions for the MSSM *)

SetOptions[OneLoop,
  NoExpand -> {USf, USfC, UCha, UChaC, VCha, VChaC, ZNeu, ZNeuC}]

USf[t_, g_][a_, b_] := USf[a, b, t, g]

Conjugate[USf[a__]] ^:= USfC[a];
Conjugate[UCha[a__]] ^:= UChaC[a];
Conjugate[VCha[a__]] ^:= VChaC[a];
Conjugate[ZNeu[a__]] ^:= ZNeuC[a]

Conjugate[USfC[a__]] ^:= USf[a];
Conjugate[UChaC[a__]] ^:= UCha[a];
Conjugate[VChaC[a__]] ^:= VCha[a];
Conjugate[ZNeuC[a__]] ^:= ZNeu[a]

MSf[a__]^(n_?EvenQ) ^= MSf2[a]^(n/2);
MCha[a__]^(n_?EvenQ) ^= MCha2[a]^(n/2);
MNeu[a__]^(n_?EvenQ) ^= MNeu2[a]^(n/2)

Mh0^(n_?EvenQ) ^= Mh02^(n/2);
MHH^(n_?EvenQ) ^= MHH2^(n/2);
MA0^(n_?EvenQ) ^= MA02^(n/2);
MG0^(n_?EvenQ) ^= MG02^(n/2);
MHp^(n_?EvenQ) ^= MHp2^(n/2);
MGp^(n_?EvenQ) ^= MGp2^(n/2)

(RealQ[#] = True)&/@
  { TB, CB, SB, CA, SA, C2A, S2A, CAB, SAB, CBA, SBA,  
    Mh0, Mh02, MHH, MHH2, MA0, MA02, MG0, MG02,
    MHp, MHp2, MGp, MGp2,
    _MSf, _MSf2, _MSNE, _MSLE1, _MSLE2,
    _MSQU1, _MSQU2, _MSQD1, _MSQD2,
    _MCha, _MCha2, _MNeu, _MNeu2 }

Null

