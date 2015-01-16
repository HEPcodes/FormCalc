(*

This is FormCalc, Version 1.4
Copyright by Thomas Hahn 1999
last modified 8 Jul 99 by Thomas Hahn

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
Print["FormCalc 1.4"];
Print["by Thomas Hahn"];
Print["last revision: 8 Jul 99"];


BeginPackage["FormCalc`"]

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

Spinor::usage = "Spinor[p, m] is the spinor with momentum p and mass m.
Spinor corresponds to the more conventional way of writing spinors by\n
   Spinor[p, m, 1] ** ...  -> \\bar u\n
   Spinor[p, m, -1] ** ... -> \\bar v\n
   ... ** Spinor[p, m, 1]  -> u\n
   ... ** Spinor[p, m, -1] -> v."

LeptonSpinor = Spinor

QuarkSpinor = Spinor

DiracSlash::usage = "DiracSlash[p] represents p_mu gamma_mu."

DiracMatrix::usage = "DiracMatrix[mu] represents the Dirac matrix with
Lorentz index mu."

ChiralityProjector::usage = "ChiralityProjector[+-1] represents the
chirality projectors omega_{+-} = (1 +- gamma_5)/2."

MetricTensor::usage = "MetricTensor[mu, nu] represents the metric tensor
with Lorentz indices mu and nu."

Index::usage = "Index[t, n] represents an index of type t with number n."

IndexDelta::usage = "IndexDelta[i1, i2] is used to force indices to be the
same. For integers i1 and i2, IndexDelta is 1 if i1 = i2 or 0 otherwise.
If i1 is a symbol, i2 in the expression multiplied by IndexDelta is
replaced by i1 and vice versa."

SumOver::usage = "SumOver[i, r] specifies that the amplitude which it is
multiplied with is to be summed in the index i over the range r."

SpinorChain::usage = "SpinorChain[Spinor[...], ..., Spinor[...]]
represents an open fermion chain."

DiracTrace::usage = "DiracTrace represents the trace over Dirac matrices
in a fermion loop."

ga::usage = "ga[li] represents the gamma matrix with Lorentz index li."

ga5::usage = "ga5 represents gamma_5."

omp::usage = "omp represents the right handed chirality projector."

omm::usage = "omm represents the left handed chirality projector."

LI::usage = "LI[n] is the nth Lorentz index in spinor chains. Summation
over all LI is implied."

FourVector::usage = "FourVector[p, mu] represents the four-vector p with
Lorentz index mu."

PolarizationVector::usage = "PolarizationVector[p, mu] represents the
polarization vector belonging to the momentum p with Lorentz index mu."

Momentum::usage = "Momentum[p] is an obsolete form to represent the
four-momentum p."

q1::usage = "q1 is the integration momentum in a one-loop amplitude."

ScalarProduct::usage = "ScalarProduct[p1, p2] is the scalar product of two
four-vectors p1 and p2. It is converted to Pair[p1, p2] in FormCalc."

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

DEN::usage = "DEN[p2, m2] stands for 1/(p2 - m2). Note that in contrast
to PropagatorDenominator p2 and m2 are the momentum and mass _squared_."

GaugeXi::usage = "GaugeXi[v] is the way FeynArts denotes the gauge
parameter for the gauge boson v. It is converted to xi[v] in FormCalc."

xi::usage = "xi[v] is the gauge parameter of gauge boson v."

OneLoop::usage = "OneLoop[amps] is the basic function for evaluating
Feynman diagrams. It takes as an argument either a single diagram (with
head FeynAmp) or a list of diagrams and returns the results in a collected
and abbreviated form."

ProcessFile::usage = "ProcessFile[amps, outtag] evaluates the diagrams in
amps by running them through OneLoop and splits the result into a bosonic
part which is written to outtag.m and a fermionic part which is written
to outtagF.m. outtag may contain a pathname together with an identifier
for the diagrams, e.g. destpath/box for box diagrams."

Comment::usage = "Comment -> \"string\" is an option of OneLoop and places
an extra comment at the beginning of the temporary file handed to FORM.
This may be useful when debugging FORM code."

AmplitudeLevel::usage = "AmplitudeLevel is an option of OneLoop. It is
used in amplitudes generated with FeynArts 2 only, and specifies the level
(Classes or Particles) at which to calculate the amplitudes. The default
setting Automatic takes the deepest level available."

DiracSimplify::usage = "DiracSimplify is an option of OneLoop. With
DiracSimplify -> True (the default), OneLoop simplifies fermionic matrix
elements by replacing in turn all momenta via momentum conservation and
then using the Dirac equation. However, this procedure can be rather slow,
so DiracSimplify -> False can be used to turn it off."

CancelEps::usage = "CancelEps is an option of OneLoop. It is used for CP
invariant processes which must not contain terms proportional to Eps.
Setting CancelEps -> True will cancel these terms at an early stage thus
speeding up the calculation (they should drop out in the end in any case).
The default is CancelEps -> False."

ChiralME::usage = "ChiralME is an option of OneLoop. It specifies whether
spinor chains are returned in terms of the chirality projectors omega_+
and omega_-, or in the vector/axial vector decomposition 1 and gamma_5.
The former is the default since it yields much more compact helicity
matrix elements."

NoExpand::usage = "NoExpand is an option of OneLoop. NoExpand -> {sym1,
sym2, ...} specifies that sums containing any of sym1, sym2, ... are
expanded in FORM."

EditCode::usage = "EditCode is an option of OneLoop and HelicityME. It
edits the temporary file passed to FORM using $Editor and is of course
used only for debugging."

RetainFile::usage = "RetainFile is an option of OneLoop and HelicityME.
When set to True, it prevents removing the temporary file which contains
the FORM input after running FORM."

DotSimplify::usage = "DotSimplify is an option of OneLoop. When set to
True, FormCalc tries to reduce the number of independent matrix elements
by eliminating one momentum via momentum conservation in dot products.
Empirically, bosonic amplitudes are shortened but fermionic ones are not.
If getting the shortest possible result is your aim, you might want to set
DotSimplify -> Automatic (the default) which is equivalent to True for
bosonic and False for fermionic amplitudes."

Classification::usage = "Classification is an option of ProcessFile.
It can take the values IndexSumsOnly, Standard (default), and Tough which
control how ProcessFile divides the diagrams into different classes.
In all cases, diagrams with different index summations are separated,
i.e. each class contains only diagrams which have the same index
summations. For Standard, the amplitude is split into bosonic and
fermionic parts, too, and for Tough also into different topologies and (if
applicable) different fermion families. This is mainly used for large
amplitudes."

IndexSumsOnly::usage = "IndexSumsOnly is a value the Classification option
of ProcessFile can take. The criterium for two diagrams to be put into
different classes is in this case that they have different index
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

ReadForm::usage = "ReadForm[\"file\"] reads file which it assumes to
contain FORM output into Mathematica.\n
ReadForm[\"!cmd\"] executes cmd and pipes its output back into
Mathematica."

PowerCountingFor::usage = "PowerCountingFor[momlist] tells ReadForm to
multiply each term in the FORM result by $O1ME^n, where n is the number
of momlist members appearing in the term. Usually, momlist is the list
of momenta, e.g. {p1, p2, k1, k2}. This function is intended for internal 
use only."

r2::usage = "r2 is a symbol used to substitute Sqrt[2] before handing an
expression to FORM. It is used internally only."

Abbreviations::usage = "Abbreviations[] returns a list of all
abbreviations introduced so far. It is typically used at the end of a
calculation to save the abbreviations with a command like
Abbreviations[] >> abbrfile."

UseAbbreviations::usage = "UseAbbreviations[abbr] registers the
abbreviations abbr with FormCalc. abbr may be either the abbreviations
themselves (a list of rules) or a filename containing them. This is useful
when continuing a former FormCalc session. UseAbbreviations requires that
the value of $O1ME that was used in the former session is also set in the
current session.";

OptimizeAbbreviations::usage = "OptimizeAbbreviations[abbr] optimizes the
set of abbreviations returned by Abbreviations[] by eliminating common
subexpressions."

Scale::usage = "Scale is a scale introduced via $O1ME to scale matrix
elements, e.g. to make them dimensionless. The variable Scale itself
appears in the Abbreviations[]."

Pair::usage = "Pair[a, b] represents the Minkovskian scalar product of
the four-vectors a and b."

Eps::usage = "Eps[a, b, c, d] represents -I times the total antisymmetric
Levi-Civita tensor: Eps[a, b, c, d] = -I a[mu] b[nu] c[rho] d[sigma]
epsilon[mu, nu, rho, sigma], where a, b, c, d are four-vectors and the
sign convention is epsilon[0, 1, 2, 3] = +1."

e::usage = "e[n] is the nth polarization vector."

k::usage = "k[n] is the nth momentum."

s::usage = "s[n] is the nth helicity reference vector."

DiagramType::usage = "DiagramType[diag] returns 2 for a self-energy
diagram, 1 for a vertex diagram, and 0 for a box. In other words, it
returns the number of denominators not containing the integration
momentum."

FermionicQ::usage = "FermionicQ[diag] gives True for a diagram containing
fermions and False otherwise."

Small::usage = "Small[sym] = 0 tells FormCalc to put sym = 0 except when
it appears in negative powers or in loop integrals."

HelicityME::usage = "HelicityME[M1, M2] returns a list of matrix
elements of the expression M1 M2^*, where typically M1 is the one-loop
result, and M2 the Born expression. The matrix elements returned contain
the helicity information as Hel[1], Hel[2], ... for each external
particle, where each Hel[i] can take the values +-1. If M1 or M2 are
amplitudes, only the fermionic matrix elements appearing in these
amplitudes are calculated. Alternatively, M1 or M2 may be set to All,
in which case all fermionic matrix elements currently defined in
Abbreviations[] are used."

AbbreviationsToUse::usage = "AbbreviationsToUse is an option of
HelicityME. It specifies which abbreviations are used to calculate the
helicity matrix elements. For instance, AbbreviationsToUse :> <<file
takes abbreviations stored in a file. Or, if the helicity matrix elements
shall be calculated without the axial part, AbbreviationsToUse :>
(Abbreviations[] /. omp | omm -> 1) will do the trick."

All::usage = "All is a possible input value for HelicityME, indicating
that all fermionic matrix element currently defined in Abbreviations[]
should be used instead of just those appearing in a particular
expression."

Hel::usage = "Hel[i] denotes the helicity of the ith external particle. It
can take the values +-1."

Mat::usage = "Mat[i, j] is the symbol FormCalc produces to denote the
helicity matrix element composed from the expression M1 M2^* of the ith
matrix element of M1 and the jth matrix element of M2."

SquaredME::usage = "SquaredME[plain, conj] returns the matrix element
plain Conjugate[conj]. This performs a nontrivial task only for fermionic
amplitudes: the product of two fermionic amplitudes\n
    M1 = a1 F1 + a2 F2 + ... and\n
    M2 = b1 F1 + b2 F2 + ... is returned as\n
    M1 M2^* = a1 b1^* Mat[F1, F1] + a2 b1^* Mat[F2, F1] + ...\n
The special case of plain === conj can be written as SquaredME[plain]
which is of course equivalent to SquaredME[plain, plain]."

RealQ::usage = "RealQ[sym] is True if sym represents a real quantity which
means in particular that Conjugate[sym] = sym."

Pick::usage = "Pick[amp, no] picks diagrams from amp where no is of the
form {5} which picks diagram # 5, {{5, 9}} which picks diagrams # 5
through 9, or All which picks all diagrams. Combinations are also allowed,
for example, Pick[amp, {6, {10, 12}, 27}] picks diagrams # 6, 10, 11, 12,
and 27 from amp."

PickLevel::usage = "PickLevel[lev][amps] selects amplitude level lev from
amps."

FeynCalcGet::usage = "FeynCalcGet[mask] reads files produced with
FeynCalc. mask is taken as input to the Mathematica function FileNames, so
it might be FeynCalcGet[\"file.m\"] or FeynCalcGet[\"*.m\"] or
FeynCalcGet[\"*.m\",\"~/feyncalcfiles\"]."

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

a0::usage = "a0[m] is used internally to represent a pre-form of the
one-point function."

A0::usage = "A0[m] is the one-point scalar Passarino-Veltman function
where m is the mass squared."

b0::usage = "b0[m, p] is used internally to represent a pre-form of the
two-point function."

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

c0::usage = "c0[p1, p2, m1, m2, m3] is used internally to represent a
pre-form of the three-point function."

C0i::usage = "C0i[id, p1, p2, p1p2, m1, m2, m3] is the generic three-point
Passarino-Veltman function which includes both scalar and tensor integrals
specified by id. For example, C0i[cc0, ...] is the scalar function C0,
C0i[cc112, ...] the tensor function C_112 etc. Call the external momenta
k1...k3, then the arguments are given as p1 = k1^2, p2 = k2^2,
p1p2 = (k1 + k2)^2, and m1...m3 are the masses squared."

d0::usage = "d0[p1, p2, p3, m1, m2, m3, m4] is used internally to
represent a pre-form of the four-point function."

D0i::usage = "D0i[id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4] is the
generic Passarino-Veltman four-point function which includes both scalar
and tensor integrals specified by id. For example, D0i[dd0, ...] is the
scalar function D0, D0i[dd1233, ...] the tensor function D_1233 etc.
Call the external momenta k1...k4, then the arguments are given as
p1 = k1^2, p2 = k2^2, p3 = k3^2, p4 = k4^2, p1p2 = (k1 + k2)^2,
p2p3 = (k2 + k3)^2, and m1...m4 are the masses squared."

e0::usage = "e0[p1, p2, p3, p4, m1, m2, m3, m4, m5] is used internally to
represent a pre-form of the five-point function."

o1::usage = "o1 is a head wrapped around the prefactor of a standard
matrix element. Its default value is o1 = Identity (doing nothing). To
get the shortest possible amplitude, you can set o1 = Simplify or
similar."

o2::usage = "o2 is a head wrapped around a linear combination of standard
matrix elements. Its default value is o2 = Identity (doing nothing). To
get the shortest possible amplitude, you can set o2 = Simplify or
similar."

$Editor::usage = "$Editor specifies the editor used in debugging FORM
code."

$TempFile::usage = "$TempFile is the name for the temporary file used to
write out FORM code in OneLoop. Unless RetainFile is set to True, it is
deleted after usage."

$OnShell::usage = "$OnShell specifies whether the external particles in
your calculation are on shell. Caution: This variable must be set *before*
loading amplitudes."

$FormCalcDir::usage = "$FormCalcDir is the path where FormCalc
(specifically, form.set and OneLoop.h[12]) resides."

$Platform::usage = "$Platform is a string that identifies the platform
FormCalc is running on. It is the value of the environment variable
HOSTTYPE in the Bourne shell (sh) and is used to distinguish the ReadForm
binaries of different platforms."

$FormCmd::usage = "$FormCmd gives the name of the actual FORM executable.
It may contain a path."

$Transversality::usage = "$Transversality specifies whether FormCalc may
assume transversality for polarization vectors, i.e. e[i].k[i] = 0."

$O1ME::usage = "Matrix elements are scaled with $O1ME for every momentum
they include, e.g. e1.k2 is abbreviated by k12/$O1ME, e1.k2 e2.k3 by
k12 k23/$O1ME^2, etc. This way, all matrix elements can be made
dimensionless. For instance, in a 2 -> 2 process, $O1ME = Sqrt[S] makes
all matrix elements of O(1)."

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
  amp[[ If[graphs === All, Range[Length[amp]],
    Flatten[ If[Head[#] === List, Range@@ #, #]&/@ graphs ]//Union] ]]

Attributes[DiagramType] = {Listable};
DiagramType[a_FeynAmp] := Exponent[a[[3]] /. DEN[__] -> DEN, DEN]

FermionicQ[a_] :=
  !FreeQ[a, Spinor | ga | omp | omm | ChiralityProjector |
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

PickLevel[_][a_FeynAmp] := a

PickLevel[lev_][a_] := PickLevel[lev]/@ a


PropagatorDenominator[p_, m_] :=
  PropagatorDenominator[-p, m] /; !FreeQ[p, -q1]

PropagatorDenominator[0, m_] = -1/m^2

FeynAmpDenominator[ p:PropagatorDenominator[q1, _], b___ ] :=
  den[p, b];
FeynAmpDenominator[ a__, p:PropagatorDenominator[q1, _], b___ ] :=
  den[p, b, a]

den[ a___, PropagatorDenominator[q_, m1_], b___,
  PropagatorDenominator[q_, m2_], c___ ] :=
  1/Factor[m1^2 - m2^2] (den[a, PropagatorDenominator[q, m1], b, c] -
    den[a, PropagatorDenominator[q, m2], b, c]);
den[ PropagatorDenominator[q1, m1_] ] :=
       I Pi^2 a0[m1^2];
den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_] ] :=
       I Pi^2 b0[p2 - p1, m1^2, m2^2];
den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_] ] :=
       I Pi^2 c0[p2 - p1, p3 - p1, m1^2, m2^2, m3^2];
den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_],
     PropagatorDenominator[p4_, m4_] ] :=
       I Pi^2 d0[p2 - p1, p3 - p1, p4 - p1, m1^2, m2^2, m3^2, m4^2];
den[ PropagatorDenominator[p1_, m1_],
     PropagatorDenominator[p2_, m2_],
     PropagatorDenominator[p3_, m3_],
     PropagatorDenominator[p4_, m4_],
     PropagatorDenominator[p5_, m5_] ] :=
       I Pi^2 e0[p2 - p1, p3 - p1, p4 - p1, p5 - p1,
                 m1^2, m2^2, m3^2, m4^2, m5^2]


Index[t_, n_] := Index[t, n] = ToSymbol[StringTake[ToString[t], 3], n];


Attributes[IndexDelta] = {Orderless}

IndexDelta[n_, n_] = 1

IndexDelta[_Integer, n_Integer] = 0


(* Fermion stuff: Dirac algebra etc *)

DiracSlash[p_] := -DiracSlash[-p] /; !FreeQ[p, -q1]

Spinor[(s_Integer:1) p_, m_] := Spinor[p, m, s]

SpinorType[-1] = "v";
SpinorType[1] = "u"

Unprotect[NonCommutativeMultiply];
Format[ Spinor[p1_, m1_, s1_] ** c___ ** Spinor[p2_, m2_, s2_] ] :=
  SequenceForm[ "(",
    ColumnForm[{"_", SpinorType[s1][p1, m1]}, Left, Above] **
      c ** SpinorType[s2][p2, m2], ")" ];
Protect[NonCommutativeMultiply]

Unprotect[Dot];
Format[ sp_Spinor . c__ ] := sp ** c
Protect[Dot]


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
  extferm = 1;
  ++fline;
  DiracMatrix[5] = omp[fline]/2 - omm[fline]/2;
  DiracMatrix[li_] = ga[fline, li];
  DiracSlash[k_] :=
    Map[If[Head[#] === Symbol, ga[fline, #], #]&, k, {-1}];
  ChiralityProjector[1] = omp[fline];
  ChiralityProjector[-1] = omm[fline];
  NonCommutativeMultiply[expr]
]


(* Boson stuff: 4-vectors etc *)

GaugeXi = xi

FourVector[ksum_Plus, li_] := FourVector[#, li]&/@ ksum

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
  PolarizationVector[p,li] = ToSymbol["e", p][li]


ScalarProduct = Pair

Attributes[Pair] = {Orderless}

	(* due to different operator priority in FORM: *)
Pair[a_, x_^n_] := Pair[a, x]^n

Pair[a_, x_Plus] := Pair[a, #]&/@ x

Pair[a_, n_?NumberQ x_] := n Pair[a, x]

	(* must say a_ (not _) here, or def will appear in
	   Abbreviations[] *)
Pair[a_, 0] = 0

Pair[e[i_], e[j_]] := Pair[e[i], e[j]] = ToSymbol["e", i, j];
Pair[e[i_], k[j_]] := Pair[e[i], k[j]] = ToSymbol["k", i, j]

Pair[s[i_], s[j_]] := Pair[s[i], s[j]] = ToSymbol["s", i, j];
Pair[s[i_], k[j_]] := Pair[s[i], k[j]] = ToSymbol["sk", i, j]

MomSquare[0] = 0

MomSquare[p_] := Pair[p, p]


Eps[___, a_, a_, ___] = 0

	(* must say a___ (not ___) here, or def will appear in
	   Abbreviations[] *)
Eps[a___, 0, ___] = 0

Eps[a___, x_Plus, b___] := Eps[a, #, b]&/@ x

Eps[a___, n_?NumberQ x_, b___] := n Eps[a, x, b]

Eps[a__] := Signature[{a}] Eps@@ Sort[{a}] /; !OrderedQ[{a}]

Eps[e[i1_], e[i2_], e[i3_], e[i4_]] := Eps[e[i1], e[i2], e[i3], e[i4]] =
  ToSymbol["e", i1, i2, i3, i4];
Eps[e[i1_], e[i2_], e[i3_], k[i4_]] := Eps[e[i1], e[i2], e[i3], k[i4]] =
  ToSymbol["e", i1, i2, i3, "k", i4];
Eps[e[i1_], e[i2_], k[i3_], k[i4_]] := Eps[e[i1], e[i2], k[i3], k[i4]] =
  ToSymbol["e", i1, i2, "k", i3, i4];
Eps[e[i1_], k[i2_], k[i3_], k[i4_]] := Eps[e[i1], k[i2], k[i3], k[i4]] =
  ToSymbol["e", i1, "k", i2, i3, i4];
Eps[k[i1_], k[i2_], k[i3_], k[i4_]] := Eps[k[i1], k[i2], k[i3], k[i4]] =
  ToSymbol["k", i1, i2, i3, i4]


ga[0] = 0

Unprotect[NonCommutativeMultiply];
___ ** 0 ** ___ = 0;
Protect[NonCommutativeMultiply]


(* Abbreviationing business *)

FME[-x_] := -FME[x]

FME[x_] :=
Block[ {c = Cases[x, ga[li_Symbol] :> li, Infinity], once},
  FME[ x /.
    MapIndexed[#1 -> LI@@ #2 &,
      once[l_] := (once[l] = Sequence[]; l); once/@ c] ] /; Length[c] =!= 0
]

FME[x_] := FME[x] = $O1ME^Count[x, ga[k[_]], Infinity] *
  (FME2[ x /. (#[[1, 1, 1]] -> #[[2]] &)/@ DownValues[Small] ] =
  Unique["F"])


SME[x_?AtomQ] = x

SME[-x_] := -SME[x]

SME[x_] := SME[x] = Unique["O"]


SMEPLUS[x_?AtomQ] = x

SMEPLUS[x_] := (SMEPLUS[x] =
Block[ {ft = FactorTerms[x]},
  If[ Head[ft] === Times, ft[[1]] SMEPLUS[ ft[[2]] ], ft ]
]) /; Union[Head/@ List@@ x] === {Times}

SMEPLUS[x_] := SMEPLUS[x] = Unique["P"]


dvhold[sym_] :=
  Cases[ DownValues[sym], (ab_ :> p_Symbol) :>
    (p -> HoldForm@@ ab Scale^Count[ab, k[_], {-2}]) /;
    FreeQ[ab, Pattern] ];

dv[sym_] :=
  Cases[ DownValues[sym], (_[_[ab_]] :> p_Symbol) :>
    (p -> ab Scale^Count[ab /. Spinor[__] -> 1, k[_], {-2}]) /;
    FreeQ[ab, Pattern] ];

Abbreviations[] :=
Block[ {ab = Join[dvhold[Pair], dvhold[Eps],
dv[FME2], dv[SME], dv[SMEPLUS]]},
  If[ $O1ME === 1, ab /. Scale -> 1,
    Prepend[ab, Scale -> $O1ME^-1]]
]


UseAbbreviations::warning =
"Warning: there are ambiguous definitions."

UseAbbreviations::scale =
"You must set $O1ME to `1` to use these abbreviations."

UseAbbreviations[file_String] := UseAbbreviations[Get[file]]

UseAbbreviations[abbr_] :=
Block[ {sc, ab},
  sc = 1/ReleaseHold[Scale /. abbr /. Scale -> 1];
  If[ sc =!= $O1ME,
    Message[UseAbbreviations::scale, sc];
    Return[$Failed] ];
  ab = DeleteCases[abbr, Scale -> _] /. Scale -> 1;
  ab /. HoldForm -> Identity;
  Apply[
    If[ AtomQ[#1],
      Which[
        !FreeQ[#2, NonCommutativeMultiply], 
          FME[#2] = sc^Count[#2, ga[k[_]], Infinity] (FME2[#2] = #1),
        !FreeQ[#2, Plus], SMEPLUS[#2] = #1,
        True, SME[#2] = #1 ],
    (* else *)
      If[#1 =!= #2, Message[UseAbbreviations::warning] ]
    ]&,
    Select[ab, FreeQ[#, HoldForm]&], 1 ]
]


IntSec[] = 0;
IntSec[a_, b_] := 0 /; Head[a] =!= Plus || Head[b] =!= Plus;
IntSec[a__] := Block[ {i = Intersection[a]}, If[Head[i] === Plus, i, 0] ]

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
      Set@@ {pl@@ #2, #1})&,
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
        iss = IntSec[repl[[i, 2]], -repl[[i + 1, 2]] ];
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
      is = ToSymbol["t", ++l];
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


(* Things to do when amplitude comes back from FORM *)

	(* for maximum simplification set e.g. o1 = o2 = Simplify *)
o1 = o2 = Identity

Global`ncm[x_] := x;
Global`ncm[x__] := NonCommutativeMultiply[x]

Global`fme[x_] := FME[x /. FromFormRules /. Global`e$ -> Eps]

Global`d$ = MetricTensor

Global`gi$[_] = Sequence[]

Global`pave4[i__Integer, p1_, p2_, p3_, m1_, m2_, m3_, m4_] := D0i[
  ToSymbol["dd", Sort[{i}]],
  MomSquare[p1], MomSquare[p1 - p2], MomSquare[p2 - p3],
  MomSquare[p3], MomSquare[p2], MomSquare[p1 - p3],
  m1, m2, m3, m4 ]

Global`pave3[i__Integer, p1_, p2_, m1_, m2_, m3_] := C0i[
  ToSymbol["cc", Sort[{i}]],
  MomSquare[p1], MomSquare[p1 - p2], MomSquare[p2],
  m1, m2, m3 ]

Global`pave2[i__Integer, p_, m1_, m2_] := B0i[
  ToSymbol["bb", Sort[{i}]], MomSquare[p], m1, m2 ]

Global`pave1[i__Integer, m_] := A0i[
  ToSymbol["aa", Sort[{i}]], m ]

Global`B0m[p_, m1_, m2_] := B0[MomSquare[p], m1, m2];
Global`B1m[p_, m1_, m2_] := B1[MomSquare[p], m1, m2];
Global`B00m[p_, m1_, m2_] := B00[MomSquare[p], m1, m2];
Global`B11m[p_, m1_, m2_] := B11[MomSquare[p], m1, m2]

B0[p_, m1_, m2_] := B0[p, m2, m1] /; !OrderedQ[{m1, m2}]

Derivative[1, 0, 0][B0] = DB0;
Derivative[1, 0, 0][B1] = DB1;
Derivative[1, 0, 0][B11] = DB11;
Derivative[1, 0, 0][B00] = DB00


a0[0] = A0[0] = 0


(* Reading and analyzing FeynArts files *)

FormMandelstam = FormDotSimplify = FormDiracSimplify = FormElimMom = ""

FromFormRules = FormVectors = {}


SetPair22[p1_, p2_, sign_, kin_] :=
  p1/: Pair[p1, p2] = sign kin - sign MomSquare[p1] - sign MomSquare[p2]

SetPair12[mom_, i_] :=
Block[ {lhs = Pair@@ Drop[mom, {i}], rhs},
  rhs = Plus@@
    Array[ If[# == i, 1, -1]/2 MomSquare[ mom[[#]] ] &, Length[mom]];
  If[MatchQ[lhs, -_], lhs = -lhs; rhs = -rhs];
  TagSet@@ {lhs[[1]], lhs, rhs}
]


FeynAmpList::noFAfile = "Hold it! This is no FeynArts amplitude."

FeynAmpList[h__][a___] :=
Block[ {proc, mom, ps, eps, dp, mdp, dot, x, x1, x2, s, t, u, amps,
PropagatorDenominator, MetricTensor = "d_"},
  proc = Process /. {h};
  If[ Head[proc] =!= Rule,
    Message[FeynAmpList::noFAfile];
    Return[$Failed] ];

  mom = #[[2]]&/@ Join[ proc[[1]], -proc[[2]] ];
  proc = Select[Join@@ proc, Head[ #[[2]] ] === Symbol &];
  Scan[Clear, ps = #[[2]]&/@ proc];
  Clear[U];
  Eps[k[1], k[2], k[3], k[4]] = .;

  dot[x1___, -p_, x2___] := -dot[x1, p, x2];
  eps = (dot[ x = ToSymbol["e", #], # ] = 0; x)&/@ ps;
  PowerCountingFor[ps];
  FormVectors = Flatten[{q1, eps, ps}];
  FromFormRules = Join[
    MapIndexed[#1 -> k@@ #2 &, ps], MapIndexed[#1 -> e@@ #2 &, eps] ];

  If[ $OnShell,
    Apply[UpSet, {MomSquare[ #[[2]] ], #[[3]]^2}&/@ proc, 1] ];

  amps = {a} /. Conjugate[ep:(Alternatives@@ eps)[__]] -> ep;

	(* introduce kinematical invariants: *)
  Switch[ Length[ps],
    2,
      amps = amps /. Rule@@ Reverse[ps],
    3,
      SetPair12[mom, 1];
      SetPair12[mom, 2];
      SetPair12[mom, 3],
    4 | 5,
      {s, t, u} =
        If[ Length[ps] === 4, 
          Eps[k[1], k[2], k[3], k[4]] = 0;
          If[ $OnShell, U/: S + T + U = Plus@@ (Last[#]^2 &)/@ proc ];
          {S, T, U},
        (* else *)
          {Sf, Tf, Uf} ];
      SetPair22[ ps[[1]], ps[[2]], 1/2, S ];
      SetPair22[ ps[[3]], ps[[1]], -1/2, T ];
      SetPair22[ ps[[4]], ps[[1]], -1/2, U ];
      SetPair22[ ps[[3]], ps[[4]], 1/2, s ];
      SetPair22[ ps[[4]], ps[[2]], -1/2, t ];
      SetPair22[ ps[[3]], ps[[2]], -1/2, u ]
  ];

  dp = If[ $Transversality,
    Apply["id " <> ToString[#1] <> " . " <> ToString[#2] <> " = 0;\n" &,
      Transpose[{eps, ps}], 1],
    {} ];
  mdp = Cases[ UpValues[#] /. Pair -> Dot, (_[p_] :> sp_) :>
    "id " <> ToString[p//InputForm] <> " = " <>
      ToString[sp//InputForm] <> ";\n" /; FreeQ[p, Pattern] ]&/@ ps;

  FormDiracSimplify = FormDotSimplify = FormElimMom = "";

  If[ Length[mom] > 2,
    x = Solve[Plus@@ mom == 0, #][[1, 1]]&/@ ps;
    FormElimMom = "id " <> ToString[ x[[-1, 1]] ] <> " = " <>
      ToString[ x[[-1, 2]]//InputForm ] <> ";\n";

    If[ FreeQ[amps, ga],
      FormDiracSimplify = StringJoin[
        ( "id ga(C?, " <>
          ToString[ #[[1]] ] <> ") = " <>
          ToString[ Map[If[Head[#] === Symbol, ga["C", #], #]&,
                      #[[2]], {-1}] ] <>
          ";\n#call DiracEquation{}\n"& )/@ x ] ];

    If[ Length[ps] > 4,
      dp = Insert[Reverse[dp], {mdp[[-1]], FormElimMom}, 2];
      mdp = Drop[mdp, -1]
    ];

    mom = Array[ {Plus@@ Delete[mom, #], -mom[[#]]}&, Length[mom] ];
    ( {x1, x2} = #;
      If[ $Transversality && AtomQ[x2],
        Scan[
          If[ (x = dot[#, x2]) =!= 0,
            FormDotSimplify = FormDotSimplify <> "id " <>
              ToString[x /. dot -> Dot] <> " = " <>
              ToString[Distribute[dot[#, x1]] /. dot -> Dot] <>
              ";\n"] &,
          eps ];
        FormDotSimplify = FormDotSimplify <> ".sort\n" ];
      x = Select[x1, AtomQ, 1];
      If[ x =!= 0, TagSet@@ Prepend[#, x]]
    )&/@ Join[mom, -mom];
  ];
  FormMandelstam = StringJoin[dp, mdp];

  PropagatorDenominator[pp_, mm_] :=
    PropagatorDenominator[pp, mm] = DEN[MomSquare[pp], mm^2];

  amps /. (1/xi[f_])^r_Rational :> xi[f]^-r
]


GList = Array[ToSymbol["Global`GM", #]&, 50]


(* DeleteIns[ins_List] := Delete[MapAt[-I # &, ins, p2], p] *)

DeleteIns[ins_List] := Delete[ins, p]

DeleteIns[ins_] := DeleteIns/@ ins


FeynAmp[g_GraphName, ___, amp_, coup_ -> ins_] :=
Block[ {t, u, p, p2, r, c = 0, coupvar = Take[GList, Length[coup]]},
  t = Transpose[{ins} /. Rule | Insertions[_] -> Sequence];
  u = Union/@ t;
  p = Position[u, {_}, {1}];
  Apply[(coupvar[[#]] = u[[#, 1]])&, p, {1}];
(*  p2 = Complement[Position[coup, G[_][_][__][__], {1}], p]; *)
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
(*
  p2 = Complement[p2, p];
  FeynAmp[ g[[-1]], g[[2]],
    amp /. Thread[coup -> MapAt[I # &, coupvar, p2]],
    Delete[coupvar, p] -> DeleteIns[ins] ]
*)
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

IList = Array[ToSymbol["Global`Ins", #]&, 500]

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

FormDecl[type_, _[v___, l_]] :=
Block[ {llen = Length[type]},
  WriteString[hh, type];
  ( WriteString[hh, #, ", "];
    If[ (llen += StringLength[ToString[#]] + 2) > 70,
      llen = 0; WriteString[hh, "\n  "] ] )&/@ {v};
  WriteString[hh, l, ";\n"];
]

FormPattern[sym_Symbol, val_] :=
  ( AppendTo[vars, sym];
    "id " <> ToString[sym] <> " = " <> ToString[val] <> ";\n" )

Arg2[_, x_] = x

FormPattern[(h_Symbol)[args__], val_] :=
Block[ {patt, lhs, c1 = 0, c2 = 0},
  lhs = h[args] /. Pattern -> Arg2 /.
    _Blank :> "Patt" <> ToString[++c1] <> "?" /.
    _BlankSequence | _BlankNullSequence :>
      StringJoin[Table["?", {++c2}]];
  vars = Join[vars, ToExpression["Patt" <> ToString[#]]&/@ Range[c1],
    Cases[lhs, _Symbol]];
  AppendTo[func, h];
  "id " <> ToString[lhs] <> " = " <> ToString[val] <> ";\n"
]


DeclareVars[expr_, indices_] :=
Block[ {theexpr, vars, indx, func = {}},
  theexpr = {expr, Last/@ Flatten[UpValues/@ FormVectors]};
  vars = Cases[theexpr, x_Symbol /; Context[x] =!= "System`",
    Infinity, Heads -> True];
  smalls = StringJoin[ Cases[DownValues[Small],
    (_[_[lhs_]] :> rhs_) :> FormPattern[lhs, rhs]] ];
  vars = Append[Complement[vars, FormVectors], Pi];
  func = Union[func, Select[vars, !FreeQ[theexpr, #[__]]&]];
  If[!FreeQ[theexpr, Conjugate], AppendTo[func, Conjugate]];
  indx = Intersection[vars, indices];
  vars = Complement[vars, func, indx];
  FormDecl["i ", indx];
  FormDecl["v ", FormVectors];
  FormDecl["s ", vars];
  FormDecl["cf ", DeleteCases[func,
    FeynAmp | Spinor | DiracMatrix | DiracSlash |
    ChiralityProjector | DiracTrace | tr]];
]


TakeIns[FeynAmp[__, _ -> ins_]] := List@@ ins

TakeIns[_] = {}


Options[OneLoop] = {
  Comment -> "",
  AmplitudeLevel -> Automatic,
  DiracSimplify -> True,
  DotSimplify -> Automatic,
  CancelEps -> False,
  ChiralME -> True,
  NoExpand -> {},
  EditCode -> False,
  RetainFile -> False }

kinobjs = Spinor | ChiralityProjector | DiracMatrix | DiracSlash |
  DiracTrace | SpinorChain | q1 | p1 | p2 | p3 | k1 | k2 | k3 | k4 |
  _String

OneLoop[___Rule] = OneLoop[_[], ___Rule] = 0

OneLoop[amps_, opt___Rule] :=
Block[ {com, lev, diracsimp, canceps, chime, dotsimp, noexp, edcode,
retain, hh, amplist, res, res3, res4, smalls,
Hide, Global`Unhide, hidec = 0, extferm = 0,
inssum = 0, ins = {}, allins = {}, iabbr = {}, names = {}, alltrace = {}},

  {com, lev, diracsimp, dotsimp, canceps, chime, noexp, edcode, retain} =
    {Comment, AmplitudeLevel, DiracSimplify, DotSimplify, CancelEps,
      ChiralME, NoExpand, EditCode, RetainFile} /.
      {opt} /. Options[OneLoop];
  If[ dotsimp === Automatic, dotsimp = !FermionicQ[amps] ];

  Print["preparing FORM code in ", $TempFile];
  hh = OpenFormTemp;

  res = Flatten[{amps}] /. {
    2^(1/2) -> r2, 2^(-1/2) -> 1/r2,
    Complex[a_, b_] -> a + "i_" b,
    Eps -> "e_",
    _SumOver -> 1 };

  noexp = Alternatives@@ Flatten[{noexp}];
  If[ Length[noexp] =!= 0,
    Hide[p__] := Plus[p] /; FreeQ[{p}, noexp] || !FreeQ[{p}, kinobjs];
    Hide[a:-_, b__] := -Hide@@ (-#&)/@ {a, b};
    Hide[p__] := Hide[p] =
      ( Global`Unhide[++hidec] = Plus[p] /. "Unhide" -> Global`Unhide;
        "Unhide"[hidec] );
    res = res /. Plus -> Hide;
  ];

  If[ !FreeQ[res, Insertions],
    If[lev === Automatic,
      lev = Union[ Cases[res,
        Insertions[t_][__] -> t, Infinity] ][[-1]] ];
    res = PickLevel[lev][res];
    allins = Union[Flatten[TakeIns/@ res]];
    ins = Select[allins, Head[#] =!= Symbol &];
    iabbr = Take[IList, Length[ins]];
    ins = Thread[ins -> iabbr];
  ];
  res3 = Select[res, Length[#] === 3 &];
  res4 = Select[res, Length[#] === 4 &];

  WriteString[hh,
    "#-\n* ", com,
    "\n#define VADecomp \"", If[chime === False, "1", "0"],
    "\"\n#define DIM \"", $Dimension /. D -> 0, "\"\n",
    If[$Dimension === D, "s D;\nd D", "d " <> ToString[$Dimension] ],
    ";\n"];

  DeclareVars[{#[[3]]&/@ res, allins}, ToExpression/@ Names["li*"]];
  WriteString[hh,
    "cf Unhide;\n",
    "f Spinor;\n",
    "#if 'VERSION_' > 1\nNT ga, omp, omm, ga5;\n",
    "#else\nf ga, omp, omm, ga5;\n#endif\n\n",
    ".global\n\n#procedure Mandelstam()\n",
    FormMandelstam,
    "#endprocedure\n\n#procedure DotSimplify()\n",
    If[TrueQ[dotsimp], FormDotSimplify, ""],
    "#endprocedure\n\n#procedure DiracSimplify()\n",
    If[TrueQ[diracsimp], FormDiracSimplify, ""],
    "#endprocedure\n\n"];

  amplist = {};
  If[ Length[res3] =!= 0,
    Scan[FormWrite, res3];
    WriteString[hh, "g amp3 = "];
    Write[hh, Plus@@ amplist, ";"];
    WriteString[hh, ".sort\n\n"];
    inssum = Global`amp3 ];
  Scan[FormWrite, res4];

  alltrace = Union[Flatten[alltrace]];
  If[ Length[alltrace] =!= 0,
    Scan[WriteString[hh, "trace4,", #, ";\n"]&, alltrace];
    If[ canceps,
      WriteString[hh, "if( count(e_, 1) > 0 );\n  discard;\nendif;\n"] ]
  ];

  WriteString[hh,
    "\n#define ExtFermions \"", extferm,
    "\"\n#if 'VERSION_' > 1\n",
    "#include ", $FormCalcDir, "OneLoop.h2\n#else\n",
    "#include ", $FormCalcDir, "OneLoop.h1\n",
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
  WriteString[hh,
    "\n", smalls,
    "\n#call TrivialSubst{}\n",
    "\nb fme, DEN, i_, Unhide, A0, B0m, B1m, B00m, B11m,\n",
    "  pave1, pave2, pave3, pave4, pave5;\n\n",
    "print +s amp;\n\n.end\n"];
  Close[hh];

  RunForm /. Reverse/@ ins
]


OpenFormTemp := OpenWrite[
  "!sed -e '/^[^#]/s/\"//g' -e 's/\\[/(/g' -e 's/]/)/g' " <>
  "-e 's/\\*\\*/\\*/g' -e 's/\\*\\\\//g' > " <> $TempFile,
  FormatType -> InputForm, PageWidth -> 73]; 

RunForm := (
  If[edcode, Pause[1]; Run[StringForm[$Editor, $TempFile]]; Pause[3]];
  WriteString["stdout", "\nrunning FORM... "];
  res = Block[ {Dot = Pair, r2 = Sqrt[2]},
    ReadForm[ToString[
      StringForm["!`1` -s `2`form.set `3`",
        $FormCmd, $FormCalcDir, $TempFile] ]] 
  ] /. FromFormRules /. Global`e$ -> Eps /.
    Global`sme -> SME /. Global`smeplus -> SMEPLUS;
  If[ res === $Failed,
    Print["warning: FORM reports an error!!!"];
    Abort[] ];
  Print["ok"];
  If[!retain, DeleteFile[$TempFile]];
  res
)


(* This is for backward compatibility with FeynCalc *)

C0[p1_, p2_, p1p2_, m1_, m2_, m3_] :=
  C0i[cc0, p1, p2, p1p2, m1, m2, m3]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_}] :=
  C0i[ToSymbol["cc", Sort[{i}]], p, m1, m2, m3]

D0[p1_, p2_, p3_, p4_, p1p2_, p2p3_, m1_, m2_, m3_, m4_] :=
  D0i[dd0, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4]

PaVe[i__Integer, {p__}, {m1_, m2_, m3_, m4_}] :=
  D0i[ToSymbol["dd", Sort[{i}]], p, m1, m2, m3, m4]

StandardMatrixElement[a_] StandardMatrixElement[b_] ^:=
  StandardMatrixElement[a b];
StandardMatrixElement[a_] := StandardMatrixElement[a /.
  StandardMatrixElement -> Identity] /; !FreeQ[a, StandardMatrixElement]

fcpair[e[i_], e[j_]] := fcsme[FCPair[e[i], e[j]]];
fcpair[e[i_], k[j_]] := $O1ME fcsme[FCPair[e[i], k[j]]];
fcpair[a__] := Pair[a]

fceps[a__] := Signature[{a}] (fceps@@ Sort[{a}]) /; !OrderedQ[{a}];
fceps[e[i1_], e[i2_], e[i3_], e[i4_]] :=
  I fcsme[FCEps[e[i1], e[i2], e[i3], e[i4]]];
fceps[e[i1_], e[i2_], e[i3_], k[i4_]] :=
  I $O1ME fcsme[FCEps[e[i1], e[i2], e[i3], k[i4]]];
fceps[e[i1_], e[i2_], k[i3_], k[i4_]] :=
  I $O1ME^2 fcsme[FCEps[e[i1], e[i2], k[i3], k[i4]]];
fceps[e[i1_], k[i2_], k[i3_], k[i4_]] :=
  I $O1ME^3 fcsme[FCEps[e[i1], k[i2], k[i3], k[i4]]]

fcsme/: fcsme[a_] fcsme[b_] := fcsme[a b]

FeynCalcGet[mask___] :=
Block[ {Global`OneLoopResult, GraphName, Global`SW, Global`CW,
Pair = fcpair, Eps = fceps, StandardMatrixElement = Identity},
  Global`SW/: Global`SW^-2 = Global`S2;
  Global`SW/: Global`SW^-4 = Global`S4;
  Global`CW/: Global`CW^-2 = Global`C2;
  Global`CW/: Global`CW^-4 = Global`C4;
  GraphName[__] = 0;
  Plus@@ ((Print[#]; Get[#]; Global`OneLoopResult[0])&)/@
    FileNames[mask] /.
    (p:S | T | U - m_)^(n_?Negative) :> DEN[p, m]^(-n) /.
    (p:S | T | U)^(n_?Negative) :> DEN[p, 0]^(-n) /.
    (m_ - p:S | T | U)^(n_?Negative) :> -DEN[p, m]^(-n)
] /. FCPair -> Pair /. FCEps -> Eps /. fcsme -> SME


FeynCalcPut[expr_, file_] :=
Block[ {PaVe, C0i, D0i, C0, D0},
  C0i[Global`cc0, args__] := C0[args];
  D0i[Global`dd0, args__] := D0[args];
  C0i[id_, p1_, p2_, p1p2_, m1_, m2_, m3_] := PaVe[
    Sequence@@ (ToExpression/@ Characters[StringDrop[ToString[id], 2]]),
    {p1, p2, p1p2}, {m1, m2, m3}];
  D0i[id_, p1_, p2_, p3_, p4_, p1p2_, p2p3_, m1_, m2_, m3_, m4_] := PaVe[
    Sequence@@ (ToExpression/@ Characters[StringDrop[ToString[id], 2]]),
    {p1, p2, p3, p4, p1p2, p2p3}, {m1, m2, m3, m4}];
  Put[expr, file]
]


(* Classifying amplitudes *)

FermionFamily[-x_] = x;
FermionFamily[_] = ""


DoTrivialSums[{ins___, relcf_}] :=
Block[ {newcf},
  newcf = relcf /. SumOver[i_, v_] :> v /; FreeQ[{amp, ins}, i];
	(* 5 newcf is a trick to make the head Times, i.e. to include
	   the case where newcf = SumOver[...] *)
  {ins, newcf, Cases[5 newcf, _SumOver]}
]

DoTrivialSums[a_, ins___] :=
Block[ {so = Cases[a, _SumOver]},
  Fold[
    If[FreeQ[ {#1, ins}, #2[[1]] ], #1 #2[[2]], #1 #2]&,
    DeleteCases[a, Alternatives@@ so],
    so ]
]

SeparateSums[ FeynAmp[n_, top_, a_, coup_ -> ins_] ] :=
Block[ {theins, so = {}, amp},
  amp = DoTrivialSums[a, ins];
  theins = DoTrivialSums/@ ins;
  While[ Length[theins] =!= 0,
    p = First/@ Position[ theins, theins[[1, -1]] ];
    AppendTo[ so,
      FeynAmp[ n, top, amp, coup -> (Drop[#, -1]&)/@ theins[[p]] ] ];
    theins = Delete[theins, List/@ p] ];
  Sequence@@ so
] /; !FreeQ[{a, ins}, SumOver]

SeparateSums[ FeynAmp[n_, top_, a_] ] :=
  FeynAmp[n, top, DoTrivialSums[a]] /; !FreeQ[a, SumOver]

SeparateSums[ a_ ] = a


ClassIndex[ FeynAmp[_, top_, amp_, rulz___] ] :=
Block[ {cname = "", so},
  If[howcls =!= IndexSumsOnly,
    If[FermionicQ[amp], cname = "F"];
    If[howcls === Tough && !FreeQ[amp, d0],
      If[cname === "F", cname = cname <> Union[Cases[
        Cases[amp, (DiracTrace | SpinorChain | Dot)[t__] :> t, Infinity],
        _. _DiracSlash + m_. :> FermionFamily[m] ]] ];
      cname = cname <> ToString[top]
    ]
  ];
  so = Union[Cases[{amp, rulz}, _SumOver, Infinity]];
  cname = cname <> ("_" <> ToString[ #[[1]] ] &)/@ so;
  sums[cname] = Times@@ so;
  cname
]


Options[ProcessFile] = {Classification -> Standard}

ProcessFile[filename_String, outtag___String, opt___Rule] :=
  ProcessFile[Get[filename], outtag, opt]

ProcessFile[amps_List, outtag___String, opt___Rule] :=
Block[ {class, sums, na, tough},
  howcls = Classification /. {opt} /. Options[ProcessFile];
  theamps = If[ FreeQ[amps, Insertions], amps,
    lev = AmplitudeLevel /. {opt} /. Options[OneLoop];
    If[lev === Automatic,
      lev = Union[ Cases[amps,
        Insertions[t_][__] -> t, Infinity] ][[-1]] ];
    PickLevel[lev][amps] ];
  class[_] = {};
  Scan[AppendTo[class[ClassIndex[#]], #]&, SeparateSums/@ theamps];
  Cases[ DownValues[class], (_[_[cname_String]] :> amp_) :>
    ( Print[""];
      Share[];
      If[ Length[{outtag}] === 0,
        Print["processing class ", na = "[" <> cname <> "]"];
        sums[cname] OneLoop[amp, Comment -> na, opt],
      (* else *)
        Print["processing class ", na = outtag <> cname <> ".m"];
        Put[sums[cname] OneLoop[amp, Comment -> na, opt], na];
        na
      ] ) ]
]


(* Helicity matrix elements *)

rho[fline_, k[n_], 0, sign_] :=
  "g_"[fline, k[n]] ** ("gi_"[fline] + sign Hel[n] "g5_"[fline])/2

rhoc[fline_, k[n_], 0, sign_] :=
  -("gi_"[fline] + sign Hel[n] "g5_"[fline]) ** "g_"[fline, k[n]]/2

rho[fline_, k[n_], m_, sign_] :=
  ("g_"[fline, k[n]] + sign m "gi_"[fline]) **
    ("gi_"[fline] + Hel[n] "g5_"[fline] ** "g_"[fline, e[n]])/2

rhoc[fline_, k[n_], m_, sign_] :=
  ("gi_"[fline] - Hel[n] "g_"[fline, e[n]] ** "g5_"[fline]) **
    (-"g_"[fline, k[n]] + sign m "gi_"[fline])/2

ToTrace[plain_, conj_] :=
Block[ {me, gi, fline = 0},
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
               gi :> "gi_"[fline],
               ga5 :> "g5_"[fline],
               omp :> "g6_"[fline]/2,
               omm :> "g7_"[fline]/2,
               (r:rho | rhoc)[x__] :> r[fline, x] } );
  tr[fline, me]
]


Options[HelicityME] = {
  AbbreviationsToUse :> Abbreviations[],
  EditCode -> False,
  RetainFile -> False }

(* this is for squared matrix elements: *)

HelicityME[plain_, conj_, opt:(_Rule | _RuleDelayed)...] :=
Block[ {abbr, fermabbr, theME, conjME, hh, edcode, retain, x, mat, res,
part, e, heli, hels, indx, smalls, o1},

  {abbr, edcode, retain} =
    {AbbreviationsToUse, EditCode, RetainFile} /.
    {opt} /. Options[HelicityME];

  fermabbr = Select[abbr, !FreeQ[#, Spinor]&];
  abbr = Select[abbr, FreeQ[#, Spinor]&];
  theME = If[ plain === All, fermabbr,
    x = plain //. abbr;
      Union[Select[fermabbr, !FreeQ[ x, #[[1]] ]&]] ];
  conjME = If[ conj === All, fermabbr,
    x = conj //. abbr;
      Union[Select[fermabbr, !FreeQ[ x, #[[1]] ]&]] ] /.
    {omp -> omm, omm -> omp, ga5 -> -ga5, ep_Eps -> -ep} /.
    n_NonCommutativeMultiply :> Reverse[n] /.
    LI[n_] :> LI[n + 10];
  mat = Outer[Mat, First/@ theME, First/@ conjME]//Flatten;
  theME = Outer[ToTrace, Last/@ theME, Last/@ conjME]//Flatten;

  Print["preparing output for FORM in ", $TempFile];
  hh = OpenFormTemp;
  WriteString[hh, "#-\n"];

  part = Cases[fermabbr, Spinor[k[n_], _, _] -> n, Infinity]//Union;
  heli = (# -> ToSymbol[ "Hel", #[[1]] ])&/@ Cases[Hel/@ part, _Hel];
  hels = Last/@ heli;
  indx = (# -> ToSymbol[ "LI", #[[1]] ])&/@
    Union[Cases[theME, _LI, Infinity]];
  theME = theME /. indx /. heli /.
    Reverse/@ FromFormRules /. Eps -> "e_";
  DeclareVars[{theME, hels}, Last/@ indx];

  WriteString[hh,
    "\n.global\n\n#procedure Simplify()\ncontract,0;\n",
    FormMandelstam,
    smalls,
    FormDotSimplify];
  FormDecl["b ", hels];
  WriteString[hh, "print +s;\n.store\n#endprocedure\n"];    

  Apply[
    ( WriteString["stdout", x = "Mat" <> ToString/@ List@@ #1 <> " "];
      WriteString[hh, "\nl ", x, "= "];
      Write[hh, #2[[2]], ";"];
      Array[ WriteString[hh, "trace4,", #, ";\n"]&, #2[[1]] ];
      WriteString[hh, "#call Simplify{}\n"] )&,
    Transpose[{mat, theME}], 1 ];

  WriteString[hh, "\n.end\n"];
  Close[hh];

  (e[#] = s[#])&/@ part;
  o1[x_] := o1[x] = Simplify[x];

  theME = RunForm;
  If[ res =!= $Failed,
    res = Thread[mat -> res] /. Reverse/@ heli;
    x = Select[Union[Position[res, #[[2]], {1}]&/@ res], Length[#] > 1 &];
    If[ Length[x] =!= 0,
      res = Fold[
        ReplacePart[#1, res[[#2[[1, 1]], 1]], Drop[#2, 1]]&,
        res, x] ]
  ];
  If[ !FreeQ[abbr, Scale], res = res /. abbr[[1]] ];
  res
]


(* this is for simple matrix elements: *)

(* HelicityME[plain_, opt___] := notyet *)


SquaredME[plain_, conj_] :=
Block[ {fc, c0, c1, coeff},
  fc = First/@ Select[Abbreviations[], !FreeQ[#, Spinor]&];
  fc = Select[fc, !FreeQ[{plain, conj}, #]&];
  If[ Length[fc] === 0, Return[plain Conjugate[conj]] ];
  {c0, c1} = If[Head[#] === Plus, List@@ #, {#}]&/@
    (Collect[{plain, conj}, fc, coeff] /. f_ coeff[c_] -> coeff[f, c]);
  Plus@@ Flatten[Outer[
    Mat[ #1[[1]], #2[[1]] ] #1[[2]] Conjugate[ #2[[2]] ] &,
    c0, c1] ]
]

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

$O1ME = 1
  (* factor which scales SMEs, see usage *)

$TempFile = "t" <> ToString[$ProcessID] <> ".frm"

EndPackage[]


(* global definitions for specific models *)

(* definitions specially for the Standard Model *)

EL/: EL^4 = 16 Pi^2 a2

MW^(n_?EvenQ) ^= MW2^(n/2);
MP^(n_?EvenQ) ^= MP2^(n/2);
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
  { EL, a2, SW, CW, S2, C2, S4, C4,
    MW, MW2, MP, MP2, MZ, MZ2, MH, MH2,
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

(RealQ[#] = True)&/@
  { TB, CB, SB, CA, SA, C2A, S2A, CAB, SAB, CBA, SBA,
    Mh0, MH0, MA0, MG0, MHp, MGp,
    _MSf, _MSf2, _MSNE, _MSLE1, _MSLE2,
    _MSQU1, _MSQU2, _MSQD1, _MSQD2,
    _MCha, _MCha2, _MNeu, _MNeu2 }

Null

