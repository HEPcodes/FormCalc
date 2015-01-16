* CalcFeynAmp.frm
* the FORM part of the CalcFeynAmp function
* this file is part of FormCalc
* last modified 29 Apr 08 th


***********************************************************************

#procedure DotSimplify
#do i = 1, `Legs'
#ifdef `k`i''
b `Patterns', Times, NoExpand, Den, GA, Spinor,
  SumOver, `SUNObjs',
  pave, cuti, `LoopInt', GA, Spinor, i_;
.sort

#call eiki
id k`i' = `k`i'';
#endif
#enddo

b `Patterns', Times, NoExpand, Den, GA, Spinor,
  SumOver, `SUNObjs',
  pave, cuti, `LoopInt', GA, Spinor, i_;
.sort

#call eiki
#call kikj
#endprocedure

***********************************************************************

#procedure DiracSimplify

#call DiracOrder

#if `OnShell' == 1

#do i = 1, `Legs'
#ifdef `k`i''
b `Tensors';
.sort

keep brackets;

id k`i' = `k`i'';
#call DiracOrder
#endif
#enddo

#endif

#endprocedure

***********************************************************************

#procedure DiracOrder
label 1;

#if `OnShell' == 1

* Apply Dirac equation to right spinor

repeat;
  id GA([om]?, ?a, [p1]?, ?b) * Spinor([p1]?, [m1]?, [s1]?) =
    ( 2*GD([om], ?a, [p1]) * distrib_(-1, 1, GD, GD, ?b) +
      sign_(nargs_(?b)) * [s1] * [m1] * GA([om], ?a, ?b) ) *
    Spinor([p1], [m1], [s1]);
  id GD(?a, [p1]?) * GD([mu]?) * GD(?b) =
    d_([p1], [mu]) * GA(?a, ?b) * TAG;
endrepeat;

* Apply Dirac equation to left spinor

repeat;
  id Spinor([p1]?, [m1]?, [s1]?) * GA([om]?{6,7}[[x]], ?a, [p1]?, ?b) =
    Spinor([p1], [m1], [s1]) * sign_(nargs_(?a)) *
    ( [s1] * [m1] * GA({7,6}[[x]], ?a, ?b) -
      2*distrib_(-1, 1, GD, GD, ?a) * GD([p1], [om], ?b) );
  id GD([mu]?) * GD(?a) * GD([p1]?, [om]?, ?b) =
    d_([p1], [mu]) * GA([om], ?a, ?b) * TAG;
endrepeat;
#endif

* Eliminate contractions within each Dirac chain using the
* formulas from M. Veltman's Gammatrica [Nucl Phys B319 (1989) 253]

id GA([om]?, [mu]?, ?a) = GA([om]) * GB([mu], ?a);
while( count(GB, 1) );
  repeat;
    id GB([LA]?, [LA]?, ?a) = d_([LA], [LA]) * GB(?a);
    also GB([LA]?, [mu]?, [LA]?, ?a) = (2 - d_([LA], [LA])) * GB([mu], ?a);
    also GB([LA]?, [mu]?, [nu]?, [LA]?, ?a) = 4*d_([mu], [nu]) * GB(?a) * TAG
#if `Dim' == D
        + Dminus4 * GB([mu], [nu], ?a)
#endif
       ;
    also GB([LA]?, [mu]?, [nu]?, [rho]?, ?b, [LA]?, ?a) =
      -sign_(nargs_(?b)) * (
        2*GB([rho], [nu], [mu], ?b, ?a) +
        2*GD([mu], [nu], [rho]) * distrib_(-1, 1, GD, GD, ?b) * GD(?a)
#if `Dim' == D
        + Dminus4 * GB([mu], [nu], [rho], ?b, ?a)
#endif
      );
    id GD(?a) * GD([mu]?) * GD(?b) * GD(?c) = GB([mu], ?a, ?b, ?c);
  endrepeat;
  id GB([mu]?, ?a) = GC([mu]) * GB(?a);
  id GB() = 1;
endwhile;

* Order the gamma matrices canonically

repeat;
  id GC([p1]?) * GC([p1]?) = [p1].[p1];
  disorder GC([mu]?) * GC([nu]?) = 2*d_([mu], [nu]) * TAG - GC([nu]) * GC([mu]);
endrepeat;
repeat id GA(?a) * GC(?b) = GA(?a, ?b);

id ifmatch->1 TAG = 1;
#endprocedure

***********************************************************************

#procedure FierzPre
id Spinor(?a) * GA(?g) * Spinor(?b) =
  CH(Spinor(?a) * GA(?g) * Spinor(?b));

* these relations are obtained by Fierzing twice

id CH(Spinor(?a) * GA(6, [mu]?, [nu]?) * Spinor(?b)) *
   CH(Spinor(?c) * GA(7, [mu]?, [nu]?) * Spinor(?d)) =
  4*CH(Spinor(?a) * GA(6) * Spinor(?b)) *
    CH(Spinor(?c) * GA(7) * Spinor(?d));

id CH(Spinor(?a) * GA([om]?, [mu]?, [nu]?) * Spinor(?b)) *
   CH(Spinor(?c) * GA([om]?, [rho]?, [mu]?, [nu]?) * Spinor(?d)) =
  4*CH(Spinor(?a) * GA([om]) * Spinor(?b)) *
    CH(Spinor(?c) * GA([om], [rho]) * Spinor(?d));
also CH(Spinor(?a) * GA([omA]?, [mu]?, [nu]?) * Spinor(?b)) *
   CH(Spinor(?c) * GA([omB]?, [rho]?, [mu]?, [nu]?) * Spinor(?d)) =
  4*CH(Spinor(?a) * GA([omA], [rho], [mu]) * Spinor(?b)) *
    CH(Spinor(?c) * GA([omB], [mu]) * Spinor(?d));

id CH(Spinor(?a) * GA(6, [rho]?, [mu]?, [nu]?) * Spinor(?b)) *
   CH(Spinor(?c) * GA(7, [sig]?, [mu]?, [nu]?) * Spinor(?d)) =
  4*CH(Spinor(?a) * GA(6, [rho]) * Spinor(?b)) *
    CH(Spinor(?c) * GA(7, [sig]) * Spinor(?d));

id CH(Spinor(?a) * GA([om]?, [mu]?, [nu]?, [rho]?) * Spinor(?b)) *
   CH(Spinor(?c) * GA([om]?, [mu]?, [nu]?, [rho]?) * Spinor(?d)) =
  16*CH(Spinor(?a) * GA([om], [mu]) * Spinor(?b)) *
     CH(Spinor(?c) * GA([om], [mu]) * Spinor(?d));
also CH(Spinor(?a) * GA([omA]?, [mu]?, [nu]?, [rho]?) * Spinor(?b)) *
   CH(Spinor(?c) * GA([omB]?, [mu]?, [nu]?, [rho]?) * Spinor(?d)) =
  4*CH(Spinor(?a) * GA([omA], [mu]) * Spinor(?b)) *
    CH(Spinor(?c) * GA([omB], [mu]) * Spinor(?d));

* The following general Fierz identity is from hep-ph/0412245.
repeat;
  once CH(Spinor(?a) * GA([omA]?, ?A) * Spinor(?d)) *
       CH(Spinor(?c) * GA([omB]?, ?B) * Spinor(?b)) =
    sum_([c], 1, 5, sum_([d], 1, 5,
      DUAL([d], MUD, NUD) * CHI([omA])*g_(1, ?A) *
      DUAL([c], MUC, NUC) * CHI([omB])*g_(1, ?B) *
      Spinor(?a) * BASIS([d], MUD, NUD) * Spinor(?b) *
      Spinor(?c) * BASIS([c], MUC, NUC) * Spinor(?d) ));
  trace4, 1;
  sum MUC, NUC, MUD, NUD;
endrepeat;
#endprocedure

***********************************************************************

#procedure FierzPost
id CH([x]?) = [x];

repeat;
  once GA([om]?, [mu]?, [nu]?, [rho]?, ?a) =
    sum_([c], 1, 5,
      DUAL([c], MUC, NUC) * CHI([om])*g_(1, [mu], [nu], [rho], ?a) *
      BASIS([c], MUC, NUC) );
  trace4, 1;
  sum MUC, NUC;
endrepeat;

contract 0;

* Chisholm's identity backwards to get rid of all Eps
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [rho]?, [LA]?) =
    sign_([om]) * sign_(nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [rho], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [rho], ?b) +
      d_([mu], [rho]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [rho]) * GA([om], ?a, [mu], ?b) );
endrepeat;

#call DiracSimplify
id D = 4;
#endprocedure

***********************************************************************

#procedure MomSquare
#ifdef `MomSum'
* Apply momentum conservation to generate as few terms as possible

id MOM([p1]?) = MOM([p1], nterms_([p1]),
  [p1] + (`MomSum'), nterms_([p1] + (`MomSum')),
  [p1] - (`MomSum'), nterms_([p1] - (`MomSum')));

symm MOM (2,1) (4,3) (6,5);
#endif

id MOM([p1]?, ?a) = [p1].[p1];
#call kikj
#endprocedure

***********************************************************************

#procedure Fewest(foo)
argument `foo';
#call Neglect
endargument;
id `foo'([x]?, [y]?) = `foo'([x], nterms_([x]), [y], nterms_([y]));
symm `foo' (2,1), (4,3);
id `foo'([x]?, 1, ?a) = [x];
id `foo'([x]?, ?a) = `foo'([x]);
#endprocedure

***********************************************************************

#procedure Factor(foo)
factarg `foo';
chainout `foo';
id `foo'([x]?number_) = [x];
id `foo'([x]?symbol_) = [x];
#endprocedure

***********************************************************************

#procedure DoInsertions
.store

#call Insertions
#call FillIns
#call Neglect

.sort

argument;
#call FillIns
endargument;
#endprocedure


***********************************************************************
*** the main program starts here

#if "`InsertAt'" == "Begin"
#call DoInsertions
#else
#call Neglect
#endif

#call Const
.sort

*----------------------------------------------------------------------

#define LoopInt "A0i, B0i, C0i, D0i, E0i, F0i"

#define CutInt "Acut, Bcut, Ccut, Dcut, Ecut, Fcut"

#define SUNObjs "SUNSum, SUNT, SUNTSum, SUNF, SUNEps"

#define Tensors "Spinor, GA, Sigma, e_, eT1,...,eT`Legs', eTc1,...,eTc`Legs'"

* variables appearing in the CalcFeynAmp input and output
s I, D, Renumber;
cf Mat, Den, Num, Pair, Eps, DiracChain, WeylChain;
cf SumOver, IndexDelta, IndexEps, `SUNObjs';
cf NoExpand, Times, Simplify;
cf `LoopInt', `CutInt';
f Spinor, DottedSpinor;
i Col1,...,Col`Legs', Ind1,...,Ind10;
nt Sigma;

* variables that make it into Mma but don't appear in the output
cf abb, pave, cuti, fme, sun;

* patterns
s [x], [y], [z], [w], [s1], [s2];
s [k1], [k2], [k1k2], [m1], [m2], [m3];
i [mu], [nu], [rho], [sig], [om], [omA], [omB], [LA];
i [i], [j], [k], [l], [I], [J], [K];
i [a], [b], [c], [d];
v [p1], [p2], [p3], [p4], [p5], [p6];
cf [f];
t [t];

* variables internal to FORM
s TAG, CUT;
i DUMMY, MUC, NUC, MUD, NUD;
cf TMP, MOM, ABB, ORD, ORD1, ORD2, CH;
t NUM, EQ, NEQ, EPS(antisymm);
nt GA, GB, GC, GD;
f WC, GR;
auto s FC;
set LOOPINT: `LoopInt';
set CUTINT: `CutInt';
set MOMS: k1,...,k`Legs';
set COLS: Col1,...,Col`Legs';

ntable BASIS(1:5, [mu]?, [nu]?);
ntable DUAL(1:5, [mu]?, [nu]?);
ntable CHI(6:7);

fill BASIS(1) = GA(6);
fill BASIS(2) = GA(7);
fill BASIS(3) = GA(6, [mu]);
fill BASIS(4) = GA(7, [mu]);
fill BASIS(5) = i_/2*(GA(6, [mu], [nu]) + GA(7, [mu], [nu]) -
                      (GA(6) + GA(7))*d_([mu], [nu]));

fill DUAL(1) = g_(1, 6_)/4;
fill DUAL(2) = g_(1, 7_)/4;
fill DUAL(3) = g_(1, 7_, [mu])/4;
fill DUAL(4) = g_(1, 6_, [mu])/4;
fill DUAL(5) = i_/4*(g_(1, [mu], [nu]) - d_([mu], [nu]));

fill CHI(6) = g6_(1)/2;
fill CHI(7) = g7_(1)/2;

.global

*----------------------------------------------------------------------

collect Times;

moduleoption polyfun=Times;
.sort

argument;
#call Square
endargument;

id Times(0) = 0;

normalize Times;

*----------------------------------------------------------------------

b q1, `LoopInt';
.sort

keep brackets;

* symmetrize the loop integrals

id A0i([p1]?, ?a) = A0i(?a);

#redefine loopn "0"
#do loopf = {`LoopInt'}
#redefine loopn "{`loopn'+1}"
#if `loopn' > 1
symm `loopf' <({`loopn'+1},1)>,...,<({`loopn'+`loopn'},`loopn')>;
id `loopf'(<[p1]?>,...,<[p`loopn']?>, ?a) = replace_(q1, 2*q1 - [p1]) *
  `loopf'(<[p2] - [p1]>,...,<[p`loopn'] - [p1]>, ?a);
#endif
#enddo

* cancel q^2's in the numerator

while( match(q1.q1) );
  once q1.q1 * A0i([m1]?) =
    1 + [m1] * A0i([m1]);
  once q1.q1 * B0i([p1]?, [m1]?, ?a) =
    replace_(q1, q1 - [p1]) * A0i(?a) +
    [m1] * B0i([p1], [m1], ?a);
#redefine loopn "0"
#do loopf = {`LoopInt'}
#if `loopn' > 1
  once q1.q1 * `loopf'(<[p1]?>,...,<[p`loopn']?>, [m1]?, ?a) =
    replace_(q1, q1 - [p1]) * `lastf'(<[p2]-[p1]>,...,<[p`loopn']-[p1]>, ?a) +
    [m1] * `loopf'(<[p1]>,...,<[p`loopn']>, [m1], ?a);
#endif
#redefine loopn "{`loopn'+1}"
#redefine lastf "`loopf'"
#enddo
endwhile;

id A0i(0) = 0;

*----------------------------------------------------------------------

#call DotSimplify

#if `CutTools' == 1

id q1.[p1]? = cuti(q1.[p1], Pair(q1, [p1]));
id e_(q1, [p1]?, [p2]?, [p3]?) =
  cuti(e_(q1, [p1], [p2], [p3]), Eps(q1, [p1], [p2], [p3]));
multiply cuti(1, 1);
repeat id cuti([x]?, [y]?) * cuti([z]?, [w]?) = cuti([x] * [z], [y] * [w]);

* need the CUT*[x] for constructing the local terms
id cuti([x]?, [y]?) = CUT*[x] + cuti([y]);

#endif

*----------------------------------------------------------------------

#if `HaveFermions' == 1

b g_, `Tensors';
.sort

keep brackets;

id g6_([i]?) = 2*GA(6);
id g7_([i]?) = 2*GA(7);
repeat id GA(?g) * g_([i]?, [mu]?) = GA(?g, [mu]);

#call DiracSimplify

#endif

*----------------------------------------------------------------------

b `Patterns', Times, NoExpand, Den,
  SumOver, `SUNObjs',
  `LoopInt', q1.q1, GA, Spinor, i_;
.sort

totensor q1, NUM;

b NUM, pave, cuti, `LoopInt', `Tensors', CUT;
.sort

keep brackets;

*----------------------------------------------------------------------

#if `Dim' == 4
* add local terms for dimred/CDR

id NUM([mu]?, [nu]?, [rho]?, [sig]?) * D0i(?a) =
  NUM([mu], [nu], [rho], [sig]) * D0i(?a) -
  5/144 * TMP(e_([mu], [nu], [rho], [sig]), [mu], [nu], [rho], [sig]) +
  1/8 * distrib_(1, 2, EQ, NEQ, [mu], [nu], [rho], [sig]);

id EQ([mu]?, [mu]?) = 1;
id NEQ([mu]?, [mu]?) = 0;
id NEQ([mu]?, [nu]?) = d_([mu], [nu]);

id NUM([mu]?, [nu]?, [rho]?) * C0i([p1]?, [p2]?, ?a) =
  NUM([mu], [nu], [rho]) * C0i([p1], [p2], ?a) +
  1/36 * TMP(e_([mu], [nu], [rho]), [mu], [nu], [rho], [p1] + [p2]);

id TMP(0, [mu]?, [nu]?, [rho]?, [sig]?) = dd_([mu], [nu], [rho], [sig]);
id TMP(?a) = 0;

id NUM([mu]?, [mu]?) * C0i(?a) = NUM([mu], [mu]) * C0i(?a) - 1/2;

id NUM([mu]?, [mu]?) * B0i([p1]?, [m1]?, [m2]?) =
  A0i([m2]) + [m1]*B0i([p1], [m1], [m2]);

id CUT * [f]?LOOPINT(?a) = 0;
id CUT = 1;

#endif


* decompose into Lorentz-covariant tensors

* The following statement introduces the g_{\mu\nu}'s in a smart way.
* Lifted from: S.A. Larin, T. van Ritbergen, and J.A.M. Vermaseren,
* The optimization of a huge FORM program,
* in: Proceedings Oberammergau 1993, ISBN 9-810-21699-8.

id NUM(?b) = sum_(DUMMY, 0, nargs_(?b), 2,
  pave(0)^DUMMY * distrib_(1, DUMMY, dd_, NUM, ?b));

id B0i([p1]?, ?a) =
  TMP(pave(1)*[p1]) *
  B0i(MOM([p1]), ?a);

id C0i([p1]?, [p2]?, ?a) =
  TMP(pave(1)*[p1] + pave(2)*[p2]) *
  C0i(MOM([p1]), MOM([p2] - [p1]), MOM([p2]), ?a);

id D0i([p1]?, [p2]?, [p3]?, ?a) =
  TMP(pave(1)*[p1] + pave(2)*[p2] + pave(3)*[p3]) *
  D0i(MOM([p1]), MOM([p2] - [p1]), MOM([p3] - [p2]), MOM([p3]),
      MOM([p2]), MOM([p3] - [p1]), ?a);

id E0i([p1]?, [p2]?, [p3]?, [p4]?, ?a) =
  TMP(pave(1)*[p1] + pave(2)*[p2] + pave(3)*[p3] + pave(4)*[p4]) *
  E0i(MOM([p1]), MOM([p2] - [p1]), MOM([p3] - [p2]), MOM([p4] - [p3]), MOM([p4]),
      MOM([p2]), MOM([p3] - [p1]), MOM([p4] - [p2]), MOM([p3]), MOM([p1] - [p4]), ?a);

id F0i([p1]?, [p2]?, [p3]?, [p4]?, [p5]?, ?a) =
  TMP(pave(1)*[p1] + pave(2)*[p2] + pave(3)*[p3] + pave(4)*[p4] + pave(5)*[p5]) *
  F0i(MOM([p1]), MOM([p2] - [p1]), MOM([p3] - [p2]), MOM([p4] - [p3]), MOM([p5] - [p4]), MOM([p5]),
      MOM([p2]), MOM([p3] - [p1]), MOM([p4] - [p2]), MOM([p5] - [p3]), MOM([p4]), MOM([p1] - [p5]),
      MOM([p3]), MOM([p4] - [p1]), MOM([p5] - [p2]), ?a);

repeat id TMP([p1]?) * NUM([mu]?, ?a) = d_([p1], [mu]) * NUM(?a) * TMP([p1]);

id TMP(?a) = 1;
id NUM() = 1;

id A0i([m1]?) * NUM(?a) = 0;

b pave, cuti, `LoopInt';
.sort

keep brackets;

id cuti(?b) * [f]?LOOPINT?CUTINT(?a) = cuti([f](?b), ?a);

chainin pave;
id pave(?b) * [f]?LOOPINT(?a) = pave([f](?b), ?a);

symm B0i 2, 3;
id [f]?LOOPINT(?a) = pave([f](0), ?a);

argument pave, cuti;
#call MomSquare
endargument;

*----------------------------------------------------------------------

#if `HaveFermions' == 1
* Dirac algebra on open fermion chains again

#if "`FermionChains'" == "Weyl"

.sort

* Chisholm's identity backwards to get rid of all Eps
*repeat id GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [rho]?, [LA]?) =
*  1/4 * sign_([om]) * sign_(nargs_(?a)) * (
*    GA([om], ?a, [mu], [nu], [rho], ?b) -
*    GA([om], ?a, [rho], [nu], [mu], ?b) );
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [rho]?, [LA]?) =
    sign_([om]) * sign_(nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [rho], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [rho], ?b) +
      d_([mu], [rho]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [rho]) * GA([om], ?a, [mu], ?b) );
endrepeat;

#elseif `Dim' == 4

.sort

* this is Chisholm's identity:
repeat;
  once GA([om]?, [mu]?, [nu]?, [rho]?, ?a) =
    sign_([om]) * GA([om], DUMMY, ?a) * e_([mu], [nu], [rho], DUMMY) +
    d_([mu], [nu]) * GA([om], [rho], ?a) -
    d_([mu], [rho]) * GA([om], [nu], ?a) +
    d_([nu], [rho]) * GA([om], [mu], ?a);
  sum DUMMY;
endrepeat;

contract 0;

#endif

b `Tensors';
.sort

keep brackets;

#call DiracSimplify

#endif

*----------------------------------------------------------------------

#if `Dim' == D
* add local terms for dimreg

b D, Dminus4, pave, CUT;
.sort

keep brackets;

id D = Dminus4 + 4;
id Dminus4 * pave(A0i(0), [m1]?) = -2*[m1];
id Dminus4 * pave(A0i(0,0), [m1]?) = -[m1]^2/2;
id Dminus4 * pave(B0i(0), ?a) = -2;
id Dminus4 * pave(B0i(1), ?a) = 1;
id Dminus4 * pave(B0i(0,0), [k1]?, [m1]?, [m2]?) =
  1/6*([k1] - 3*[m1] - 3*[m2]);
id Dminus4 * pave(B0i(1,1), ?a) = -2/3;
id Dminus4 * pave(B0i(0,0,1), [k1]?, [m1]?, [m2]?) =
  -1/12*([k1] - 2*[m1] - 4*[m2]);
id Dminus4 * pave(B0i(1,1,1), ?a) = 1/2;
id Dminus4 * pave(C0i(0,0), ?a) = -1/2;
id Dminus4 * pave(C0i(0,0,[i]?), ?a) = 1/6;
id Dminus4 * pave(C0i(0,0,0,0), [k1]?, [k2]?, [k1k2]?, [m1]?, [m2]?, [m3]?) =
  1/48*([k1] + [k2] + [k1k2] - 4*([m1] + [m2] + [m3]));
id Dminus4 * pave(C0i(0,0,[i]?,[i]?), ?a) = -1/12;
id Dminus4 * pave(C0i(0,0,[i]?,[j]?), ?a) = -1/24;
id Dminus4 * pave(D0i(0,0,0,0), ?a) = -1/12;
id Dminus4 * pave(D0i(0,0,0,0,[i]?), ?a) = 1/48;
id Dminus4 = 0;

id CUT * pave(?a) = 0;
id CUT = 1;

#endif

*----------------------------------------------------------------------

#call DotSimplify

.sort

id Den([p1]?, [m1]?) = Den(MOM([p1]), [m1]);
argument Den;
#call MomSquare
endargument;

*----------------------------------------------------------------------

#if `HaveFermions' == 1 && "`FermionChains'" == "Weyl"

id Spinor(?a) * GA(?g) * Spinor(?b) =
  WeylChain(DottedSpinor(?a), ?g, Spinor(?b));

id GA(?g) = WeylChain(DottedSpinor(), ?g, Spinor());

repeat;
  once WeylChain([s1]?, [x]?, ?a, [LA]?, ?b) *
       WeylChain([s2]?, [y]?, ?c, [LA]?, ?d) =
    WC(sign_(nargs_(?a, ?c) + [x] + [y]), [s1], [x], ?a) *
    WC(?b) * WC([s2], [y], ?c) * WC(?d);

* Fierz 1: <A|sig_mu|B> <C|sigbar^mu|D> = 2 <A|D> <C|B>
  id ifmatch->2 WC(-1, ?a) * WC(?b) * WC(?c) * WC(?d) =
    2*WeylChain(?a, ?d) * WeylChain(?c, ?b);

* Fierz 2: <A|sig(bar)_mu|B> <C|sig(bar)^mu|D> = 2 <A|eps|C> <B|eps|D>
  id WC(1, ?a) * WC(?b, [s1]?) * WC([s2]?, [x]?, ?c) * WC(?d) =
    2*WeylChain(?a, reverse_(?c), -1, [s2]) *
      WeylChain([s1], 6 + mod_([x] + nargs_(?c, ?b), 2), -1, reverse_(?b), ?d);

  id WeylChain(?a, -1, -1, ?b) = -WeylChain(?a, ?b);

  label 2;

* due to the canonical ordering of the Dirac chains this
* is the only(?) case we need of Fierz on the same chain:
  repeat id WeylChain(?a, [LA]?, [LA]?, ?b) = 4*WeylChain(?a, ?b);
endrepeat;

id WeylChain([s1]?, ?g, [s2]?) = abb(fme(WeylChain([s1], ?g, [s2])))
#if "`Scale'" != "1"
  * MOM(?g)
#endif
  ;

#endif

*----------------------------------------------------------------------

b Times, NoExpand;
moduleoption polyfun=abb;

#if "`InsertAt'" == "Default"
#call DoInsertions
#else
.sort
#endif

#call Const
.sort

collect Times;

argument Times;
#call Neglect
endargument;

argument;
#call Square
endargument;

id Times(0) = 0;

.sort

*----------------------------------------------------------------------
* index handling

repeat;
  once SumOver([i]?, [x]?, Renumber) =
    TMP(DUMMY) * SumOver(DUMMY, [x]) * replace_([i], DUMMY);
  sum DUMMY;
endrepeat;

id IndexEps([i]?, [j]?, [k]?) = EPS([i], [j], [k]);

repeat;
  id EPS([I]?, [J]?, [K]?) * EPS([I]?, [J]?, [K]?) *
    SumOver([I]?, [x]?) * SumOver([J]?, [y]?) * SumOver([K]?, [z]?) =
    6;
  id EPS([I]?, [J]?, [k]?) * EPS([I]?, [J]?, [c]?) *
    SumOver([I]?, [x]?) * SumOver([J]?, [y]?) =
    2*IndexDelta([k], [c]);
  id EPS([I]?, [j]?, [k]?) * EPS([I]?, [b]?, [c]?) *
    SumOver([I]?, [x]?) =
    IndexDelta([j], [b])*IndexDelta([k], [c]) -
    IndexDelta([j], [k])*IndexDelta([b], [c]);
  repeat;
    id IndexDelta([I]?, [I]?) = 1;
    symm IndexDelta;
    once ifmatch->1 IndexDelta([i]?, [J]?) * SumOver([J]?, [x]?) =
      replace_([J], [i]);
    once IndexDelta([I]?, [j]?) * SumOver([I]?, [x]?) =
      replace_([I], [j]);
    label 1;
  endrepeat;
endrepeat;

id IndexDelta([x]?int_, [y]?int_) = delta_([x], [y]);

id TMP([x]?int_) = 1;
repeat id TMP([I]?)^2 = TMP([I]);

#do i = 1, 10
once TMP([I]?) = replace_([I], Ind`i');
#enddo

id EPS([i]?, [j]?, [k]?) = IndexEps([j], [k], [i]);

*----------------------------------------------------------------------

#if `HaveSUN' == 1
* simplification of SU(N) structures

* The algorithm implemented here is an extension of the one given in
* J.A.M. Vermaseren, The use of computer algebra in QCD,
* in: Proceedings Schladming 1996, ISBN 3-540-62478-3.

* The idea is to transform all SU(N) objects to generators, SUNT.
* In the output, only two types of objects can appear:
* - chains of SUNTs (with external colour indices), or
* - traces of SUNTs.
* A chain of SUNTs is denoted by SUNT(a, b, ..., i, j), where
* a, b, ... are gluon indices and i and j are colour indices.
* SUNT(i, j) is the special case of the identity in colour space.
* A trace over SUNTs is marked by both colour indices being zero,
* i.e. SUNT(a, b, ..., 0, 0).

b `SUNObjs';
.sort

keep brackets;

if( count(SUNF, 1) );

  repeat;
    once SUNF(?a, [a]?, [b]?, [c]?, [d]?) =
      SUNF(?a, [a], [b], DUMMY) * SUNF(DUMMY, [c], [d]) * SUNSum(DUMMY);
    sum DUMMY;
  endrepeat;

* f^{abc} = 2 i Tr(T^c T^b T^a - T^a T^b T^c)

  id SUNF([a]?, [b]?, [c]?) =
    2*i_*(SUNT([c], [b], [a], 0, 0) - SUNT([a], [b], [c], 0, 0));

endif;


repeat;
  once SUNT(?a, 0, 0) = SUNT(?a, DUMMY, DUMMY) * SUNSum(DUMMY);
  sum DUMMY;
endrepeat;

repeat;
  once SUNT(?a, [a]?, [b]?, [i]?, [j]?) =
    SUNT(?a, [a], [i], DUMMY) * SUNT([b], DUMMY, [j]) * SUNSum(DUMMY);
  sum DUMMY;
endrepeat;


* T^a_{ij} T^a_{kl} =
*   1/2 (delta_{il} delta_{jk} - 1/N delta_{ij} delta_{kl})

id SUNT([a]?, [i]?, [j]?) * SUNT([a]?, [k]?, [l]?) * SUNSum([a]?, ?a) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);

id SUNTSum([i]?, [j]?, [k]?, [l]?) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);


id SUNEps([i]?, [j]?, [k]?) = EPS([i], [j], [k]);


* cleaning up, step 1: get rid of the deltas

repeat;
  id EPS([I]?, [j]?, [k]?) * EPS([I]?, [b]?, [c]?) *
    SUNSum([I]?, [x]?) = 
    SUNT([j], [b])*SUNT([k], [c]) -
    SUNT([j], [k])*SUNT([b], [c]);
  repeat;
    id SUNT([I]?, [I]?) * SUNSum([I]?, ?a) = `SUNN';
    symm SUNT:2 1, 2;
    once ifmatch->1 SUNT([I]?, [j]?) * SUNSum([I]?, ?a) = replace_([I], [j]);
    once SUNT([i]?, [J]?) * SUNSum([J]?, ?a) = replace_([J], [i]);
    label 1;
  endrepeat;
endrepeat;

id SUNT([x]?int_, [y]?int_) = delta_([x], [y]);
id SUNT([a]?, [i]?, [i]?) * SUNSum([i]?, ?a) = 0;

id EPS([i]?, [j]?, [k]?) = sun(SUNEps([j], [k], [i]));

* cleaning up, step 2: bead up the SUNTs into chains

repeat;
  once SUNT(?a, [a]?, [i]?, [j]?) = TMP(?a, [a], [i], [j]);
  repeat;
    id TMP(?a, [i]?, [j]?) * SUNT(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
    id SUNT(?a, [i]?, [j]?) * TMP(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
  endrepeat;

  id TMP(?a, [i]?, [i]?) * SUNSum([i]?, ?b) = TMP(?a, 0, 0);

* special case of Tr(T^a T^b) = 1/2 delta_{ab}
*  id TMP([a]?, [a]?, 0, 0) = 1/2;
  id TMP([x]?int_, [y]?int_, 0, 0) = 1/2*delta_([x], [y]);

#if "`FermionOrder'" == "Colour"
  id TMP(?a, [i]?COLS[[x]], [j]?COLS[[y]]) =
    sun(SUNT(?a, [i], [j])) * ORD(MOMS[[x]], MOMS[[y]]);
#endif

  id TMP(?a) = sun(SUNT(?a));
endrepeat;

#if "`FermionOrder'" == "Colour"
id SUNT(?a, [i]?COLS[[x]], [j]?COLS[[y]]) =
  sun(SUNT(?a, [i], [j])) * ORD(MOMS[[x]], MOMS[[y]]);
#endif

id SUNT(?a) = sun(SUNT(?a));

repeat id sun([x]?) * sun([y]?) = sun([x] * [y]);

id SUNSum([i]?, [x]?) = [x];

#endif

*----------------------------------------------------------------------

#if `HaveFermions' == 1 && "`FermionChains'" != "Weyl"

#if "`FermionOrder'" != "None"

* apply 4D Fierz identities

#redefine Dim "4"

#if "`FermionOrder'" != "Fierz"

b `Tensors', ORD;
.sort

keep brackets;

#if "`FermionOrder'" == "Colour"

#elseif "`FermionOrder'" == "Automatic"

* use lexicographical ordering
id Spinor([p1]?, ?a) = ORD([p1]) * Spinor([p1], ?a);
chainin ORD;

#else

#define order "0"
#do i = {`FermionOrder'}
#redefine order "`order',k{`i'}"
#enddo

multiply ORD(`order');
id ORD([i]?, ?r) = ORD(?r);

#endif

id Spinor(?a) * GA(?g) * Spinor(?b) =
  CH(Spinor(?a) * GA(?g) * Spinor(?b));

while( count(ORD, 1) );
* charge conjugation to get first spinor in front
  once ORD([p1]?, ?r) * CH(Spinor([p2]?, [m2]?, [s2]?) *
                           GA([om]?, ?g) *
                           Spinor([p1]?, [m1]?, [s1]?)) =
    ORD([p1], ?r) *
    CH(Spinor([p1], [m1], -[s1]) *
       GR(sign_(nargs_(?g)), [om], reverse_(?g)) *
       Spinor([p2], [m2], -[s2]));

* charge conjugation to get second spinor in back
  once ORD([p3]?, [p2]?, ?r) * CH(Spinor([p2]?, [m2]?, [s2]?) *
                                  GA([om]?, ?g) *
                                  Spinor([p1]?, [m1]?, [s1]?)) =
    ORD([p3], [p2], ?r) *
    CH(Spinor([p1], [m1], -[s1]) *
       GR(sign_(nargs_(?g)), [om], reverse_(?g)) *
       Spinor([p2], [m2], -[s2]));

  argument CH;
    id GR(-1, [om]?{6,7}[[x]], ?g) = GA({7,6}[[x]], ?g);
    id GR(1, ?g) = -GA(?g);
  endargument;

* Fierz to get second spinor together with first
  once ORD([p1]?, [p2]?, ?r) *
    CH(Spinor([p1]?, ?a) * GA([omA]?, ?A) * Spinor([p4]?, ?d)) *
    CH(Spinor([p3]?, ?c) * GA([omB]?, ?B) * Spinor([p2]?, ?b)) =
    ORD([p2], [p4], ?r) * sum_([c], 1, 5, sum_([d], 1, 5,
      DUAL([d], MUD, NUD) * CHI([omA])*g_(1, ?A) *
      DUAL([c], MUC, NUC) * CHI([omB])*g_(1, ?B) *
      CH(Spinor([p1], ?a) * BASIS([d], MUD, NUD) * Spinor([p2], ?b)) *
      CH(Spinor([p3], ?c) * BASIS([c], MUC, NUC) * Spinor([p4], ?d)) ));
  trace4, 1;
  sum MUC, NUC, MUD, NUD;

  id ORD([p1]?, [p2]?, ?r) = ORD(?r);
  id ORD() = 1;
endwhile;

#call FierzPost

#else

b `Tensors';
.sort

keep brackets;

#call FierzPre
#call FierzPost

b `Tensors';
.sort

keep brackets;

#call FierzPre
#call FierzPost

#endif

.sort

#endif

id GA([om]?, [mu]?, [nu]?) =
  GA([om]) * (d_([mu], [nu]) + Sigma([mu], [nu]));

antisymm Sigma;

#call DiracOrder
#call kikj

id D = 4;

b `Tensors';
.sort

keep brackets;

#if "`FermionChains'" == "VA"
id GA([om]?, ?a) = GA(1, ?a)/2 + sign_([om]) * GA(5, ?a)/2;
#endif

id Spinor(?a) * GA(?g) * Sigma(?s) * Spinor(?b) =
  ABB(1, DiracChain(Spinor(?a), ?g, Sigma(?s), Spinor(?b)), ?g, ?s);

id Spinor(?a) * GA(?g) * Spinor(?b) =
  ABB(1, DiracChain(Spinor(?a), ?g, Spinor(?b)), ?g);

id GA(?g) * Sigma(?s) = ABB(1, DiracChain(?g, Sigma(?s)), ?g, ?s);

id GA(?g) = ABB(1, DiracChain(?g), ?g);

id Spinor(?a) * Spinor(?b) =
  ABB(1, DiracChain(Spinor(?a), Spinor(?b)));

#endif

*----------------------------------------------------------------------

.sort

id [p1]?.[p2]? = ABB(0, [p1].[p2], [p1], [p2]);

id e_([mu]?, [nu]?, [rho]?, [sig]?) =
  ABB(0, e_([mu], [nu], [rho], [sig]), [mu], [nu], [rho], [sig]);

id [t]?(?a) = ABB(0, [t](?a), ?a);

repeat;
  once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b) *
       ABB([s2]?, [y]?, ?c, [mu]?, ?d) =
    ABB([s1] + [s2], [x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], DUMMY);
  also ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    ABB([s1], [x], ?a, ?b, ?c) * replace_([mu], DUMMY);
  sum DUMMY;
endrepeat;

id ABB(0, [x]?, ?a) = abb([x])
#if "`Scale'" != "1"
  * MOM(?a)
#endif
  ;
id ABB([i]?, [x]?, ?a) = fme([x])
#if "`Scale'" != "1"
  * MOM(?a)
#endif
  ;

#if "`Scale'" != "1"
$pow = 0;
argument MOM;
  $pow = $pow + count_(<k1,1>,...,<k`Legs',1>);
endargument;
id MOM(?a) = 1;
if( $pow != 0 ) multiply pow(`Scale', $pow/2);
#endif

id i_ = Times(i_);

argument Times;
id i_ = I;
endargument;

*----------------------------------------------------------------------

moduleoption polyfun=Times;
.sort

normalize Times;

moduleoption polyfun=abb;
.sort

normalize abb;
id abb(1) = 1;

* the Mat(...) are kept at the almost-outermost level (only SumOver
* comes before), i.e. the amplitude is of the form Sum[c[i] Mat[i], i]
* -> need this for the calculation of the squared amplitude

#if "`FermionChains'" != "Weyl"
repeat id fme([x]?) * fme([y]?) = fme([x] * [y]);
id fme([x]?) = Mat(fme([x]));
#endif

id sun([x]?) = Mat(sun([x]));

repeat id Mat([x]?) * Mat([y]?) = Mat([x] * [y]);

b SumOver, Mat, Den, pave, cuti, abb, Times;
.sort

collect Simplify;

moduleoption polyfun=Simplify;
.sort

normalize Simplify;

#call Factor(Simplify)
#call InvSimplify(Simplify)

.sort

cf rest;

id Simplify(?a) = rest(Simplify(?a));
id abb(?a) = rest(abb(?a));

moduleoption polyfun=rest;
.sort

normalize rest;
id rest(1) = 1;
id rest(?a) = dum_(?a);

moduleoption polyfun=Times;
.sort

factarg Times;

b SumOver, Mat, Den, pave, cuti;
print;

.end

