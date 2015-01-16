* CalcFeynAmp.frm
* the FORM part of the CalcFeynAmp function
* this file is part of FormCalc
* last modified 6 Mar 03 th


#procedure DotSimplify(momsubst, moresimp)
argument abb;
`momsubst'
#call eiki
`moresimp'
endargument;
sort;
#endprocedure


#procedure DiracSimplify(momsubst, moresimp)
`momsubst'
#call DiracOrder

while( count(TAG, 1) );
id TAG = 1;
#call DiracOrder
endwhile;

b GA, Spinor;
.sort

keep brackets;
#endprocedure


#procedure DiracOrder
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
#if `Dim' == 0
        + Dminus4 * GB([mu], [nu], ?a)
#endif
       ;
    also GB([LA]?, [mu]?, [nu]?, [rho]?, ?b, [LA]?, ?a) =
      -sign_(nargs_(?b)) * (
        2*GB([rho], [nu], [mu], ?b, ?a) +
        2*GD([mu], [nu], [rho]) * distrib_(-1, 1, GD, GD, ?b) * GD(?a)
#if `Dim' == 0
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

#endprocedure


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


#procedure Factor(foo)
factarg `foo';
repeat id `foo'([x]?, [y]?, ?a) = `foo'([x]) * `foo'([y], ?a);
id `foo'([x]?number_) = [x];
id `foo'([x]?symbol_) = [x];
#endprocedure


*** the main program starts here ***

ab `Const';
.sort

* variables appearing in the CalcFeynAmp input and output
s sqrt2, Pi, CW, SW, CW2, SW2;
s EL, Alfa, Alfa2, GS, Alfas, Alfas2;
s External;
cf Mat, Eps, DiracChain, SumOver, Delta, NoExpand;
cf Times, Simplify, Den, A0i, B0i, C0i, D0i, E0i;
cf SUNT, SUNTSum, SUNF;

* variables that make it into Mma but don't appear in the output
s scale;
cf abb, pave, fme, sun;

* variables internal to FORM
s TAG;
i DUMMY;
cf TMP, MOM;
t NUM, EQ, NEQ;
nt GA, GB, GC, GD;
auto s FC;

* patterns
s [x], [y], [m1], [m2], [s1];
i [mu], [nu], [rho], [sig], [om], [LA];
i [i], [j], [k], [l];
i [a], [b], [c], [d];
v [p1], [p2], [p3], [p4];
cf [f];

.global


collect Times;

moduleoption polyfun=Times;
.sort

argument Times;
#call TrivialSubst
endargument;
normalize Times;

id Times(0) = 0;

#if `FermionChains' == 1

b g_, Spinor;
.sort

keep brackets;

id g6_([om]?) = 2*GA(6);
id g7_([om]?) = 2*GA(7);
repeat id GA(?a) * g_([om]?, [mu]?) = GA(?a, [mu]);

#if `OnShell' == 1
#call MomSimplify(DiracSimplify)
#else
#call DiracSimplify(,)
#endif

#endif

b Times, Den, SumOver, A0i, B0i, C0i, D0i, E0i, q1.q1, Spinor, GA;
.sort

collect abb, abb;

* cancel q^2's in the numerator

while( match(q1.q1) );
  once q1.q1 * A0i([m1]?) =
    1 + [m1] * A0i([m1]);
  once q1.q1 * B0i([p1]?, [m1]?, ?a) =
    replace_(q1, q1 - [p1]) * A0i(?a) +
    [m1] * B0i([p1], [m1], ?a);
  once q1.q1 * C0i([p1]?, [p2]?, [m1]?, ?a) =
    replace_(q1, q1 - [p1]) * B0i([p2] - [p1], ?a) +
    [m1] * C0i([p1], [p2], [m1], ?a);
  once q1.q1 * D0i([p1]?, [p2]?, [p3]?, [m1]?, ?a) =
    replace_(q1, q1 - [p1]) * C0i([p2] - [p1], [p3] - [p1], ?a) +
    [m1] * D0i([p1], [p2], [p3], [m1], ?a);
  once q1.q1 * E0i([p1]?, [p2]?, [p3]?, [p4]?, [m1]?, ?a) =
    replace_(q1, q1 - [p1]) * D0i([p2] - [p1], [p3] - [p1], [p4] - [p1], ?a) +
    [m1] * E0i([p1], [p2], [p3], [p4], [m1], ?a);
endwhile;

id A0i(0) = 0;

term;
#call MomSimplify(DotSimplify)
endterm;

id abb([x]?) = [x];

totensor q1, NUM;

b NUM, A0i, B0i, C0i, D0i, E0i;
.sort

keep brackets;

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

#endif


* decompose into Lorentz-covariant tensors

* The following statement introduces the g_{\mu\nu}'s in a tricky way.
* Lifted from: S.A. Larin, T. van Ritbergen, and J.A.M. Vermaseren,
* The optimization of a huge FORM program,
* in: Proceedings Oberammergau 1993, ISBN 9-810-21699-8.

id NUM(?b) = sum_(DUMMY, 0, nargs_(?b), 2,
  pave(0)^DUMMY * distrib_(1, DUMMY, dd_, NUM, ?b));

id A0i([m1]?) * NUM(?a) = 0;

id B0i([p1]?, ?a) = TMP(pave(1)*[p1]) * B0i(MOM([p1]), ?a);

id C0i([p1]?, [p2]?, ?a) = TMP(pave(1)*[p1] + pave(2)*[p2]) *
  C0i(MOM([p1]), MOM([p2] - [p1]), MOM([p2]), ?a);

id D0i([p1]?, [p2]?, [p3]?, ?a) =
  TMP(pave(1)*[p1] + pave(2)*[p2] + pave(3)*[p3]) *
  D0i(MOM([p1]), MOM([p2] - [p1]), MOM([p3] - [p2]), MOM([p3]),
      MOM([p2]), MOM([p3] - [p1]), ?a);

id E0i([p1]?, [p2]?, [p3]?, [p4]?, ?a) =
  TMP(pave(1)*[p1] + pave(2)*[p2] + pave(3)*[p3] + pave(4)*[p4]) *
  E0i(MOM([p1]), MOM([p2] - [p1]), MOM([p3] - [p2]), MOM([p4] - [p3]), MOM([p4]),
      MOM([p2]), MOM([p3] - [p1]), MOM([p4] - [p2]), MOM([p3]), MOM([p1] - [p4]), ?a);

repeat id TMP([p1]?) * NUM([mu]?, ?a) = d_([p1], [mu]) * NUM(?a) * TMP([p1]);

id TMP(?a) = 1;
id NUM() = 1;

b pave, A0i, B0i, C0i, D0i, E0i;
.sort

keep brackets;

if( count(pave, 1) );
  repeat id pave(?a) * pave(?b) = pave(?a, ?b);
  id pave(?b) * [f]?{A0i,B0i,C0i,D0i,E0i}(?a) = pave([f](?b), ?a);
else;
  symm B0i 2, 3;
  id [f]?{A0i,B0i,C0i,D0i,E0i}(?a) = pave([f](0), ?a);
endif;

argument pave;
#call MomSquare
endargument;


#if `FermionChains' == 1
* Dirac algebra on open fermion chains again

#if `OnShell' == 1
#call MomSimplify(DiracSimplify)
#else
#call DiracSimplify(,)
#endif

#if `Dim' == 4
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

#endif


#if `Dim' == 0
* add local terms for dimreg

b D, Dminus4, pave;
.sort

keep brackets;

id D = Dminus4 + 4;
id Dminus4 * pave(A0i(0), [m1]?) = -2*[m1];
id Dminus4 * pave(B0i(0), ?a) = -2;
id Dminus4 * pave(B0i(1), ?a) = 1;
id Dminus4 * pave(B0i(0,0), [x]?, [m1]?, [m2]?) = [x]/6 - [m1]/2 - [m2]/2;
id Dminus4 * pave(B0i(1,1), ?a) = -2/3;
id Dminus4 * pave(C0i(0,0), ?a) = -1/2;
id Dminus4 * pave(C0i(0,0,[i]?), ?a) = 1/6;
id Dminus4 * pave(D0i(0,0,0,0), ?a) = -1/12;
id Dminus4 = 0;

#endif


b Times, Den, SumOver, pave, GA, Spinor, i_;
.sort

collect abb, abb;

term;
#call MomSimplify(DotSimplify)
endterm;

id Den([p1]?, [m1]?) = Den(MOM([p1]), [m1]);
argument Den;
#call MomSquare
endargument;

#if `Scaled' == 1
argument abb, GA;
$pow = count_(<k1,1>,...,<k`Legs',1>);
endargument;
multiply scale^$pow;
#endif

id abb([x]?) = [x];


#if `FermionChains' == 1

repeat;
  once e_([LA]?, [mu]?, [nu]?, [rho]?) * GA(?a, [LA]?, ?b) =
    GA(?a, DUMMY, ?b) * fme(Eps(DUMMY, [mu], [nu], [rho])) *
    EQ([mu]) * EQ([nu]) * EQ([rho]);
  sum DUMMY;
endrepeat;

b GA, Spinor, fme, EQ;
.sort

keep brackets;

id GA([om]?, ?a) = GA([om], ?a) * EQ(?a);
repeat id EQ([mu]?, [nu]?, ?a) = EQ([mu]) * EQ([nu], ?a);
repeat;
  once EQ([mu]?) * EQ([mu]?) = replace_([mu], DUMMY);
  sum DUMMY;
endrepeat;
id EQ(?a) = 1;

#if `VADecompose' == 1
id GA([om]?, ?a) = GA(1, ?a)/2 + sign_([om]) * GA(5, ?a)/2;
#endif

id Spinor(?a) * GA(?b) * Spinor(?c) =
  fme(DiracChain(Spinor(?a), ?b, Spinor(?c)));
id Spinor(?a) * Spinor(?b) =
  fme(DiracChain(Spinor(?a), Spinor(?b)));
id GA(?a) = fme(DiracChain(?a));
repeat id fme([x]?) * fme([y]?) = fme([x] * [y]);

#endif

.sort

id [p1]?.[p2]? = abb([p1].[p2]);
id e_([mu]?, [nu]?, [rho]?, [sig]?) = abb(e_([mu], [nu], [rho], [sig]));

b Times;
moduleoption polyfun=abb;
.store

#call Insertions

ab `Const';
.sort

collect Times;

argument Times;
#call TrivialSubst
endargument;

.sort

repeat;
  repeat once Delta([mu]?, [LA]?) * SumOver([LA]?, [nu]?) = replace_([LA], [mu]);
  repeat once Delta([LA]?, [mu]?) * SumOver([LA]?, [nu]?) = replace_([LA], [mu]);
  id Delta([x]?int_, [y]?int_) = delta_([x], [y]);
endrepeat;


#if `SUNobjs' == 1
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

b SUNT, SUNTSum, SUNF, SumOver;
.sort

keep brackets;

id [f]?{SUNT,SUNTSum,SUNF}(?a) = [f](?a) * TMP(?a);
repeat;
  id TMP([i]?, ?a) * SumOver([i]?, [j]?, External) = TMP(?a);
  also TMP([i]?, ?a) = TMP(?a);
endrepeat;
id TMP() = 1;


if( count(SUNF, 1) );

  repeat;
    once SUNF(?a, [a]?, [b]?, [c]?, [d]?) =
      SUNF(?a, [a], [b], DUMMY) * SUNF(DUMMY, [c], [d]) * SumOver(DUMMY);
    sum DUMMY;
  endrepeat;

* f^{abc} = 2 i Tr(T^c T^b T^a - T^a T^b T^c)

  id SUNF([a]?, [b]?, [c]?) =
    2*i_*(SUNT([c], [b], [a], 0, 0) - SUNT([a], [b], [c], 0, 0));

endif;


repeat;
  once SUNT(?a, 0, 0) = SUNT(?a, DUMMY, DUMMY) * SumOver(DUMMY);
  sum DUMMY;
endrepeat;

repeat;
  once SUNT(?a, [a]?, [b]?, [i]?, [j]?) =
    SUNT(?a, [a], [i], DUMMY) * SUNT([b], DUMMY, [j]) * SumOver(DUMMY);
  sum DUMMY;
endrepeat;


* T^a_{ij} T^a_{kl} =
*   1/2 (delta_{il} delta_{jk} - 1/N delta_{ij} delta_{kl})

id SUNT([a]?, [i]?, [j]?) * SUNT([a]?, [k]?, [l]?) * SumOver([a]?, ?a) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);

id SUNTSum([i]?, [j]?, [k]?, [l]?) =
  1/2 * SUNT([i], [l]) * SUNT([j], [k]) -
  1/2/(`SUNN') * SUNT([i], [j]) * SUNT([k], [l]);


* cleaning up, step 1: get rid of the deltas

repeat;
  id SUNT([i]?, [i]?) * SumOver([i]?, ?a) = `SUNN';
  once SUNT([i]?, [j]?) * SumOver([i]?, ?a) = replace_([i], [j]);
  id SUNT([i]?, [i]?) * SumOver([i]?, ?a) = `SUNN';
  once SUNT([j]?, [i]?) * SumOver([i]?, ?a) = replace_([i], [j]);
endrepeat;

id SUNT([i]?, [i]?) = 1;
id SUNT([a]?, [i]?, [i]?) * SumOver([i]?, ?a) = 0;

symm SUNT:2 1, 2;

* cleaning up, step 2: bead up the SUNTs into chains

repeat;
  once SUNT(?a, [a]?, [i]?, [j]?) = TMP(?a, [a], [i], [j]);
  repeat;
    id TMP(?a, [i]?, [j]?) * SUNT(?b, [j]?, [k]?) * SumOver([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
    id SUNT(?a, [i]?, [j]?) * TMP(?b, [j]?, [k]?) * SumOver([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
  endrepeat;

  id TMP(?a, [i]?, [i]?) * SumOver([i]?, ?b) = TMP(?a, 0, 0);

* special case of Tr(T^a T^b) = 1/2 delta_{ab}
  id TMP([a]?, [a]?, 0, 0) = 1/2;

  id TMP(?a) = sun(SUNT(?a));
endrepeat;

id SUNT(?a) = sun(SUNT(?a));

repeat id sun([x]?) * sun([y]?) = sun([x] * [y]);

#endif

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

id fme([x]?) = Mat(fme([x]));
id sun([x]?) = Mat(sun([x]));
repeat id Mat([x]?) * Mat([y]?) = Mat([x] * [y]);

b SumOver, Mat, Den, pave, i_, Times, abb, NoExpand;
.sort

collect Simplify;

moduleoption polyfun=Simplify;
.sort

normalize Simplify;

#call Factor(Simplify)
#call MandelSimplify(Simplify)

.sort

cf rest;

id Simplify(?a) = rest(Simplify(?a));
id abb(?a) = rest(abb(?a));

moduleoption polyfun=rest;
.sort

polyfun Times;

normalize rest;
id rest(1) = 1;
id rest(?a) = dum_(?a);

factarg Times;

print;
b SumOver, Mat, Den, pave, i_;

.end

