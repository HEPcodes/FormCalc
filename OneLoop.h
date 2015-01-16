* OneLoop.h
* the main simplification of FormCalc amplitudes in FORM 3
* this file is part of FormCalc
* last modified 12 Dec 01 th


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

while(count(again, 1));
id again = 1;
#call DiracOrder
endwhile;

b ga, Spinor;
.sort
keep brackets;
#endprocedure


#procedure DiracOrder
#if `OnShell' == 1
repeat;
  once ga(OM?, ?a, P1?, ?b) * Spinor(P1?, M1?, X1?) =
    ( 2*g4(OM, ?a, P1) * distrib_(-1, 1, g4, g4, ?b) +
      sign_(nargs_(?b))*X1*M1 * ga(OM, ?a, ?b) ) * Spinor(P1, M1, X1);
  id g4(?a, P1?) * g4(J1?) * g4(?b) =
    d_(P1, J1) * ga(?a, ?b) * again;
endrepeat;
repeat;
  once Spinor(P1?, M1?, X1?) * ga(OM?{6,7}[X2], ?a, P1?, ?b) =
    Spinor(P1, M1, X1) * sign_(nargs_(?a)) *
    ( X1*M1 * ga({7,6}[X2], ?a, ?b) -
      2*distrib_(-1, 1, g4, g4, ?a) * g4(P1, OM, ?b) );
  id g4(J1?) * g4(?a) * g4(P1?, OM?, ?b) =
    d_(P1, J1) * ga(OM, ?a, ?b) * again;
endrepeat;
#endif

id ga(OM?, J1?, ?a) = ga(OM) * g2(J1, ?a);
while(count(g2, 1));
  repeat;
    id g2(LA?, LA?, ?a) = d_(LA, LA) * g2(?a);
    id g2(LA?, J1?, LA?, ?a) = (2 - d_(LA, LA)) * g2(J1, ?a);
    id g2(LA?, J1?, J2?, LA?, ?a) = 4*d_(J1, J2) * g2(?a) * again
#if `Dim' == 0
        + [D-4] * g2(J1, J2, ?a)
#endif
       ;
    id g2(LA?, J1?, J2?, J3?, ?b, LA?, ?a) =
      -sign_(nargs_(?b)) * (
        2*g2(J3, J2, J1, ?b, ?a) +
        2*g4(J1, J2, J3) * distrib_(-1, 1, g4, g4, ?b) * g4(?a)
#if `Dim' == 0
        + [D-4] * g2(J1, J2, J3, ?b, ?a)
#endif
      );
    id g4(?a) * g4(J1?) * g4(?b) * g4(?c) = g2(J1, ?a, ?b, ?c);
  endrepeat;
  id g2(J1?, ?a) = g3(J1) * g2(?a);
  id g2 = 1;
endwhile;
repeat;
  id g3(P1?) * g3(P1?) = P1.P1;
  disorder g3(J1?) * g3(J2?) = 2*d_(J1, J2) * again - g3(J2) * g3(J1);
endrepeat;
repeat id ga(?a) * g3(?b) = ga(?a, ?b);
#endprocedure


#procedure MomSquare
#ifdef `MomSum'
id sq(P1?) = sq(P1, `MomSum' + P1, `MomSum' - P1);
id sq(P1?, P2?, P3?) =
  sq(P1, P2, P3, nterms_(P1), nterms_(P2), nterms_(P3));
symm sq (4,1) (5,2) (6,3);
#endif
id sq(P1?, ?a) = P1.P1;
#call kikj
#endprocedure


*** the actual FORM code starts here ***

ab `Const';
.sort

s sqrt2, CW, SW, CW2, SW2;
s EL, Alfa, Alfa2, GS, Alfas, Alfas2, scale;
cf Mat, Eps, DiracChain, fme, sun, SumOver, Delta;
cf num, abb, kin, Den, pave, A0i, B0i, C0i, D0i, E0i;
cf SUNT, SUNTSum, SUNF;
nt ga;

s X1, X2;
i J1, J2, J3, J4, LA, OM;

.global

s M1, M2, X1, X2, [D-4], again;
v P1, P2, P3, P4;
cf PP, sq;
t QQ, eq, neq;
nt g2, g3, g4;

polyfun num;
collect num;

.sort

polyfun;

argument num;
#call TrivialSubst
endargument;
normalize num;

id num(0) = 0;

#if `FermionChains' == 1

b g_, Spinor;
.sort
keep brackets;

id g6_(J1?) = 2*ga(6);
id g7_(J1?) = 2*ga(7);
repeat id ga(?a) * g_(J1?, J2?) = ga(?a, J2);

#if `OnShell' == 1
#call MomConserv(DiracSimplify)
#else
#call DiracSimplify(,)
#endif

#endif

b num, Den, SumOver, A0i, B0i, C0i, D0i, E0i, q1.q1, ga, Spinor;
.sort

collect abb, abb;

* cancel q^2's in the numerator

while( match(q1.q1) );
  once q1.q1 * A0i(M1?) =
    1 + M1 * A0i(M1);
  once q1.q1 * B0i(P1?, M1?, ?a) =
    replace_(q1, q1 - P1) * A0i(?a) +
    M1 * B0i(P1, M1, ?a);
  once q1.q1 * C0i(P1?, P2?, M1?, ?a) =
    replace_(q1, q1 - P1) * B0i(P2 - P1, ?a) +
    M1 * C0i(P1, P2, M1, ?a);
  once q1.q1 * D0i(P1?, P2?, P3?, M1?, ?a) =
    replace_(q1, q1 - P1) * C0i(P2 - P1, P3 - P1, ?a) +
    M1 * D0i(P1, P2, P3, M1, ?a);
  once q1.q1 * E0i(P1?, P2?, P3?, P4?, M1?, ?a) =
    replace_(q1, q1 - P1) * D0i(P2 - P1, P3 - P1, P4 - P1, ?a) +
    M1 * E0i(P1, P2, P3, P4, M1, ?a);
endwhile;

id A0i(0) = 0;

term;
#call MomConserv(DotSimplify)
endterm;

id abb(X1?) = X1;

totensor q1, QQ;

b QQ, A0i, B0i, C0i, D0i, E0i;
.sort
keep brackets;


#if `Dim' == 4
* add local terms for dimred/CDR

id QQ(J1?, J2?, J3?, J4?) * D0i(?a) =
  QQ(J1, J2, J3, J4) * D0i(?a) -
  5/144 * PP(e_(J1, J2, J3, J4), J1, J2, J3, J4) +
  1/8 * distrib_(1, 2, eq, neq, J1, J2, J3, J4);

id eq(J1?, J1?) = 1;
id neq(J1?, J1?) = 0;
id neq(J1?, J2?) = d_(J1, J2);

id QQ(J1?, J2?, J3?) * C0i(P1?, P2?, ?a) =
  QQ(J1, J2, J3) * C0i(P1, P2, ?a) +
  1/36 * PP(e_(J1, J2, J3), J1, J2, J3, P1 + P2);

id PP(0, J1?, J2?, J3?, J4?) = dd_(J1, J2, J3, J4);
id PP(?a) = 0;

id QQ(J1?, J1?) * C0i(?a) = QQ(J1, J1) * C0i(?a) - 1/2;

id QQ(J1?, J1?) * B0i(P1?, M1?, M2?) = A0i(M2) + M1*B0i(P1, M1, M2);

#endif


* decompose into Lorentz-covariant tensors

id QQ(?b) = sum_(X1, 0, nargs_(?b), 2,
  pave(0, 0)^(X1/2) * distrib_(1, X1, dd_, QQ, ?b));

id A0i(M1?) * QQ(?a) = 0;

id B0i(P1?, ?a) = PP(pave(1)*P1) * B0i(sq(P1), ?a);

id C0i(P1?, P2?, ?a) = PP(pave(1)*P1 + pave(2)*P2) *
  C0i(sq(P1), sq(P2 - P1), sq(P2), ?a);

id D0i(P1?, P2?, P3?, ?a) =
  PP(pave(1)*P1 + pave(2)*P2 + pave(3)*P3) *
  D0i(sq(P1), sq(P2 - P1), sq(P3 - P2), sq(P3), sq(P2), sq(P3 - P1), ?a);

id E0i(P1?, P2?, P3?, P4?, ?a) =
  PP(pave(1)*P1 + pave(2)*P2 + pave(3)*P3 + pave(4)*P4) *
  E0i(sq(P1), sq(P2 - P1), sq(P3 - P2), sq(P4 - P3), sq(P4),
      sq(P2), sq(P3 - P1), sq(P4 - P2), ?a);

repeat id PP(P1?) * QQ(J1?, ?a) = d_(P1, J1) * QQ(?a) * PP(P1);

id PP(?a) = 1;
id QQ = 1;

b pave, A0i, B0i, C0i, D0i, E0i;
.sort
keep brackets;

if(count(pave, 1));
  repeat id pave(?a) * pave(?b) = pave(?a, ?b);
  id pave(?b) * PP?{A0i,B0i,C0i,D0i,E0i}(?a) = pave(PP(?b), ?a);
else;
  symm B0i 2, 3;
  id PP?{A0i,B0i,C0i,D0i,E0i}(?a) = pave(PP(0), ?a);
endif;

argument pave;
#call MomSquare
endargument;


#if `FermionChains' == 1
* Dirac algebra on open fermion chains again

#if `OnShell' == 1
#call MomConserv(DiracSimplify)
#else
#call DiracSimplify(,)
#endif

#if `Dim' == 4
* this is Chisholm's identity:
  repeat;
    once ga(OM?, J1?, J2?, J3?, ?a) =
      sign_(OM) * ga(OM, LA, ?a) * e_(J1, J2, J3, LA) +
      d_(J1, J2) * ga(OM, J3, ?a) -
      d_(J1, J3) * ga(OM, J2, ?a) +
      d_(J2, J3) * ga(OM, J1, ?a);
    sum LA;
  endrepeat;
  contract 0;
#endif

#endif


#if `Dim' == 0
* add local terms for dimreg

b D, [D-4], pave;
.sort
keep brackets;

id D = [D-4] + 4;
id [D-4] * pave(A0i(0), M1?) = -2*M1;
id [D-4] * pave(B0i(0), ?a) = -2;
id [D-4] * pave(B0i(1), ?a) = 1;
id [D-4] * pave(B0i(0,0), X1?, M1?, M2?) = X1/6 - M1/2 - M2/2;
id [D-4] * pave(B0i(1,1), ?a) = -2/3;
id [D-4] * pave(C0i(0,0), ?a) = -1/2;
id [D-4] * pave(C0i(0,0,J1?), ?a) = 1/6;
id [D-4] * pave(D0i(0,0,0,0), ?a) = -1/12;
id [D-4] = 0;

#endif


b num, Den, SumOver, pave, ga, Spinor, i_;
.sort

collect abb, abb;

term;
#call MomConserv(DotSimplify)
endterm;

id Den(P1?, M1?) = Den(sq(P1), M1);
argument Den;
#call MomSquare
endargument;

#if `Scaled' == 1
argument abb, ga;
$pow = count_(<k1,1>,...,<k`Legs',1>);
endargument;
multiply scale^$pow;
#endif

id abb(X1?) = X1;


#if `FermionChains' == 1

repeat;
  once e_(J1?, J2?, J3?, J4?) * ga(?a, J1?, ?b) =
    ga(?a, LA, ?b) * fme(Eps(LA, J2, J3, J4)) * eq(J2) * eq(J3) * eq(J4);
  sum LA;
endrepeat;

b ga, Spinor, fme, eq;
.sort
keep brackets;

id ga(OM?, ?a) = ga(OM, ?a) * eq(?a);
repeat id eq(J1?, J2?, ?a) = eq(J1) * eq(J2, ?a);
repeat;
  once eq(J1?) * eq(J1?) = replace_(J1, LA);
  sum LA;
endrepeat;
id eq(?a) = 1;

#if `VADecompose' == 1
id ga(OM?, ?a) = ga(1, ?a)/2 + sign_(OM) * ga(5, ?a)/2;
#endif

id Spinor(?a) * ga(?b) * Spinor(?c) =
  fme(DiracChain(Spinor(?a), ?b, Spinor(?c)));
id Spinor(?a) * Spinor(?b) =
  fme(DiracChain(Spinor(?a), Spinor(?b)));
id ga(?a) = fme(DiracChain(?a));
repeat id fme(X1?) * fme(X2?) = fme(X1 * X2);

#endif

.sort

polyfun abb;

id P1?.P2? = abb(P1.P2);
id e_(J1?, J2?, J3?, J4?) = abb(e_(J1, J2, J3, J4));

b num;
.store

polyfun;

#call Insertions

ab `Const';
.sort

polyfun num;
collect num;

argument num;
#call TrivialSubst
endargument;


repeat;
  repeat once Delta(J1?, LA?) * SumOver(LA?, J2?) = replace_(LA, J1);
  repeat once Delta(LA?, J1?) * SumOver(LA?, J2?) = replace_(LA, J1);
  id Delta(X1?int_, X2?int_) = delta_(X1, X2);
endrepeat;


#if `SUNobjs' == 1
* simplification of SU(N) structures

* The algorithm implemented here is an extension of the one given in
*   J.A.M. Vermaseren, The use of computer algebra in QCD,
*   in: H. Latal, W. Schweiger, Proceedings Schladming 1996,
*       Springer Verlag, ISBN 3-540-62478-3.

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

d `SUNN';
i COL1, COL2, COL3, COL4;
i GLU1, GLU2, GLU3, GLU4, GLU5;
cf bead;

keep brackets;

id bead?{SUNT,SUNTSum,SUNF}(?a) = bead(?a) * sun(?a);
repeat;
  id sun(COL1?, ?a) * SumOver(COL1?, COL2?, X1?) = sun(COL1, ?a);
  id sun(COL1?, ?a) = sun(?a);
endrepeat;
id sun(?a) = 1;


if(count(SUNF, 1));

  repeat;
    once SUNF(?a, GLU1?, GLU2?, GLU3?, GLU4?) =
      SUNF(?a, GLU1, GLU2, GLU5) * SUNF(GLU5, GLU3, GLU4) * SumOver(GLU5);
    sum GLU5;
  endrepeat;

* f^{abc} = 2 i Tr(T^c T^b T^a - T^a T^b T^c)

  id SUNF(GLU1?, GLU2?, GLU3?) =
    2*i_*(SUNT(GLU3, GLU2, GLU1, 0, 0) - SUNT(GLU1, GLU2, GLU3, 0, 0));

endif;


repeat;
  once SUNT(?a, 0, 0) = SUNT(?a, COL1, COL1) * SumOver(COL1);
  sum COL1;
endrepeat;

repeat;
  once SUNT(?a, GLU1?, GLU2?, COL1?, COL2?) =
    SUNT(?a, GLU1, COL1, COL3) * SUNT(GLU2, COL3, COL2) * SumOver(COL3);
  sum COL3;
endrepeat;


* T^a_{ij} T^a_{kl} =
*   1/2 (delta_{il} delta_{jk} - 1/N delta_{ij} delta_{kl})

id SUNT(GLU1?, COL1?, COL2?) * SUNT(GLU1?, COL3?, COL4?) *
    SumOver(GLU1?, ?a) =
  1/2 * SUNT(COL1, COL4) * SUNT(COL2, COL3) -
  1/2/`SUNN' * SUNT(COL1, COL2) * SUNT(COL3, COL4);

id SUNTSum(COL1?, COL2?, COL3?, COL4?) =
  1/2 * SUNT(COL1, COL4) * SUNT(COL2, COL3) -
  1/2/`SUNN' * SUNT(COL1, COL2) * SUNT(COL3, COL4);


* cleaning up, step 1: get rid of the deltas

repeat;
  id SUNT(COL1?, COL1?) * SumOver(COL1?, ?a) = d_(COL1, COL1);
  once SUNT(COL1?, COL2?) * SumOver(COL1?, ?a) = replace_(COL1, COL2);
  id SUNT(COL1?, COL1?) * SumOver(COL1?, ?a) = d_(COL1, COL1);
  once SUNT(COL2?, COL1?) * SumOver(COL1?, ?a) = replace_(COL1, COL2);
endrepeat;

id SUNT(COL1?, COL1?) = 1;
id SUNT(GLU1?, COL1?, COL1?) * SumOver(COL1?, ?a) = 0;

symm SUNT:2 1, 2;

* cleaning up, step 2: bead up the SUNTs into chains

repeat;
  once SUNT(?a, GLU1?, COL1?, COL2?) = bead(?a, GLU1, COL1, COL2);
  repeat;
    id bead(?a, COL1?, COL2?) * SUNT(?b, COL2?, COL3?) *
        SumOver(COL2?, ?c) =
      bead(?a, ?b, COL1, COL3);
    id SUNT(?a, COL1?, COL2?) * bead(?b, COL2?, COL3?) *
        SumOver(COL2?, ?c) =
      bead(?a, ?b, COL1, COL3);
  endrepeat;
  id bead(?a, COL1?, COL1?) * SumOver(COL1?, ?b) = bead(?a, 0, 0);
* special case of Tr(T^a T^b) = 1/2 delta_{ab}
  id bead(GLU1?, GLU1?, 0, 0) = 1/2;
  id bead(?a) = sun(SUNT(?a));
endrepeat;

id SUNT(?a) = sun(SUNT(?a));

repeat id sun(X1?) * sun(X2?) = sun(X1 * X2);

#endif

.sort

polyfun abb;
normalize num;

.sort

polyfun;
normalize abb;
id abb(1) = 1;

* the Mat(...) are kept at the almost-outermost level (only SumOver
* comes before), i.e. the amplitude is of the form Sum[c[i] Mat[i], i]
* -> need this for the calculation of the squared amplitude

id fme(X1?) = Mat(fme(X1));
id sun(X1?) = Mat(sun(X1));
repeat id Mat(X1?) * Mat(X2?) = Mat(X1 * X2);

b SumOver, Mat, Den, pave, i_, num, abb;
.sort

polyfun kin;
collect kin;

.sort

cf rest;
polyfun rest;

normalize kin;
id kin(1) = 1;
id kin(?a) = rest(kin(?a));
id abb(?a) = rest(abb(?a));

.sort

polyfun num;

normalize rest;
id rest(1) = 1;
id rest(?a) = dum_(?a);

print;
b SumOver, Mat, Den, pave, i_;

.end

