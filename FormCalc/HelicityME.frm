* HelicityME.frm
* the FORM part of the HelicityME function
* this file is part of FormCalc
* last modified 7 Jun 16 th


#procedure Emit
id DiracChain(Spinor(?p), [x]?pos_, ?g, Spinor(?q)) =
  DiracChain(Spinor(?p), 1, [x], ?g, Spinor(?q));

* un-antisymmetrize the Dirac chains if necessary
also DiracChain(Spinor(?p), [x]?, ?g, Spinor(?q)) =
  Spinor(?p) * GF(-[x]) * sum_(KK, 0, nargs_(?g), 2,
    sign_(KK/2) * distrib_(-1, KK, DD, GA, ?g)) * Spinor(?q);

id DD() = 1;
id DD([mu]?, [nu]?) = d_([mu], [nu]);
repeat;
  once DD(?a) = g_(1, ?a)/4;
  trace4, 1;
endrepeat;

id Spinor(?p) * GF([om]?) * GA(?g) * Spinor(?q) =
  DiracChain(Spinor(?p), 1, [om], ?g, Spinor(?q));

* The explicit 1 above is make each chirality projector count as
* two, such that sign_(nargs_(.)) effectively ignores the projector.

.sort

repeat;
  id DiracChain(?a, Spinor(?p)) * DiracChain(Spinor(?p), ?b) =
    DiracChain(?a, RHO(?p), ?b);

* If the spinors at the ends don't match directly, i.e.
*   <s2| g1 g2... |s1> <s2| ... |>,
* we use charge conjugation to reverse the first chain to have the
* |s2>'s side by side for substituting the projector (|s2><s2|).
* Inserting 1 = C C^-1 results in
*   <s2|C (C^-1 g1 C) (C^-1 g2 C) ... C^-1|s1>
*   = <anti-s2| (-g1)^T (-g2)^T ... |anti-s1>
*   = (-1)^(# gammas) <anti-s2| (... g2 g1)^T |anti-s1>
*   = (-1)^(# gammas) <anti-s1| ... g2 g1 |anti-s2>
* Thus follow the rules:
*   a) reverse the chain and exchange u <-> v,
*   b) gamma_mu -> -gamma_mu.
*   c) add a global minus sign to compensate for the
*      change in the permutation of the external fermions.
* For more details see the Denner/Eck/Hahn/Kueblbeck paper.
* Note that RHO and RHOC are counted as gamma matrices towards
* the overall sign; this is corrected in the RHOC substitution
* later.

  id DiracChain(Spinor([p2]?, ?p), ?a, Spinor([p1]?, [m1]?, [s1]?)) *
     DiracChain(Spinor([p2]?, ?q), ?b) =
    -sign_(nargs_(?a)) *
    DiracChain(Spinor([p1], [m1], -[s1]),
      reverse_(?a)*replace_(RHO, RHOC, RHOC, RHO),
      RHO([p2], ?q), ?b);
endrepeat;

repeat;
  once DiracChain(Spinor(?p), ?g, Spinor(?p)) = RHO(?p) * GF(?g);
  chainout GF;

  id GF([mu]?) = CHI([mu]);
  id CHI([mu]?) = g_(1, [mu]);
  id GF([x]?) = [x];

  argument RHO, RHOC;
#call Neglect
  endargument;

  id RHO([p1]?MOMS[[x]], 0, [s1]?) =
    [s1]/4*(g6_(1)*HEL([x], [s1]) - g7_(1)*HEL([x], -[s1])) *
      g_(1, [p1]);

  id RHOC([p1]?MOMS[[x]], 0, [s1]?) =
    g_(1, [p1]) *
      [s1]/4*(g6_(1)*HEL([x], [s1]) - g7_(1)*HEL([x], -[s1]));

  id RHO([p1]?MOMS[[x]], [m1]?, [s1]?) =
    (g_(1) + HEL([x], 0)*g_(1, 5_, EPSS[[x]]))/2 *
      (g_(1, [p1]) + [s1]*[m1]*g_(1));

  id RHOC([p1]?MOMS[[x]], [m1]?, [s1]?) =
    (g_(1, [p1]) - [s1]*[m1]*g_(1)) *
      (g_(1) - HEL([x], 0)*g_(1, EPSS[[x]], 5_))/2;

  trace4, 1;
endrepeat;

id D = Dminus4 + 4;

contract;
id D = Dminus4Eps + 4;

#call eiei
#call DotSimplify
#call Abbreviate

b helM;
.sort

#call CollectTerms

.sort

#write "%X"

b helM;
print;
.end
#endprocedure

***********************************************************************

cf DiracChain, Mat;
s D, Dminus4, Dminus4Eps, `Invariants';

i KK;
t DD, Pol;
nt GA;
f RHO, RHOC, GF;
cf TMP, ABB;
auto s ARG;
s TAG, ETAG, QTAG;
set MOMS: k1,...,k`Legs';
set EPSS: e1,...,e`Legs';
set INVS: `Invariants';

i [om], [mu], [nu], [ro], [si], [i];
v [p1], [p2];
s [m1], [m2], [s1], [s2], [x], [y], [n];
t [t];

extrasymbols array subM;
cf abbM, fermM, dotM, helM, addM, mulM, powM;

ntable CHI(0:7);
fill CHI(0) = g_(1);
fill CHI(1) = g_(1);
fill CHI(4) = -g5_(1);
fill CHI(5) = g5_(1);
fill CHI(6) = g6_(1)/2;
fill CHI(7) = g7_(1)/2;

