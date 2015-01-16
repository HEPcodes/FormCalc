* CalcFeynAmp.frm
* the FORM part of the CalcFeynAmp function
* this file is part of FormCalc
* last modified 10 Oct 12 th


***********************************************************************

#procedure Contract
repeat once e_([i]?, [j]?, [k]?, [LA]?)*e_([I]?, [J]?, [K]?, [LA]?) =
#if "`Dim'" != "4"
  (1 + Dminus4) *
#endif
  ( d_([i], [I]) * (d_([j], [J])*d_([k], [K]) - d_([j], [K])*d_([k], [J])) +
    d_([i], [J]) * (d_([j], [K])*d_([k], [I]) - d_([j], [I])*d_([k], [K])) +
    d_([i], [K]) * (d_([j], [I])*d_([k], [J]) - d_([j], [J])*d_([k], [I])) );
#endprocedure

***********************************************************************

#procedure DotSimplify(cmd)
#if "`MomElim'" == "Automatic"

.sort

#call eiki

ab `Vectors';
.sort
keep brackets;

#$term = 0;
$term = $term + 1;
mul TMP($term, termsinbracket_);
id TMP([x]?, 1) = 1;

b TMP, `Vectors', e_;
.sort
keep brackets;

id TMP([x]?, ?r) = <TMP([x], 0)> + ... + <TMP([x], `Legs')>;

#do i = 1, `Legs'
#ifdef `k`i''
if( match(TMP([x]?, {`Legs'-`i'})) ) id k`i' = `k`i'';
#endif
#enddo

`cmd'
#call eiki
#call kikj

b TMP;
.sort
keep brackets;

id TMP([x]?, [y]?) = TMP([x], termsinbracket_, [y]);

b TMP;
.sort
keep brackets;

#$term = 0;
$prev = $term;
id TMP([x]?$term, ?a) = deltap_([x], $prev);

#else

b `Vectors', e_;
.sort
keep brackets;

#ifdef `k`MomElim''
#call eiki
id k`MomElim' = `k`MomElim'';
#endif

`cmd'
#call eiki
#call kikj

#endif
#endprocedure

***********************************************************************

#procedure ChainOrder(contr)
label 1;

#if `OnShell' == 1
* Apply Dirac equation to right spinor
repeat;
  id GA([om]?, ?a, [p1]?, ?b) * Spinor([p1]?, [m1]?, [s1]?, ?s1) =
    ( 2*GD([om], ?a, [p1]) * distrib_(-1, 1, GD, GD, ?b) +
      sign_(nargs_(?b)) * [s1]*[m1] * GA([om], ?a, ?b) ) *
    Spinor([p1], [m1], [s1], ?s1);
  id GD(?a, [p1]?) * GD([mu]?) * GD(?b) =
    d_([p1], [mu]) * GA(?a, ?b) * TAG;
endrepeat;

* Apply Dirac equation to left spinor
repeat;
  id Spinor([p1]?, [m1]?, [s1]?, ?s1) * GA([om]?{6,7}[[n]], ?a, [p1]?, ?b) =
    Spinor([p1], [m1], [s1], ?s1) * sign_(nargs_(?a)) *
    ( [s1]*[m1] * GA({7,6}[[n]], ?a, ?b) -
      2*distrib_(-1, 1, GD, GD, ?a) * GD([p1], [om], ?b) );
  id GD([mu]?) * GD(?a) * GD([p1]?, [om]?, ?b) =
    d_([p1], [mu]) * GA([om], ?a, ?b) * TAG;
endrepeat;
#endif

* Eliminate contractions within each Dirac chain using the
* formulas from M. Veltman's Gammatrica [Nucl Phys B319 (1989) 253]
id GA([om]?, [mu]?, ?a) = GA([om]) * GB([mu], ?a);
while( count(GB, 1) );
#if `contr'
  repeat;
    id GB([LA]?, [LA]?, ?a) = d_([LA], [LA]) * GB(?a);
    also GB([LA]?, [mu]?, [LA]?, ?a) = (2 - d_([LA], [LA])) * GB([mu], ?a);
    also GB([LA]?, [mu]?, [nu]?, [LA]?, ?a) =
#if "`Dim'" != "4"
      Dminus4 * GB([mu], [nu], ?a) +
#endif
      4*d_([mu], [nu]) * GB(?a) * TAG;
    also GB([LA]?, [mu]?, [nu]?, [ro]?, ?b, [LA]?, ?a) =
      -sign_(nargs_(?b)) * (
#if "`Dim'" != "4"
        Dminus4 * GB([mu], [nu], [ro], ?b, ?a) +
#endif
        2*GB([ro], [nu], [mu], ?b, ?a) +
        2*GD([mu], [nu], [ro]) * distrib_(-1, 1, GD, GD, ?b) * GD(?a) );
    id GD(?a) * GD([mu]?) * GD(?b) * GD(?c) = GB([mu], ?a, ?b, ?c);
  endrepeat;
#endif
  id GB([mu]?, ?a) = GC([mu]) * GB(?a);
  id GB() = 1;
endwhile;

* Order the gamma matrices canonically
repeat;
  id GC([p1]?) * GC([p1]?) = [p1].[p1];
  disorder GC([mu]?) * GC([nu]?) = 2*d_([mu], [nu]) * TAG - GC([nu]) * GC([mu]);
endrepeat;
chainin GC;
id GA(?a) * GC(?b) = GA(?a, ?b);

id ifmatch->1 TAG = 1;
#endprocedure

***********************************************************************

#procedure ChainSimplify(contr)
#call ChainOrder(`contr')

#if `OnShell' == 1

#do i = 1, `Legs'
#ifdef `k`i''
b `Tensors';
.sort
keep brackets;

id k`i' = `k`i'';
#call ChainOrder(`contr')
#endif
#enddo

#endif
#endprocedure

***********************************************************************

#procedure FierzBefore
#if `Evanescent' == 1
b Spinor;
.sort
keep brackets;

id Spinor(?a) * GA(?g) * Spinor(?b) =
  Spinor(?a) * GA(?g) * Spinor(?b) *
  Evanescent(DiracChain(Spinor(?a), ?g, Spinor(?b)));
repeat id Evanescent([x]?) * Evanescent([y]?) = Evanescent([x] * [y]);
#endif
#endprocedure

***********************************************************************

#procedure FierzPre(ord)
b `Tensors', ORD;
.sort
keep brackets;

id Spinor([p1]?, ?a) * GA(?g) * Spinor([p2]?, ?b) =
  CH(Spinor([p1], ?a), ?g, Spinor([p2], ?b))
#if `ord'
  * ORD([p1]) * ORD([p2]);

repeat id ORD([p1]?) * ORD(?a, [p1]?, ?b) = ORD(?a, [p1], ?b);

chainin ORD
#endif
  ;

* these relations are obtained by Fierzing twice

id CH(Spinor(?a), 6, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), 7, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), 6, Spinor(?b)) *
    CH(Spinor(?c), 7, Spinor(?d));

id CH(Spinor(?a), [om]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), [om]?, [ro]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), [om], Spinor(?b)) *
    CH(Spinor(?c), [om], [ro], Spinor(?d));
also CH(Spinor(?a), [omA]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), [omB]?, [ro]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), [omA], [ro], [mu], Spinor(?b)) *
    CH(Spinor(?c), [omB], [mu], Spinor(?d));

id CH(Spinor(?a), 6, [ro]?, [mu]?, [nu]?, Spinor(?b)) *
   CH(Spinor(?c), 7, [si]?, [mu]?, [nu]?, Spinor(?d)) =
  4*CH(Spinor(?a), 6, [ro], Spinor(?b)) *
    CH(Spinor(?c), 7, [si], Spinor(?d));

id CH(Spinor(?a), [om]?, [mu]?, [nu]?, [ro]?, Spinor(?b)) *
   CH(Spinor(?c), [om]?, [mu]?, [nu]?, [ro]?, Spinor(?d)) =
  16*CH(Spinor(?a), [om], [mu], Spinor(?b)) *
     CH(Spinor(?c), [om], [mu], Spinor(?d));
also CH(Spinor(?a), [omA]?, [mu]?, [nu]?, [ro]?, Spinor(?b)) *
   CH(Spinor(?c), [omB]?, [mu]?, [nu]?, [ro]?, Spinor(?d)) =
  4*CH(Spinor(?a), [omA], [mu], Spinor(?b)) *
    CH(Spinor(?c), [omB], [mu], Spinor(?d));
#endprocedure

***********************************************************************

#procedure FierzPost
id CH([x]?, ?g, [y]?) = [x] * GA(?g) * [y];

id D = 4;
#call ChainSimplify(1)

.sort

id D = 4;
#endprocedure

***********************************************************************
* The following general Fierz identity is from hep-ph/0412245.

#procedure FierzIdentity(lhs, rhs)
once ifnomatch->2 `lhs'
  CH(Spinor([p1]?, ?a), [omA]?, ?A, [s2]?) *
  CH([s1]?, [omB]?, ?B, Spinor([p2]?, ?b)) = `rhs' sum_(JJ, 1, 5,
    CHI([omA])*g_(1, ?A) *
    DUAL(JJ, N21_?, N22_?) *
    CHI([omB])*g_(1, ?B) * sum_(KK, 1, 5,
      DUAL(KK, N31_?, N32_?) *
      Spinor([p1], ?a) * BASIS(KK, 3, N31_?, N32_?) * Spinor([p2], ?b)) *
    [s1] * BASIS(JJ, 2, N21_?, N22_?) * [s2]);

trace4, 1;

id g6_([i]?) = 2*GA(6);
id g7_([i]?) = 2*GA(7);
repeat id GA(?g) * g_([i]?, [mu]?) = GA(?g, [mu]);

* Chisholm's identity backwards to get rid of all e_
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
    sign_([om] + nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [ro], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [ro], ?b) +
      d_([mu], [ro]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [ro]) * GA([om], ?a, [mu], ?b) );
endrepeat;

renumber;

#call ChainOrder(1)
#endprocedure

***********************************************************************

#procedure FierzUnordered
#call FierzPre(0)

#do rep = 1, 1
#call FierzIdentity(,)

label 2;
if( count(CH, 1) ) redefine rep "0";

b CH;
.sort
keep brackets;
#enddo

#call FierzPost
#endprocedure

***********************************************************************

#procedure FierzOrdered
#call FierzPre(1)

#do rep = 1, 1
id ifmatch->2 ORD([p1]?, [p2]?, ?p) *
     CH(Spinor([p1]?, ?m1), ?g, Spinor([p2]?, ?m2)) =
  ORD(?p) * CH(Spinor([p1], ?m1), ?g, Spinor([p2], ?m2));

* charge conjugation to get first spinor in front
* the rules for this are exactly as in HelicityME
id ORD([p1]?, [p3]?, ?p) *
   CH(Spinor([p2]?, [m2]?, [s2]?), [x]?, ?g,
      Spinor([p1]?, [m1]?, [s1]?)) =
  -ORD([p1], [p3], ?p) * sign_(nargs_(?g)) *
  CH(Spinor([p1], [m1], -[s1]),
     (sign_(nargs_(?g))*(2*[x] - 13) + 13)/2, reverse_(?g),
     Spinor([p2], [m2], -[s2]));

* charge conjugation to get second spinor in back
id ORD([p3]?, [p2]?, ?p) *
   CH(Spinor([p2]?, [m2]?, [s2]?), [x]?, ?g,
      Spinor([p1]?, [m1]?, [s1]?)) =
  -ORD([p3], [p2], ?p) * sign_(nargs_(?g)) *
  CH(Spinor([p1], [m1], -[s1]),
     (sign_(nargs_(?g))*(2*[x] - 13) + 13)/2, reverse_(?g),
     Spinor([p2], [m2], -[s2]));

* Fierz to get second spinor together with first
#call FierzIdentity(ORD([p1]?\, [p2]?\, ?p) *, ORD(?p) *)

id Spinor(?a) * GA(?g) * Spinor(?b) = CH(Spinor(?a), ?g, Spinor(?b));

label 2;

id ORD() = 1;
if( count(ORD, 1) ) redefine rep "0";

b CH, ORD;
.sort
keep brackets;
#enddo

#call FierzPost
#endprocedure

***********************************************************************

#procedure DiracFinal

#if `Antisymmetrize' == 0

id GA([om]?, ?g) = CC([om]) * GD(?g);

#else

* introduce antisymmetrized Dirac chains
id GA([om]?, ?g) = CC([om]) *
  sum_(KK, 0, nargs_(?g), 2, distrib_(-1, KK, DD, GD, ?g));

id DD() = 1;
id DD([mu]?, [nu]?) = d_([mu], [nu]);
repeat;
  once DD(?a) = g_(1, ?a)/4;
  trace4, 1;
endrepeat;
id D = 4;

antisymm GD;

#endif

b CC, GD, Evanescent, `Tensors';
.sort
keep brackets;

#if "`FermionChains'" == "VA"
id CC([om]?) = CC(1)/2 + sign_([om]) * CC(5)/2;
#endif

#if `Antisymmetrize' == 1
id CC([w]?) * GD([mu]?, [nu]?, ?r) = CC(-[w]) * GD([mu], [nu], ?r);
#endif

id Spinor(?a) * CC([w]?) * GD(?g) * Spinor(?b) =
  ABB(1, DiracChain(Spinor(?a), [w], ?g, Spinor(?b)), ?g);

id CC([w]?) * GD(?g) = ABB(1, DiracChain([w], ?g), ?g);

repeat id ABB(1, [x]?, ?g) * ABB(1, [y]?, ?h) =
  ABB(1, [x]*[y], ?g, ?h);

#if 1
id ABB(1, DiracChain(?a, -6, [mu]?, [nu]?, ?b) *
          DiracChain(?c, -7, [mu]?, [nu]?, ?d), ?g) = 0;
#endif

#if "`FermionOrder'" != "None" && `Evanescent' == 1
ab ABB;
.sort
collect TMP;
normalize TMP;
id Evanescent([x]?) * TMP([y]?) = Evanescent([x], [y]) + [y];
id TMP([y]?) = [y];
argument Evanescent;
  id ABB([x]?, [y]?, ?r) = [y];
  id DiracChain(?a) = DiracChain(?a) * TMP(?a);
  chainout TMP;
  id TMP([mu]?index_) = ORD([mu]);
  id ORD([mu]?fixed_) = 1;
  id TMP(?a) = 1;
  id ORD([mu])^[n]? = ORD([mu]);
  repeat;
    once ORD([mu]?) = replace_([mu], N100_?);
    renumber;
  endrepeat;
endargument;
id Evanescent([x]?, [x]?) = 0;
#endif

#call Abbreviate
#endprocedure

***********************************************************************

#procedure Abbreviate
#call DotSimplify()

.sort

*----------------------------------------------------------------------

if( count(cutM, 1) );

id q1.[p1]? = qfM(Pair(q1, [p1]));
id e_(q1, [p1]?, [p2]?, [p3]?) = qfM(Eps(q1, [p1], [p2], [p3]));
id abbM(fmeM(WeylChain(?a, q1, ?b))) = qfM(WeylChain(?a, q1, ?b));

repeat id qfM([x]?) * qfM([y]?) = qfM([x] * [y]);

id cutM([n]?, ?a) * HelDelta([x]?) = cutM([n], [x], ?a);
also cutM([n]?, ?a) = cutM([n], 1, ?a);

endif;

id HelDelta(?h) = 1;

.sort

*----------------------------------------------------------------------

id [p1]?.[p2]? = ABB(0, [p1].[p2], [p1], [p2]);

id e_([mu]?, [nu]?, [ro]?, [si]?) =
  ABB(0, e_([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]);

id d_([mu]?, [nu]?) = ABB(0, d_([mu], [nu]), [mu], [nu]);

id [t]?(?a) = ABB(0, [t](?a), ?a);

id [p1]?([mu]?) = ABB(0, [p1]([mu]), [p1]);

repeat;
  once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b) *
       ABB([s2]?, [y]?, ?c, [mu]?, ?d) =
    ABB([s1] + [s2], [x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], N100_?);
  also once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    ABB([s1], [x], ?a, ?b, ?c) * replace_([mu], N100_?);
  renumber;
endrepeat;

id ABB(0, [x]?, ?a) = abbM([x])
#if "`Scale'" != "1"
  * MOM(?a)
#endif
  ;
id ABB([i]?, [x]?, ?a) = fmeM([x])
#if "`Scale'" != "1"
  * MOM(?a)
#endif
  ;

#if "`FermionChains'" != "Weyl"
repeat id fmeM([x]?) * fmeM([y]?) = fmeM([x] * [y]);
argument fmeM;
id e_([mu]?, [nu]?, [ro]?, [si]?) = Eps([mu], [nu], [ro], [si]);
endargument;
id fmeM([x]?) = Mat(fmeM([x]));
#endif

#if "`Scale'" != "1"
chainout MOM;
id MOM([p1]?MOMS) = SCALE;
id MOM(?a) = 1;
id SCALE^[x]? = powM(`Scale', [x]/2);
id powM([x]?, [y]?pos_) = [x]^[y];
#endif

b addM, mulM;
moduleoption polyfun=abbM;
.sort

normalize abbM;

id abbM(1) = 1;

*----------------------------------------------------------------------

#ifdef `PIPES_'
b Den, intM;
.sort

keep brackets;

#setexternal `PIPE1_'

#call MapAll(ToMma)
#toexternal "#"

drop;
.sort

#fromexternal
.global

#call MapAll(FromMma)
#endif
#endprocedure

***********************************************************************

#procedure ToMma(expr)
#toexternal "%E,", `expr'
#endprocedure

***********************************************************************

#procedure FromMma(expr)
G `expr' =
#fromexternal
  ;
#endprocedure

***********************************************************************

#procedure MomSquare
id MOM(0) = 0;

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
id `foo'([x]?, [y]?) = `foo'([x], nterms_([x])*2 - 1, [y], nterms_([y])*2);
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
*id `foo'([x]?symbol_^[n]?) = [x]^[n];
#endprocedure

***********************************************************************

#procedure DoInsertions
.store

#call Insertions
#call Neglect

.sort
#endprocedure


***********************************************************************
*** the main program starts here

#if "`InsertionPolicy'" == "Begin"
#call DoInsertions
#define Inserted
#else
#call Neglect
#endif

#call eiki

#call Const
.sort

*----------------------------------------------------------------------

#define SUNObjs "SUNSum, SUNT, SUNTSum, SUNF, SUNEps"

#define Tensors "Spinor, GA, e_, eT1,...,eT`Legs', eTc1,...,eTc`Legs'"

#define Vectors "k1,...,k`Legs', e1,...,e`Legs', ec1,...,ec`Legs'"

* variables appearing in the CalcFeynAmp input and output
s I, D, Dminus4, Gamma5Test, Finite, MuTilde, MuTildeSq, Renumber;
cf Mat, Den, Pair, Eps, DiracChain, WeylChain, Evanescent, FormSimplify;
cf SumOver, PowerOf, IGram, List;
cf HelDelta, IndexDelta, IndexEps, `SUNObjs', SUNTr;
f Spinor, g5M, g6M, g7M;
i Col1,...,Col`Legs', Ind1,...,Ind10;

* variables that make it into Mma but don't appear in the output
cf addM, mulM, powM, abbM, fmeM, sunM;
cf intM, paveM, cutM, numM, qfM, qcM;
s dm4M;

* patterns
s [x], [y], [z], [w], [n], [h];
s [k1], [k2], [k1k2];
s <[m0]>,...,<[m20]>, [mk];
s <[s0]>,...,<[s20]>;
v <[p0]>,...,<[p20]>, [pk];
i <[i0]>,...,<[i20]>;
i [mu], [nu], [ro], [si], [LA];
i [om], [omA], [omB];
i [i], [j], [k], [l], [I], [J], [K];
i [a], [b], [c], [d];
cf [f];
t [t];

* variables internal to FORM
s TAG, CUTRAT, SCALE, JJ, KK;
cf TMP, MOM, ABB, ORD, NOW, CH, SIGN(antisymm);
cf NEQ, NN, FF, DROP, JGRAM, D1, D2, E1, E2;
cf JOIN0, JOIN1(symm);
t NUM, EQ, DD, EPS(antisymm);
nt GA, GB, GC, GD;
f CC, WC;
auto s ARG;
set MOMS: k1,...,k`Legs';
set COLS: Col1,...,Col`Legs';

ntable BASIS(1:5, [i]?, [mu]?, [nu]?);
ntable DUAL(1:5, [mu]?, [nu]?);
ntable CHI(6:7);

*#define sig(i,mu,nu) "i_/2*(g_(`~i',`~mu',`~nu') - g_(`~i',`~nu',`~mu'))"
#define sig(i,mu,nu) "i_*(g_(`~i',`~mu',`~nu') - g_(`~i')*d_(`~mu',`~nu'))"

fill BASIS(1) = g6_([i])/2;
fill BASIS(2) = g7_([i])/2;
fill BASIS(3) = g_([i], 6_, [mu])/2;
fill BASIS(4) = g_([i], 7_, [mu])/2;
fill BASIS(5) = 1/2*(g6_([i]) + g7_([i]))/2*`sig([i], [mu], [nu])';

* DUAL includes the 1/2 to cancel the 2 = Tr BASIS(i) DUAL(i)
fill DUAL(1) = 1/2*g6_(1)/2;
fill DUAL(2) = 1/2*g7_(1)/2;
fill DUAL(3) = 1/2*g_(1, 7_, [mu])/2;
fill DUAL(4) = 1/2*g_(1, 6_, [mu])/2;
fill DUAL(5) = 1/2*1/2*`sig(1, [mu], [nu])';

fill CHI(6) = g6_(1)/2;
fill CHI(7) = g7_(1)/2;

.global

*----------------------------------------------------------------------

collect mulM;

moduleoption polyfun=mulM;
.sort

#call Square
argument;
#call Square
endargument;

id mulM(0) = 0;

normalize mulM;

id mulM(1) = 1;

id Eps([mu]?, [nu]?, [ro]?, [si]?) = e_([mu], [nu], [ro], [si]);


if( count(dirM, 1) ) ;

repeat;
  repeat id dirM([x]?, [i]?, [j]?) * dirM([y]?, [j]?, [k]?) =
    dirM([x]*[y], [i], [k]);
  id dirM([x]?, [i]?, [j]?) * dirM([y]?, [k]?, [j]?) =
    dirM([x]*CC([y]), [i], [k]);
  id dirM([x]?, [j]?, [i]?) * dirM([y]?, [j]?, [k]?) =
    dirM(CC([x])*[y], [i], [k]);
  argument dirM, 1;
    argument CC;
      id g_(sM, [mu]?) = -CC([mu]);
      chainin CC;
      id CC(?a) = CC(reverse_(?a));
      id g5_(sM) * CC(?a) = g_(sM, ?a, 5_);
      id g6_(sM) * CC(?a) = g_(sM, ?a, 6_);
      id g7_(sM) * CC(?a) = g_(sM, ?a, 7_);
      id CC(?a) = g_(sM, ?a);
    endargument;
    id CC([x]?) = [x];
  endargument;
endrepeat;

id dirM([x]?, [i]?) * dirM([y]?, [i]?, [j]?) *
     dirM(Spinor(?k, [s2]?)*gi_(sM), [j]?) =
  dirM([x]*[y]*Spinor(?k, -[s2]), [i], [j]);

$fline = 0;
id dirM([x]?, [i]?, [i]?) = -CH([x]);
while( count(CH, 1) );
  $fline = $fline + 1;
  once CH([x]?) = TMP([x]*replace_(sM, $fline));
endwhile;

$fline = 9;
while( count(dirM, 1) );
  $fline = $fline + 1;
  once dirM([x]?, ?i) = TMP([x]*replace_(sM, $fline)) * ORD(?i);
endwhile;

id TMP([x]?) = [x];

if( count(ORD, 1) );
  redefine HaveFermions "1";
  chainin ORD;
  id ORD(?a) = SIGN(?a)*sign_(nargs_(?a)/2);
  mul replace_(SIGN, ORD);
  id ORD(?a) = 1;
endif;

endif;

id g6M([i]?{<10}) = gi_([i]) + g5M([i]);
id g7M([i]?{<10}) = gi_([i]) - g5M([i]);
repeat;
  once g5M([i]?{<10}) = e_(N100_?, N101_?, N102_?, N103_?) *
    g_([i], N100_?, N101_?, N102_?, N103_?)/24;
  renumber;
endrepeat;

mul replace_(g5M, g5_, g6M, g6_, g7M, g7_);

trace4, 1;
trace4, 2;

.sort

#call Contract

*----------------------------------------------------------------------

b g_, `Tensors';
.sort

#if `HaveFermions' == 1

keep brackets;

id gi_([i]?) = g6_([i])/2 + g7_([i])/2;
id g5_([i]?) = g6_([i])/2 - g7_([i])/2;
id g_([i]?, [mu]?) = TMP([i]) * g_([i], [mu]);
repeat id TMP([i]?) * TMP([i]?) = TMP([i]);
id g_([i]?, [mu]?) * TMP([i]?) = g_([i], 6_, [mu])/2 + g_([i], 7_, [mu])/2;

id g6_([i]?) = 2*GA(6)
#if `Gamma5Test' == 1
  + 2*(GA(6) - GA(7)) * Gamma5Test * Dminus4
#endif
  ;
id g7_([i]?) = 2*GA(7)
#if `Gamma5Test' == 1
  - 2*(GA(6) - GA(7)) * Gamma5Test * Dminus4
#endif
  ;
repeat id GA(?g) * g_([i]?, [mu]?) = GA(?g, [mu]);

#call ChainSimplify(1)

#endif

*----------------------------------------------------------------------

#if `SimplifyQ2' == 1

#if `OPP' < 100
* join denominators, as in:
* N2/(D0 D1) + N3/(D0 D1 D2) -> (N2 D2 + N3)/(D0 D1 D2)

b intM;
.sort
collect TMP;

id intM(<[s1]?>,...,<[s`OPP']?>) =
  TAG * JOIN0(<[s1]>,...,<[s`OPP']>);

#do n = {`OPP'+1}, `NPointMax'
id intM(<[s1]?>,...,<[s`n']?>) = TAG * JOIN1(<[s1]>,...,<[s`n']>);
b TAG;
.sort

collect ORD;
splitarg ORD;

repeat id ORD(?a, TMP([x]?) * JOIN0(<[s2]?>,...,<[s`n']?>), ?b,
    TMP([y]?) * JOIN1(Den([p1]?, [m1]?), <[s2]?>,...,<[s`n']?>), ?c) =
  ORD(?a, ?b, TMP([x]*(MOM([p1]) - [m1]) + [y])*
    JOIN1(Den([p1], [m1]), <[s2]>,...,<[s`n']>), ?c);

repeat id ORD([x]?, ?r) = [x] + ORD(?r);
id ORD() = 0;

id TAG * JOIN0(?a) = JOIN0(?a);
mul replace_(JOIN1, JOIN0);

.sort
#enddo

id TAG = 1;
argument TMP;
#call MomSquare
endargument;
id TMP([x]?) = [x];

#endif

b q1, intM;
.sort
keep brackets;

* cancel q^2's in the numerator

repeat;
  once q1.q1 * intM(?a, Den(q1, [m1]?), ?b) =
    TAG * intM(?a, ?b) + [m1] * intM(?a, Den(q1, [m1]), ?b);
  once TAG * q1.q1 * intM(?a, Den([p1]?!{q1}, 0), ?b) =
    replace_(q1, 2*q1 - [p1]) * q1.q1 * intM(?a, Den([p1], 0), ?b);
  also once TAG * q1.q1 * intM(Den([p1]?!{q1}, [m1]?), ?a) =
    replace_(q1, 2*q1 - [p1]) * q1.q1 * intM(Den([p1], [m1]), ?a);
  id TAG = 1;
endrepeat;

id intM() = 0;

#endif

*----------------------------------------------------------------------

.sort

mul replace_(JOIN0, intM);

id intM(Den([p1]?, 0)) = 0;

id Den([p1]?, [m1]?) * [p1]?.[p1]? = 1 + [m1]*Den([p1], [m1]);

*----------------------------------------------------------------------

b q1, intM;
.sort
keep brackets;

argument intM;
id Den([p1]?, [m1]?) = Den([m1])*MOM([p1]);
normalize (0) MOM;
endargument;

#if `SortDen' == 1
symm intM;
#endif

#define MaxPaVe "`OPP'"
#define RatTag ""
#if `OPPRat' == 1
#redefine MaxPaVe "100"
#redefine RatTag "* CUTRAT"
#endif

once intM(Den([m1]?)*MOM([p1]?)) = ORD(0) *
  replace_(q1, 2*q1 - [p1]) * intM(Den(0,[p1],[m1]));
#do n = 2, `NPointMax'
also once intM(<Den([m1]?)*MOM([p1]?)>,...,<Den([m`n']?)*MOM([p`n']?)>) =
  replace_(q1, 2*q1 - [p1]) * 
  intM(<Den(0,[p1],[m1])>*...*<Den({`n'-1},[p`n'],[m`n'])>) * (
#if `n' < `MaxPaVe'
    NN(`n') *
    ORD(<paveM(1)*([p2]-[p1])>+...+<paveM({`n'-1})*([p`n']-[p1])>)
#endif
#if `n' >= `OPP'
    `RatTag'
  + cutM(<[p2]-[p1]>,...,<[p`n']-[p1]>, <[m1]>,...,<[m`n']>)
#endif
  );
#enddo

*----------------------------------------------------------------------

b q1, NUM, ORD, NN, Den, intM, cutM, `Tensors', D, Dminus4, CUTRAT;
.sort
keep brackets;

if( count(ORD, 1) );
  totensor q1, NUM;
elseif( count(cutM, 1) );
  id q1.q1 = q1.q1 - qfM(MuTildeSq);
#if `OPPQSlash' == 1
  id GA([om]?{6,7}[[n]], ?g, q1, ?r) = GA([om], ?g, q1, ?r) +
    i_ * sign_([om] + nargs_(?g)) * GA({7,6}[[n]], ?g, ?r) * MuTilde;
  id MuTilde^2 = qfM(MuTildeSq);
  id MuTilde = 0;
#endif
endif;

#if "`Dim'" == "4"
* add local terms for dimred/CDR as given in Appendix B of
* hep-ph/9806451 (note: 1/(16 pi^2) already included in intM)

#if 0
id NUM([mu]?, [nu]?, [ro]?, [si]?) * NN(4) * intM([x]?) =
  NUM([mu], [nu], [ro], [si]) * NN(4) * intM([x]) -
  5/144 * NEQ([mu], [nu], [ro], [si]) * Finite +
  1/8 * distrib_(1, 2, EQ, NEQ, [mu], [nu], [ro], [si]) * Finite;
#endif

#if 0
also NUM([mu]?, [nu]?, [ro]?) * NN(3) *
       intM(<Den(0,[p0]?,[m0]?)>*...*<Den(2,[p2]?,[m2]?)>) =
  NUM([mu], [nu], [ro]) * NN(3) *
    intM(<Den(0,[p0],[m0])>*...*<Den(2,[p2],[m2])>) +
  1/36 * NEQ([mu], [nu], [ro], [p2] - [p0]) * Finite;
#endif

#if 0
also NUM([mu]?, [nu]?) * NN(3) * intM([x]?) =
  NUM([mu], [nu]) * NN(3) * intM([x]) -
  1/8 * NEQ([mu], [nu]) * Finite;
#endif

id EQ([mu]?, [mu]?) = 1;
id EQ(?a) = 0;

symm NEQ;
id NEQ(?a, [mu]?, [mu]?, ?b) = 0;
id NEQ(?a) = dd_(?a);

id CUTRAT * intM(?a) = 0;
id CUTRAT = 1;
#endif


* decompose into Lorentz-covariant tensors

* The following statement introduces the g_{\mu\nu}'s in a smart way.
* Lifted from: S.A. Larin, T. van Ritbergen, and J.A.M. Vermaseren,
* The optimization of a huge FORM program,
* in: Proceedings Oberammergau 1993, ISBN 9-810-21699-8.

id NUM(?i) = sum_(KK, 0, nargs_(?i), 2,
  paveM(0)^KK * distrib_(1, KK, dd_, NUM, ?i));

id ORD(0) * NUM([mu]?, ?i) = 0;
repeat id ORD([p1]?) * NUM([mu]?, ?i) = ORD([p1]) * d_([p1], [mu]) * NUM(?i);

id ORD(?p) = 1;
id NUM() = 1;

chainin paveM;

*----------------------------------------------------------------------

id D = Dminus4 + 4;

#if "`PaVeReduce'" != "False"

id NN([i]?) * Dminus4 = Dminus4;

#do rep = 1, 1

id NN([n]?) * paveM(?i) * intM([x]?) =
  TMP(NN([n]) * paveM(?i) * intM([x]) * [x]);

b TMP;
.sort
keep brackets;

argument TMP;

* symmetrize the coefficients for N > 4
* hep-ph/0509141 Eq. (6.14+15)
if( match(NN([n]?{>4})) );
  id paveM(0,0,[i1]?,?i) = paveM([i1],0,0,?i);
  id paveM([i]?,[j]?,?i) = TMP([i],[j],TAG,?i) + paveM([j],[i],?i);
  repeat id TMP([i]?,?i,TAG,[j]?,?j) =
    TMP([i],?i,[j],TAG,?j) + paveM([j],?i,[i],?j);
  id TMP(?i,TAG) = paveM(?i);
endif;

* hep-ph/0509141 Eq. (7.13)
id NN(6) * paveM([i1]?,?i) = NN(6) *
  deltap_([i1], 0) * sum_(KK, 1, 5, IGram([i1],KK) * DROP(KK, ?i));

* hep-ph/0509141 Eq. (6.13)
also NN(5) * paveM(0,0,?i) = NN(5) * (
  Dminus4 * paveM(0,0,?i) +
  sum_(KK, 1, 4, JGRAM(KK,0) * DROP(KK, 0,0,?i)) );

* hep-ph/0509141 Eq. (6.12)
also NN(5) * paveM([i1]?,?i) = NN(5) * (
  Dminus4 * paveM([i1],?i) +
  sum_(KK, 0, 4, JGRAM([i1],KK) * DROP(KK, ?i)) -
  2*sum_(KK, 1, 4, E1([i1], KK) * distrib_(1, 1, E2, paveM, ?i)) );

* hep-ph/0509141 Eq. (5.10)
also NN([n]?) * paveM(0,0,?i) = NN([n])/(3 + nargs_(0,0,?i) - [n]) * (
  -Dminus4 * paveM(0,0,?i) +
  FF(0) * paveM(?i) +
  sum_(KK, 1, [n] - 1, FF(KK) * paveM(KK,?i))/2 -
  DROP(0, ?i)/2 );

* hep-ph/0509141 Eq. (5.11+8)
also NN([n]?) * paveM([i1]?,?i) = NN([n]) *
  sum_(KK, 1, [n] - 1, IGram([i1],KK) * (
    DROP(KK, ?i) -
    FF(KK) * paveM(?i) -
    2*D1(KK) * distrib_(1, 1, D2, paveM, ?i) ));

id D1([k]?) * D2([i2]?) * paveM(?i) =
  delta_([k], [i2]) * paveM(0,0,?i);
id E1([i1]?, [k]?) * E2([i2]?) * paveM(?i) =
  JGRAM([i1],[k], 0,[i2]) * DROP([k], 0,0,?i);

* hep-ph/0509141 Eq. (2.28)
id NN([n]?) * JGRAM([s1]?,0) = -NN([n]) *
  sum_(KK, 1, [n], IGram([s1],KK) * FF(KK));
id NN([n]?) * JGRAM([s1]?,[s2]?) = NN([n]) * (
  2*FF(0) * IGram([s1],[s2]) +
  sum_(JJ, 1, [n], sum_(KK, 1, [n],
    IGram([s1],KK, [s2],JJ) * FF(KK) * FF(JJ))) );
* hep-ph/0509141 Eq. (2.29)
id NN([n]?) * JGRAM([s1]?,[s2]?, 0,[s4]?) =
  -NN([n]) * sum_(KK, 1, [n], IGram([s1],[s2], KK,[s4]));

id FF(0) * Den(0, [p0]?, [m0]?) = [m0] * Den(0, [p0], [m0]);
id FF([k]?) * Den(0, [p0]?, [m0]?) * Den([k]?, [pk]?, [mk]?) =
  (MOM([pk] - [p0]) - [mk] + [m0]) *
  Den(0, [p0], [m0]) * Den([k], [pk], [mk]);

id IGram([s1]?,[s2]?) * intM([x]?) =
  sign_([s1] + [s2]) *
  IGram(1, DROP([s1]) * [x]) *
  IGram(1, DROP([s2]) * [x]) *
  IGram(2, [x]);
also IGram([s1]?,[s2]?, [s3]?,[s4]?) * intM([x]?) =
  sign_([s1] + [s2] + [s3] + [s4]) *
  sig_([s1] - [s3]) * IGram(1, DROP([s1]) * DROP([s3]) * [x]) *
  sig_([s4] - [s2]) * IGram(1, DROP([s2]) * DROP([s4]) * [x]) *
  IGram(2, [x]);
id intM(?x) = 1;

argument IGram;
id DROP([k]?) * Den([k]?, ?p) = 1;
endargument;

id IGram(?i, Den(?p)) = IGram(?i);
#do n = 2, `NPointMax'
id IGram(?i, <Den([i1]?,[p1]?,[m1]?)>*...*<Den([i`n']?,[p`n']?,[m`n']?)>) =
  IGram(?i, <[p2]-[p1]>,...,<[p`n']-[p1]>);
#enddo
id IGram(1, ?n1) * IGram(1, ?n2) * IGram(2, ?d) =
  IGram(MOM(?n1) * MOM(?n2), MOM(?d)^2)/2;

#call MomSquare

argument IGram;
id MOM() = 1;
#do n = 2, `NPointMax'
id MOM(<[p2]?>,...,<[p`n']?>) = e_(<[p2]>,...,<[p`n']>);
#enddo
contract;
endargument;

id IGram([x]?, [p1]?.[p1]?) * [p1]?.[p1]? = [x];

#if "`PaVeReduce'" == "True"
argument IGram;
#call kikj
#call InvSimplify(FormSimplify)
endargument;
#endif

id IGram([x]?, [y]?) = [x] * IGram([y]);
factarg IGram;
chainout IGram;
id IGram(0) = IGram(0);
also IGram([x]?number_) = 1/[x];
also IGram([x]?symbol_) = 1/[x];

id NN([n]?) * DROP([k]?, ?i) = NN([n] - 1) * paveM() *
  (deltap_([k], 0) * DROP([k], ?i) - DROP(0, ?i));

* hep-ph/0509141 Eq. (2.8)
repeat;
  id paveM(?n) * DROP(0, 1,?i) * NN([n]?) =
    (-paveM(?n) - sum_(KK, 1, [n] - 1, paveM(?n, KK))) *
    DROP(0, ?i) * NN([n]);
  also paveM(?n) * DROP(0, [h]?,?i) =
    paveM(?n, [h] - theta_([h] - 1)) * DROP(0, ?i);
endrepeat;

repeat id paveM(?n) * DROP([n]?, [h]?,?i) =
  deltap_([n], [h]) * paveM(?n, [h] - theta_([h] - [n])) *
  DROP([n], ?i);

id DROP([k]?) * Den([k]?, ?q) = 1;

#do n = 1, `NPointMax'
id NN(`n') * <Den([i1]?,?p1)>*...*<Den([i`n']?,?p`n')> =
  NN(`n') * intM(<Den(0,?p1)>*...*<Den({`n'-1},?p`n')>);
#enddo

id NN(1) = 1;
id NN([i]?) * Dminus4 = Dminus4;
id paveM() = 1;

endargument;

id TMP([x]?) = [x];

if( count(paveM, 1, NN, 1) == 2 ) redefine rep "0";

.sort

#enddo

#endif

*----------------------------------------------------------------------

b paveM, intM, cutM, NN, Den;
.sort
keep brackets;

id NN(?i) = 1;

id paveM(?i) * intM(Den(0, [p1]?, [m1]?)) =
  theta_(sign_(nargs_(?i))) *
    ([m1]/2)^(nargs_(?i)/2)/fac_(nargs_(?i)/2 + 1) *
    (paveM(1, List(0), [m1]) + [m1]*sum_(KK, 2, nargs_(?i)/2 + 1, 1/KK));
also intM(Den(0, [p1]?, [m1]?)) = paveM(1, List(0), [m1]);
#do n = 2, `NPointMax'
also intM(<Den(0,[p1]?,[m1]?)>*...*<Den({`n'-1},[p`n']?,[m`n']?)>) =
  intM(`n',
#do i = 1, {`n'/2}
    <MOM([p{`i'+1}]-[p1])>,...,<MOM([p`n']-[p{`n'-`i'}])>,
#if {2*`i'} < `n'
    <MOM([p1]-[p{`n'-`i'+1}])>,...,<MOM([p`i']-[p`n'])>,
#endif
#enddo
    <[m1]>,...,<[m`n']>);
#enddo

id Den([p1]?, ?m) = Den(MOM([p1]), ?m);

argument intM, Den;
#call MomSquare
endargument;

id intM([n]?, ?r) * cutM(?k) = cutM([n], List(?r), ?k);
also intM([n]?, ?r) * paveM(?i) = paveM([n], List(?i), ?r);
symm intM:4 3, 4;
id intM([n]?, ?r) = paveM([n], List(0), ?r);

*----------------------------------------------------------------------

#call DotSimplify(#call Contract)

*----------------------------------------------------------------------

#if `HaveFermions' == 1
* Dirac algebra on open fermion chains again

#if "`FermionChains'" == "Weyl"

.sort

* Chisholm's identity backwards to get rid of all Eps
*repeat id GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
*  1/4 * sign_([om]) * sign_(nargs_(?a)) * (
*    GA([om], ?a, [mu], [nu], [ro], ?b) -
*    GA([om], ?a, [ro], [nu], [mu], ?b) );
repeat;
  once GA([om]?, ?a, [LA]?, ?b) * e_([mu]?, [nu]?, [ro]?, [LA]?) =
    sign_([om] + nargs_(?a)) * (
      GA([om], ?a, [mu], [nu], [ro], ?b) -
      d_([mu], [nu]) * GA([om], ?a, [ro], ?b) +
      d_([mu], [ro]) * GA([om], ?a, [nu], ?b) -
      d_([nu], [ro]) * GA([om], ?a, [mu], ?b) );
endrepeat;

*#elseif "`Dim'" == "4"
#elseif 0

.sort

* this is Chisholm's identity:
repeat;
  once GA([om]?, [mu]?, [nu]?, [ro]?, ?a) =
    sign_([om]) * GA([om], N100_?, ?a) * e_([mu], [nu], [ro], N100_?) +
    d_([mu], [nu]) * GA([om], [ro], ?a) -
    d_([mu], [ro]) * GA([om], [nu], ?a) +
    d_([nu], [ro]) * GA([om], [mu], ?a);
  renumber;
endrepeat;

#call Contract

#endif

b `Tensors';
.sort

keep brackets;

#call ChainSimplify(1)

#endif

*----------------------------------------------------------------------

#if "`Dim'" != 4
b D, Dminus4, paveM, cutM, CUTRAT;
.sort

keep brackets;

id D = Dminus4 + 4;

#if "`Dim'" == "D"

#if `OPPRat' == 0
id Dminus4 * cutM(?a) = dm4M * cutM(?a);
#else
id Dminus4 * cutM(?a) = 0;
#endif

* add local terms for dimreg
also Dminus4 * paveM(1, List(0), [m1]?) = -2*[m1]*Finite;
also Dminus4 * paveM(1, List(0,0), [m1]?) = -[m1]^2/2*Finite;
also Dminus4 * paveM(2, List(0), ?a) = -2*Finite;
also Dminus4 * paveM(2, List(1), ?a) = Finite;
also Dminus4 * paveM(2, List(0,0), [k1]?, [m1]?, [m2]?) =
  1/6*([k1] - 3*[m1] - 3*[m2])*Finite;
also Dminus4 * paveM(2, List(1,1), ?a) = -2/3*Finite;
also Dminus4 * paveM(2, List(0,0,1), [k1]?, [m1]?, [m2]?) =
  -1/12*([k1] - 2*[m1] - 4*[m2])*Finite;
also Dminus4 * paveM(2, List(1,1,1), ?a) = 1/2*Finite;
also Dminus4 * paveM(3, List(0,0), ?a) = -1/2*Finite;
also Dminus4 * paveM(3, List(0,0,[i]?), ?a) = 1/6*Finite;
also Dminus4 * paveM(3, List(0,0,0,0), [k1]?, [k2]?, [k1k2]?, [m1]?, [m2]?, [m3]?) =
  1/48*([k1] + [k2] + [k1k2] - 4*([m1] + [m2] + [m3]))*Finite;
also Dminus4 * paveM(3, List(0,0,[i]?,[i]?), ?a) = -1/12*Finite;
also Dminus4 * paveM(3, List(0,0,[i]?,[j]?), ?a) = -1/24*Finite;
also Dminus4 * paveM(4, List(0,0,0,0), ?a) = -1/12*Finite;
also Dminus4 * paveM(4, List(0,0,0,0,[i]?), ?a) = 1/48*Finite;

also Dminus4 = 0;

#redefine Dim "4"

#endif

id CUTRAT * paveM(?a) = 0;
id CUTRAT = 1;

#endif

*----------------------------------------------------------------------

#if `HaveFermions' == 0

#call Abbreviate

#else

#if "`FermionChains'" == "Weyl"

#if `OPP' < 100
id Spinor([p1]?MOMS[[n]], 0, [s1]?) * GA([om]?, ?g) =
  HelDelta(mulM(HEL([n]) + [s1]*sign_([om]))) *
  Spinor([p1], 0, [s1]) * GA([om], ?g);
id GA([om]?, ?g) * Spinor([p1]?MOMS[[n]], 0, [s1]?) =
  HelDelta(mulM(HEL([n]) - [s1]*sign_([om] + nargs_(?g)))) *
  GA([om], ?g) * Spinor([p1], 0, [s1]);

id HelDelta(mulM([x]?number_)) = delta_([x]);

moduleoption polyfun=HelDelta;
.sort

normalize HelDelta;
#endif

b Spinor;
.sort

keep brackets;

id Spinor(?a) * GA(?g) * Spinor(?b) =
  CH(Spinor(?a, 2, 0), ?g, Spinor(?b, 1, 0));

repeat;
  once CH([s1]?, [x]?, ?a, [LA]?, ?b) *
       CH([s2]?, [y]?, ?c, [LA]?, ?d) =
    WC(sign_(nargs_(?a, ?c) + [x] + [y]), [s1], [x], ?a) *
    WC(?b) * WC([s2], [y], ?c) * WC(?d);

* Fierz 1: <A|sig_mu|B> <C|sigbar^mu|D> = 2 <A|D> <C|B>
  id WC(-1, ?a) * WC(?b) * WC(?c) * WC(?d) =
    2 * CH(?a, ?d) * CH(?c, ?b);

* Fierz 2: <A|sig(bar)_mu|B> <C|sig(bar)^mu|D> = 2 <A|eps|C> <B|eps|D>
  also WC(1, ?a) * WC(?b, Spinor(?s1, [s1]?)) *
         WC(Spinor(?s2, [s2]?), [x]?, ?c) * WC(?d) =
    2 * CH(?a, reverse_(?c), Spinor(?s2, 1 - [s2])) *
      CH(Spinor(?s1, 1 - [s1]),
        7 - mod_([x] + nargs_(?b, ?c), 2),
        reverse_(?b), ?d);

* due to the canonical ordering of the Dirac chains this
* is the only(?) case we need of Fierz on the same chain:
  repeat id CH(?a, [LA]?, [LA]?, ?b) = 4*CH(?a, ?b);
endrepeat;

id CH([s1]?, ?g, [s2]?) = [s1] * GA(?g) * [s2];

#call ChainSimplify(0)
*id CH(?a, [p1]?, [p1]?, ?b) = [p1].[p1] * CH(?a, ?b);

#call kikj
#call eiki

id Spinor(?s1) * GA(?g) * Spinor(?s2) =
  abbM(fmeM(WeylChain(Spinor(?s1), ?g, Spinor(?s2))))
#if "`Scale'" != "1"
  * MOM(?g)
#endif
  ;

#call Abbreviate

#else

#switch "`FermionOrder'"

#case "None"
.sort
#call DiracFinal
#break

#case "Fierz"
* Fierz twice for simplification
#call FierzBefore
#call FierzUnordered
#call FierzUnordered
#call DiracFinal
#break

#case "Automatic"
* lexicographical ordering
#call FierzBefore
#call FierzOrdered
#call DiracFinal
#break

#case "Colour"
* postponed until after SUNT simplification
#break

#default
#do i = {`FermionOrder'}
#ifndef `order'
#define order "k{`i'}"
#else
#redefine order "`order',k{`i'}"
#endif
#enddo
mul ORD(`order');
#call FierzBefore
#call FierzOrdered
#call DiracFinal

#endswitch

#endif

#endif

*----------------------------------------------------------------------

#ifdef `Inserted'
.sort
#else
#call DoInsertions
#define Inserted
#endif

id paveM(1, [x]?, 0) = 0;

id Den(0, [x]?) = -Den([x], 0);
id Den([x]?, 0) * [x]? = 1;

#call Const
.sort

collect mulM;

argument mulM;
#call Neglect
endargument;

argument;
#call Square
endargument;

id mulM(0) = 0;

*----------------------------------------------------------------------
* index handling

.sort

repeat;
  once SumOver([i]?, ?a, Renumber) =
    TMP(N100_?) * SumOver(N100_?, ?a) * replace_([i], N100_?);
  renumber;
endrepeat;

.sort

renumber 1;

id IndexEps([i]?, [j]?, [k]?) = EPS([i], [j], [k]);

repeat;
  id EPS([I]?, [J]?, [K]?) * EPS([I]?, [J]?, [K]?) *
    SumOver([I]?, 3) * SumOver([J]?, 3) * SumOver([K]?, 3) = 6;
  id EPS([I]?, [J]?, [k]?) * EPS([I]?, [J]?, [c]?) *
    SumOver([I]?, 3) * SumOver([J]?, 3) = 2*IndexDelta([k], [c]);
  id EPS([I]?, [j]?, [k]?) * EPS([I]?, [b]?, [c]?) *
    SumOver([I]?, 3) =
    IndexDelta([j], [b])*IndexDelta([k], [c]) -
    IndexDelta([j], [c])*IndexDelta([k], [b]);
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
      SUNF(?a, [a], [b], N100_?) * SUNF(N100_?, [c], [d]) * SUNSum(N100_?);
    renumber;
  endrepeat;

* f^{abc} = 2 i Tr(T^c T^b T^a - T^a T^b T^c)

  id SUNF([a]?, [b]?, [c]?) =
    2*i_*(SUNT([c], [b], [a], 0, 0) - SUNT([a], [b], [c], 0, 0));

endif;


repeat;
  once SUNT(?a, 0, 0) = SUNT(?a, N100_?, N100_?) * SUNSum(N100_?);
  renumber;
endrepeat;

repeat;
  once SUNT(?a, [a]?, [b]?, [i]?, [j]?) =
    SUNT(?a, [a], [i], N100_?) * SUNT([b], N100_?, [j]) * SUNSum(N100_?);
  renumber;
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
    SUNT([j], [c])*SUNT([k], [b]);
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

id EPS([i]?, [j]?, [k]?) = sunM(SUNEps([j], [k], [i]));

* cleaning up, step 2: bead up the SUNTs into chains

repeat;
  once SUNT(?a, [a]?, [i]?, [j]?) = TMP(?a, [a], [i], [j]);
  repeat;
    id TMP(?a, [i]?, [j]?) * SUNT(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
    id SUNT(?a, [i]?, [j]?) * TMP(?b, [j]?, [k]?) * SUNSum([j]?, ?c) =
      TMP(?a, ?b, [i], [k]);
  endrepeat;

  id TMP(?a, [i]?, [i]?) * SUNSum([i]?, ?b) = SUNTr(?a);
  id SUNTr(?a) = TMP(?a, 0, 0);

* special case of Tr(T^a T^b) = 1/2 delta_{ab}
*  id TMP([a]?, [a]?, 0, 0) = 1/2;
  id TMP([x]?int_, [y]?int_, 0, 0) = 1/2*delta_([x], [y]);

  id TMP(?a) = sunM(SUNT(?a));
endrepeat;

symm SUNT;

id SUNT(?a) = sunM(SUNT(?a));

#if "`FermionOrder'" == "Colour"
id sunM(SUNT(?a, [i]?COLS[[x]], [j]?COLS[[y]])) =
  sunM(SUNT(?a, [i], [j])) * ORD(MOMS[[x]], MOMS[[y]]);
#endif

repeat id sunM([x]?) * sunM([y]?) = sunM([x] * [y]);

id SUNSum([i]?, [x]?) = [x];

* the Mat(...) are kept at the almost outermost level (only SumOver
* comes before), i.e. the amplitude is of the form Sum[c[i] Mat[i], i];
* this is required for the calculation of the squared amplitude

id sunM([x]?) = Mat(sunM([x]));

#if "`FermionOrder'" == "Colour"
#call FierzBefore
#call FierzOrdered
#call DiracFinal
#endif

#endif

*----------------------------------------------------------------------

.sort

mul replace_(i_, I);
id I = mulM(I);

repeat id Mat([x]?) * Mat([y]?) = Mat([x] * [y]);

*----------------------------------------------------------------------

moduleoption polyfun=mulM;
.sort

normalize mulM;

moduleoption polyfun=abbM;
.sort

normalize abbM;
id abbM(1) = 1;

*----------------------------------------------------------------------

b SumOver, PowerOf, Mat, Den, IGram,
  abbM, mulM, paveM, cutM, qfM, Dminus4, dm4M;
.sort

collect FormSimplify;

moduleoption polyfun=FormSimplify;
.sort

normalize FormSimplify;
#call Factor(FormSimplify)
#call InvSimplify(FormSimplify)

id FormSimplify(?a) = TMP(FormSimplify(?a));
id abbM(?a) = TMP(abbM(?a));

moduleoption polyfun=TMP;
.sort

*----------------------------------------------------------------------

normalize TMP;
id TMP(1) = 1;
id TMP(?a) = dum_(?a);

moduleoption polyfun=mulM;
.sort

argument mulM;
#call Square
endargument;

normalize mulM;
#call Factor(mulM)

*----------------------------------------------------------------------

#if `OPP' < 100

b cutM, paveM, qfM, Dminus4, dm4M, SumOver, PowerOf, Mat, Den;
.sort

collect FormSimplify;
normalize FormSimplify;

b cutM, paveM, qfM, Dminus4, dm4M, SumOver, PowerOf, Mat;
.sort

collect qcM;
normalize qcM;

id SumOver([i]?, ?a) = SumOver([i], ?a) * ORD([i]);
chainin ORD;
id qcM([x]?) * ORD(?i) = qcM([x], List(?i)) * ORD(?i);

id qfM([x]?) = [x];

b cutM, paveM, SumOver, PowerOf, Mat, ORD;
.sort

collect numM, numM;
normalize numM;

id numM(qcM(?a)) = qcM(?a) * numM(1);
also numM([x]?) * ORD(?i) = numM([x], List(?i));

id ORD(?i) = 1;

id cutM([n]?, ?a) * numM(?x) = cutM([n], numM(?x), ?a);

id numM([x]?, ?i) = [x];

#endif

*----------------------------------------------------------------------

b SumOver, PowerOf, Mat, Den, IGram, paveM, cutM;
print;

.end

