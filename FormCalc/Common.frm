* Common.frm
* FORM procedures common to CalcFeynAmp, HelicityME, PolarizationSum
* this file is part of FormCalc
* last modified 19 Mar 18 th


#procedure CommonDecl
cf MetricTensor, Eps, Mat;
t Pol;

extrasymbols array subM;
cf abbM, fermM, dotM, addM, mulM;

nt GA;
f GF;
cf TMP, ABB;
auto s ARG;
s `Invariants', TAG, ETAG, QTAG;
set INVS: `Invariants';
set LOOPMOM: `LoopMomenta';

i [mu], [nu], [ro], [si], [i];
s <[m0]>,...,<[m20]>;
s <[s0]>,...,<[s20]>;
v <[p0]>,...,<[p20]>, [q1];
s [x], [y], [n];
t [t];
#endprocedure

***********************************************************************

#procedure EiKi(e, k)
id d_(`e', `k') = 0;
id Pol(`e',?a, `k',iM?) = 0;
id Pol(?a,`e', iM?,`k') = 0;
#endprocedure

***********************************************************************

#procedure Fewest(foo)
argument `foo';
#call Neglect
endargument;
id `foo'([x]?, [y]?) = `foo'([x], nterms_([x])*2 - 1, [y], nterms_([y])*2);
symm `foo' (2,1), (4,3);
id `foo'([x]?, ?a) = `foo'([x]);
#endprocedure

***********************************************************************

#procedure Factor(foo)
factarg `foo';
id `foo'(?x) = TMP(`foo'(?x));
argument TMP;
chainout `foo';
makeinteger `foo';
id `foo'([x]?) = `foo'(nterms_([x]), [x]);
id `foo'(1, [x]?) = [x];
id `foo'([n]?, [x]?) = `foo'([x]);
endargument;
makeinteger TMP;
id TMP(?a) = mulM(?a);
#endprocedure

***********************************************************************

#if (("`MomElim'" == "Automatic") && (isdefined(k1)))

#procedure MomConserv(foo)
id `foo'([x]?) = `foo'(TMP([x], [x]));
argument `foo';
argument TMP, 1;
#call kikj
#call Square
endargument;
id TMP([x]?, [y]?) = TMP(nterms_([x]), [x], [y]);

#do rep = 1, 2
#do i = 1, `Legs'
id TMP([n]?, [x]?, [y]?, ?a) = TMP([x], [n], [x], [y]);
argument TMP, 1;
id k`i' = `k`i'';
#call eiki
endargument;
id TMP([x]?, ?a) = TMP([x], [x], ?a);
argument TMP, 1;
#call kikj
#call Square
endargument;
id TMP(0, ?a) = 0;
id TMP([x]?, ?a) = TMP(nterms_([x]), [x], ?a);
symm TMP (1,2,3) (4,5,6);
#enddo
#enddo

id TMP([n]?, [x]?, [y]?, ?a) = [x];
endargument;

#call InvSimplify(`foo')
id `foo'(0) = 0;
#endprocedure

#else

#procedure MomConserv(foo)
argument `foo';
#ifdef `k`MomElim''
id k`MomElim' = `k`MomElim'';
#call eiki
#endif
#call kikj
#call Square
endargument;

id `foo'(0) = 0;
#call InvSimplify(`foo')
id `foo'(0) = 0;
#endprocedure

#endif

***********************************************************************

#procedure DotSimplify
#call eiki

id GA(?g) = GF(?g);

id [q1]?LOOPMOM = [q1] * QTAG;
id e_([mu]?, [nu]?, [ro]?, [si]?) =
  e_([mu], [nu], [ro], [si]) * TMP([mu], [nu], [ro], [si]) * ETAG;
id [t]?(?i) = [t](?i) * TMP(?i);
chainout TMP;
id TMP([p1]?) = 1;
id TMP([mu]?)^2 = 1;
id TMP([mu]?) = TAG;
id ETAG^[n]?{>1} = ETAG;

ab k1,...,k`Legs';
.sort
on oldFactArg;

collect dotM, dotM, 50;
makeinteger dotM;

id ETAG = 1;
id QTAG = TAG;

b dotM;
.sort
keep brackets;

#call MomConserv(dotM)

#if `DotExpand' == 1

id dotM([x]?) = [x];

.sort
off oldFactArg;

id TAG = 1;

#else

b dotM;
.sort
keep brackets;

factarg dotM;
chainout dotM;
id dotM([n]?number_) = [n];
makeinteger dotM;

ab `Vectors', `Invariants', dotM;
.sort
off oldFactArg;

collect dotM, dotM, 50;

#call InvSimplify(dotM)
id dotM(0) = 0;

repeat id TAG * dotM([x]?) = TAG * [x];
id TAG = 1;

*makeinteger dotM;
*id dotM(dotM(?x)) = dotM(?x);

argument dotM;
id dotM([x]?) = dotM(nterms_([x]), [x]);
id dotM(1, [x]?) = [x];
id dotM([n]?, [x]?) = dotM([x]);
argument dotM;
toPolynomial;
endargument;
toPolynomial;
endargument;

makeinteger dotM;
id dotM(1) = 1;
id dotM([x]?^[n]?) = dotM([x])^[n];
id dotM([x]?INVS) = [x];

toPolynomial onlyfunctions dotM;

.sort

#endif

id GF(?g) = GA(?g);
#endprocedure

***********************************************************************

#procedure Abbreviate
id [p1]?.[p2]? = ABB(0, [p1].[p2], [p1], [p2]);

id e_([mu]?, [nu]?, [ro]?, [si]?) =
  ABB(0, Eps([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]);

id d_([mu]?, [nu]?) = ABB(0, MetricTensor([mu], [nu]), [mu], [nu]);

id Pol(?a) = ABB(0, Pol(?a), ?a);

id [p1]?([mu]?) = ABB(0, [p1]([mu]), [p1]);

repeat;
  once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b) *
       ABB([s2]?, [y]?, ?c, [mu]?, ?d) =
    ABB([s1] + [s2], [x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], N100_?);
  renumber;
  once ABB([s1]?, [x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    ABB([s1], [x], ?a, ?b, ?c) * replace_([mu], N100_?);
  renumber;
endrepeat;

id ABB(0, [x]?, ?a) = abbM([x]);
id ABB([i]?, [x]?, ?a) = fermM([x]);

#if "`FermionChains'" != "Weyl"
repeat id fermM([x]?) * fermM([y]?) = fermM([x] * [y]);
argument fermM;
toPolynomial;
endargument;
id fermM([x]?) = Mat(fermM([x]));
#endif

argument abbM, Mat;
toPolynomial;
endargument;

b addM, mulM;
moduleoption polyfun=abbM;
.sort

b abbM;
.sort
on oldFactArg;
keep brackets;

id abbM([x]?) = abbM(nterms_([x]), [x]);
id abbM(1, [x]?) = TMP([x]);
also abbM([n]?, [x]?) = abbM([x]);

factarg abbM;
chainout abbM;
id TMP([x]?) = abbM([x]);

makeinteger abbM;
id abbM(1) = 1;

b abbM;
.sort
off oldFactArg;
keep brackets;

toPolynomial onlyfunctions abbM;
#endprocedure

***********************************************************************

#procedure CollectTerms
collect dotM;

moduleoption polyfun=dotM;
.sort
on oldFactArg;

#call Factor(dotM)

.sort

argument mulM;
toPolynomial;
endargument;

moduleoption polyfun=mulM;
.sort

#call Factor(mulM)

b mulM;
.sort
off oldFactArg;
keep brackets;

argument mulM;
toPolynomial;
endargument;

id mulM([x]?symbol_) = [x];

toPolynomial onlyfunctions mulM;
#endprocedure

