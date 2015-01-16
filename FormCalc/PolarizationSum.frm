* PolarizationSum.frm
* the FORM part of the PolarizationSum function
* this file is part of FormCalc
* last modified 17 Apr 13 th


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
id `foo'(?x) = mulM(`foo'(?x));
argument mulM;
factarg `foo';
chainout `foo';
makeinteger `foo';
id `foo'([x]?) = `foo'(nterms_([x]), [x]);
id `foo'(1, [x]?) = [x];
id `foo'([n]?, [x]?) = `foo'([x]);
endargument;
makeinteger mulM;
#endprocedure

***********************************************************************

#if "`MomElim'" == "Automatic"
#define MomRange "1, `Legs'"
#elseif `MomElim'
#define MomRange "`MomElim', `MomElim'"
#endif

#procedure DotSimplify
#call eiki

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

#do rep1 = 1, 1
b dotM;
.sort
keep brackets;

#ifdef `MomRange'
id dotM([x]?) = dotM(nterms_([x]), [x]);

#do rep2 = 1, 2
#do i = `MomRange'
#ifdef `k`i''
id dotM([n]?, [x]?) = dotM([n], [x]) * NOW([x]);
argument NOW;
id k`i' = `k`i'';
#call eiki
endargument;

id NOW(0) = 0;
id NOW([x]?) = dotM(nterms_([x]), [x]);
once dotM(?a) = dotM(?a);
also dotM(?a) = 1;
#endif
#enddo
#enddo

id dotM([n]?, [x]?) = dotM([x]);
#endif

argument dotM;
#call kikj
#call Square
endargument;
#call InvSimplify(dotM)
id dotM(0) = 0;
#enddo

#if `DotExpand' == 1

id dotM([x]?) = [x];

.sort
off oldFactArg;

id TAG = 1;

#else

factarg dotM;
chainout dotM;
makeinteger dotM;
id dotM(1) = 1;

ab `Vectors', `Invariants', dotM;
.sort
off oldFactArg;

collect dotM;

#call InvSimplify(dotM)
id dotM(0) = 0;

repeat id TAG * dotM([x]?) = TAG * [x];
id TAG = 1;

makeinteger dotM;
id dotM(dotM(?x)) = dotM(?x);

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

b dotM;
.sort
keep brackets;

toPolynomial;

.sort

#endif
#endprocedure

***********************************************************************

#procedure Abbreviate
#call DotSimplify

id [p1]?.[p2]? = abbM([p1].[p2], [p1], [p2]);

id e_([mu]?, [nu]?, [ro]?, [si]?) =
  abbM(e_([mu], [nu], [ro], [si]), [mu], [nu], [ro], [si]);

id d_([mu]?, [nu]?) = abbM(d_([mu], [nu]), [mu], [nu]);

id [t]?(?a) = abbM([t](?a), ?a);

id [p1]?([mu]?) = abbM([p1]([mu]), [p1]);

repeat;
  once abbM([x]?, ?a, [mu]?!fixed_, ?b) *
       abbM([y]?, ?c, [mu]?, ?d) =
    abbM([x]*[y], ?a, ?b, ?c, ?d) * replace_([mu], N100_?);
  also once abbM([x]?, ?a, [mu]?!fixed_, ?b, [mu]?, ?c) =
    abbM([x], ?a, ?b, ?c) * replace_([mu], N100_?);
  renumber;
endrepeat;

id abbM([x]?, ?a) = abbM([x]);

#call Square

moduleoption polyfun=abbM;
.sort

makeinteger abbM;
id abbM(1) = 1;
#endprocedure

***********************************************************************

#procedure Prepare
#call Square
#call eiei

#call ConstBracket
.sort

collect mulM;

id Conjugate([f]?LOOPINT(?x)) = Conjugate([f](?x));
also Conjugate(?x) = mulM(Conjugate(?x));

moduleoption polyfun=mulM;
.sort

#call Factor(mulM)

b mulM;
.sort
keep brackets;

argument mulM;
toPolynomial;
endargument;

toPolynomial;

.sort
#endprocedure

***********************************************************************

#procedure PolSum(i, m, d)
b e`i', ec`i', z`i', zc`i', eT`i', eTc`i';
.sort
d `d';
keep brackets;

id e`i' = ET(?);
id ec`i' = ETC(?);
id z`i' = ET(?);
id zc`i' = ETC(?);
mul DF;

#if `m' == "0"
* massless case

id DF * ET([mu]?) * ETC([nu]?) = -d_([mu], [nu]) +
  (d_(eta`i', [mu])*d_(k`i', [nu]) +
   d_(eta`i', [nu])*d_(k`i', [mu]))/(eta`i'.k`i');
* The eta are gauge-dependent vectors.
* Their appearance in the result is supposed to alert
* the user to the presence of gauge-dependent terms.
* The eta must fulfill eta.eta = e.eta = 0 and k.eta != 0.
* Instead of imposing eta.eta = 0 one can add
* - d_(k`i', [mu])*d_(k`i', [nu])*(eta`i'.eta`i')/(eta`i'.k`i')^2

also DF * eT`i'([mu]?, [nu]?) * eTc`i'([ro]?, [si]?) =
  d_([mu], [ro])*d_([nu], [si]) +
  d_([mu], [si])*d_([nu], [ro]) -
  d_([mu], [nu])*d_([ro], [si]);

also DF = `Dim' - 2;

#else
* massive case

id DF * ET([mu]?) * ETC([nu]?) = -d_([mu], [nu]) +
  k`i'([mu])*k`i'([nu])/(`m')^2;

also DF = `Dim' - 1;

#call Square
#endif

.sort
d `Dim';

#call eiei
#call eiki
#call kikj
#call Neglect
#endprocedure

***********************************************************************

#procedure Emit
.sort

id D = Dminus4 + 4;

contract;

* Cycling the momenta at this point (rather than in DotSimplify)
* is necessary to obtain all possible terms coming from the
* cancellation of an etaM.kN before sending the etaM to zero.
#do i = `MomRange'
#ifdef `k`i''
id k`i' = `k`i'';
.sort;
#endif
#enddo

#if `GaugeTerms' == 0
#do i = 1, `Legs'
id eta`i' = 0;
#enddo
#endif

id D = Dminus4Eps + 4;

#call Abbreviate

#if 0
b SumOver, Conjugate, Den, abbM;
.sort

collect dotM;
#call Factor(dotM)
#call InvSimplify(dotM)
#endif

b SumOver, Conjugate, Den;
.sort

collect mulM;
makeinteger mulM;

argument mulM;
toPolynomial;
endargument;

b Conjugate, mulM;
.sort
keep brackets;

toPolynomial;

.sort

#write "%X"

b SumOver, Den;
print +s;
.end
#endprocedure

***********************************************************************

#define LoopInt "A0, A00, B0, B1, B00, B11, B001, B111, A0i, B0i, C0i, D0i, E0i, F0i"

cf SumOver, Den, Conjugate;
cf `LoopInt';
s D, Dminus4, Dminus4Eps, `Invariants';

s DF, TAG, ETAG;
t ET, ETC;
cf TMP, NOW;
auto s ARG;
set LOOPINT: `LoopInt';
set INVS: `Invariants';

i [mu], [nu], [ro], [si];
v [p1], [p2];
s [x], [y], [n];
t [t];
cf [f];

extrasymbols array subM;
cf abbM, dotM, mulM, powM;

