* PolarizationSum.frm
* the FORM part of the PolarizationSum function
* this file is part of FormCalc
* last modified 16 Mar 18 th


#procedure Prepare
b `Vectors', Mat, SumOver;
.sort

collect mulM, mulM;

argument mulM;
argument Conjugate;
toPolynomial onlyfunctions;
endargument;
toPolynomial onlyfunctions;
endargument;

.sort

#call Factor(mulM)

.sort

toPolynomial onlyfunctions mulM;

.sort
drop;
#endprocedure

***********************************************************************

#procedure PolSum(i, m, d)
b e`i', ec`i', z`i', zc`i', Pol;
.sort
d `d';
keep brackets;

id e`i' = ET(?);
id ec`i' = ETC(?);
id z`i' = ET(?);
id zc`i' = ETC(?);
mul TAG;

#if `m' == "0"
* massless case

id TAG * ET([mu]?) * ETC([nu]?) = -d_([mu], [nu])
#if "`GaugeTerms'" != "Off"
  + (d_(eta`i', [mu])*d_(k`i', [nu]) +
     d_(eta`i', [nu])*d_(k`i', [mu]))/(eta`i'.k`i')
#endif
  ;
* The eta are gauge-dependent vectors.
* Their appearance in the result is supposed to alert
* the user to the presence of gauge-dependent terms.
* The eta must fulfill eta.eta = e.eta = 0 and k.eta != 0.
* Instead of imposing eta.eta = 0 one can add
* - d_(k`i', [mu])*d_(k`i', [nu])*(eta`i'.eta`i')/(eta`i'.k`i')^2

also TAG * Pol(e`i', [mu]?, [nu]?) * Pol(ec`i', [ro]?, [si]?) =
  d_([mu], [ro])*d_([nu], [si]) +
  d_([mu], [si])*d_([nu], [ro]) -
  d_([mu], [nu])*d_([ro], [si]);

also TAG = `Dim' - 2;

#else
* massive case

id TAG * ET([mu]?) * ETC([nu]?) = -d_([mu], [nu]) +
  k`i'([mu])*k`i'([nu])/(`m')^2;

also TAG = `Dim' - 1;

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
#do i = {`MomRange'}
#ifdef `k`i''
id k`i' = `k`i'';
.sort;
#endif
#enddo

#if "`GaugeTerms'" == "False"
#do i = 1, `Legs'
id eta`i' = 0;
#enddo
#endif

#call EtaSubst

id D = Dminus4Eps + 4;

#call DotSimplify
#call Abbreviate

.sort

#write "%X"

b `Vectors', Mat, SumOver;
print +s;
.end
#endprocedure

***********************************************************************

#call CommonDecl

s D, Dminus4, Dminus4Eps;
cf SumOver, Den, Mat, Conjugate;
cf [f];
t ET, ETC;

#define LoopInt "A0, A00, B0, B1, B00, B11, B001, B111, A0i, B0i, C0i, D0i, E0i, F0i"
cf `LoopInt';
set LOOPINT: `LoopInt';

