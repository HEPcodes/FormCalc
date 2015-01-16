* PolarizationSum.frm
* the FORM part of the PolarizationSum function
* this file is part of FormCalc
* last modified 17 Nov 07 th


#procedure PolSum(i, m)
#$dim = 4;
if( count(z`i',1, zc`i',1) ) $dim = Dminus4;

b z`i', zc`i', e`i', ec`i', eT`i', eTc`i';
.sort
d `$dim';

keep brackets;

id e`i' = E(?);
id ec`i' = EC(?);
id z`i' = E(?);
id zc`i' = EC(?);

#if `m' == "0"
* massless case

multiply 2;
id E([mu]?) * EC([nu]?) = 1/2*( -d_([mu], [nu])
#if `GaugeTerms' == 1
* The eta are gauge-dependent vectors.
* Their appearance in the result is supposed to alert
* the user to the presence of gauge-dependent terms.
* The eta must fulfill e.eta = 0 and k.eta != 0.
    - d_(k`i', [mu])*d_(k`i', [nu])*(eta`i'.eta`i')/(eta`i'.k`i')^2
    + (d_(eta`i', [mu])*d_(k`i', [nu]) +
       d_(eta`i', [nu])*d_(k`i', [mu]))/(eta`i'.k`i')
#endif
  );

id eT`i'([mu]?, [nu]?) * eTc`i'([rho]?, [sig]?) = 1/2*(
  d_([mu], [rho])*d_([nu], [sig]) +
  d_([mu], [sig])*d_([nu], [rho]) -
  d_([mu], [nu])*d_([rho], [sig]) );

#else
* massive case

multiply 3;
id E([mu]?) * EC([nu]?) = 1/3*( -d_([mu], [nu]) +
  k`i'([mu])*k`i'([nu])/(`m')^2 );

#endif

.sort

#call eiei
#call eiki
#call kikj
#call Neglect
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
repeat id `foo'([x]?, [y]?, ?a) = `foo'([x]) * `foo'([y], ?a);
id `foo'([x]?number_) = [x];
id `foo'([x]?symbol_) = [x];
#endprocedure

***********************************************************************

#procedure Emit
contract 0;

#do i = 1, `Legs'
#ifdef `k`i''
id k`i' = `k`i'';
.sort;
#endif
#enddo

#call kikj
#call Neglect

.sort

#if "`Scale'" != "1"
$pow = count_(<k1,1>,...,<k`Legs',1>);
if( $pow != 0 ) multiply pow(`Scale', $pow/2);
#endif

#call Square

id [p1]?.[p2]? = abb([p1].[p2]);
id 1/[p1]?.[p2]? = abb(1/[p1].[p2]);
id e_([mu]?, [nu]?, [rho]?, [sig]?) = abb(e_([mu], [nu], [rho], [sig]));
id [t]?(?a) = abb([t](?a));

b abb, `Bracket';
.sort

collect Simplify, Simplify;
normalize Simplify;

#call Factor(Simplify)
#call InvSimplify(Simplify)

id Simplify(0) = 0;

.sort

moduleoption polyfun=abb;
.sort

normalize abb;
id abb(1) = 1;

b abb, `Bracket';
print;
.end
#endprocedure

***********************************************************************

#define Bracket "Den, A0, A00, B0, B1, B00, B11, B001, B111, A0i, B0i, C0i, D0i, E0i, F0i"

s Dminus4;
i [mu], [nu], [rho], [sig];
v [p1], [p2];
s [x], [y];
t [t];

cf abb, Simplify;
cf `Bracket';
t E, EC;

auto s FC;

.global

