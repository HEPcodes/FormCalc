* PolarizationSum.frm
* the FORM part of the PolarizationSum function
* this file is part of FormCalc
* last modified 16 Oct 05 th


#procedure PolSum(i, m)
b e`i', ec`i';
.sort

keep brackets;

id e`i' = E(?);
id ec`i' = EC(?);

#if `m' == "0"
* massless case

multiply 2;
id E([mu]?) * EC([nu]?) = 1/2*( -d_([mu], [nu])
#if `GaugeTerms' == 1
* Note: we really ought to use new vectors eta here,
* the e are re-used for convenience only.
* Their appearance in the result is supposed to alert
* the user to the presence of gauge-dependent terms.
* The eta must fulfill e.eta = 0 and k.eta != 0.
    - k`i'([mu])*k`i'([nu])*(e`i'.e`i')*inv(k`i'.e`i')^2
    + (e`i'([mu])*k`i'([nu]) + e`i'([nu])*k`i'([mu]))*inv(k`i'.e`i')
#endif
  );

#else
* massive case

multiply 3;
id E([mu]?) * EC([nu]?) = 1/3*( -d_([mu], [nu]) +
  k`i'([mu])*k`i'([nu])/(`m')^2 );

#endif

.sort
#endprocedure

***********************************************************************

#procedure Shortest(foo)
argument `foo';
#call Small
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

#call Small

.sort

#if "`Scale'" != "1"
$pow = count_(<k1,1>,...,<k`Legs',1>);
argument inv;
  $pow = $pow - count_(<k1,1>,...,<k`Legs',1>);
endargument;
if( $pow != 0 ) multiply pow(`Scale', $pow/2);
#endif

#call Square

id [p1]?.[p2]? = abb([p1].[p2]);
id inv([x]?) = abb(1/[x]);
id e_([mu]?, [nu]?, [rho]?, [sig]?) = abb(e_([mu], [nu], [rho], [sig]));

b abb, A0, B0, B1, B00, B11, C0i, D0i, E0i, Conjugate, Den;
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

b abb, A0, B0, B1, B00, B11, C0i, D0i, E0i, Conjugate, Den;
print;
.end
#endprocedure

***********************************************************************

i [mu], [nu], [rho], [sig];
v [p1], [p2];
s [x], [y];

cf inv, abb, Simplify;
cf A0, B0, B1, B00, B11, C0i, D0i, E0i, Conjugate, Den;
t E, EC;

auto s FC;

.global

