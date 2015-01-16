* HelicityME.frm
* the FORM part of the HelicityME function
* this file is part of FormCalc
* last modified 16 Jul 07 th


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

#call eiei

#do i = 1, `Legs'
#ifdef `k`i''
#call eiki
id k`i' = `k`i'';
.sort;
#endif
#enddo

#call eiki
#call kikj

#call Neglect

.sort

#if "`Scale'" != "1"
$pow = count_(<k1,1>,...,<k`Legs',1>);
if( $pow != 0 ) multiply pow(`Scale', $pow/2);
#endif

#call Square

id [p1]?.[p2]? = abb([p1].[p2]);
id e_([mu]?, [nu]?, [rho]?, [sig]?) = abb(e_([mu], [nu], [rho], [sig]));

moduleoption polyfun=abb;
.sort

normalize abb;
id abb(1) = 1;

b `Hels', abb;
.sort

collect Simplify, Simplify;
normalize Simplify;

#call Factor(Simplify)
#call InvSimplify(Simplify)

b `Hels';
print;
.store
#endprocedure

***********************************************************************

i [mu], [nu], [rho], [sig];
v [p1], [p2];
s [x], [y];

cf abb, Simplify;

auto s FC;

.global

