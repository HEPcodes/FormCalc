* HelicityME.frm
* the FORM part of the HelicityME function
* this file is part of FormCalc
* last modified 20 Nov 02 th


i [mu], [nu], [rho], [sig];
v [p1], [p2];
s [x], [y];

s scale;
cf abb, Simplify;

auto s FC;

.global


#procedure DotSimplify(momsubst, moresimp)
`momsubst'
#call eiki
`moresimp'
.sort
#endprocedure


#procedure Factor(foo)
factarg `foo';
repeat id `foo'([x]?, [y]?, ?a) = `foo'([x]) * `foo'([y], ?a);
id `foo'([x]?number_) = [x];
id `foo'([x]?symbol_) = [x];
#endprocedure


#procedure Emit
contract 0;
#call eiei
#call MomSimplify(DotSimplify)

#call TrivialSubst

#if `Scaled'
$pow = count_(<k1,1>,...,<k`Legs',1>);
multiply scale^$pow;
#endif

id [p1]?.[p2]? = abb([p1].[p2]);
id e_([mu]?, [nu]?, [rho]?, [sig]?) = abb(e_([mu], [nu], [rho], [sig]));

moduleoption polyfun=abb;
.sort

normalize abb;
id abb(1) = 1;

b `Hels', abb;
.sort

collect Simplify;
normalize Simplify;

#call Factor(Simplify)
#call MandelSimplify(Simplify)

b `Hels';
print;
.store
#endprocedure

