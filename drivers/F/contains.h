* contains.h
* inline versions of the util functions
* this file is part of FormCalc
* last modified 26 Oct 15 th


#ifndef CONTAINS_H
#define CONTAINS_H

#define SpiLV(iL,eL) (1-2*eL)*Vec(1+eL,1+eL,iL), Vec(2-eL,1+eL,iL)
#define SpiLB(iL,eL) (1-2*eL)*Vec(1+eL,2-eL,iL), Vec(2-eL,2-eL,iL)
#define SpiRV(eR,iR) Vec(1+eR,1+eR,iR), (1-2*eR)*Vec(2-eR,1+eR,iR)
#define SpiRB(eR,iR) Vec(1+eR,2-eR,iR), (1-2*eR)*Vec(2-eR,2-eR,iR)

#else

	contains

	function Pair(a, b)
	HelType Pair
	integer a, b
	Pair = .5D0*(
     &    Vec(1,1,a)*Vec(2,2,b) + Vec(2,2,a)*Vec(1,1,b) -
     &    Vec(1,2,a)*Vec(2,1,b) - Vec(2,1,a)*Vec(1,2,b) )
	end function

	function Eps_(a, b, c, d)
	HelType Eps_
	integer a, b, c, d
	Eps_ =
     &    (Vec(1,1,a)*Vec(2,2,b) - Vec(2,2,a)*Vec(1,1,b))*
     &    (Vec(2,1,c)*Vec(1,2,d) - Vec(1,2,c)*Vec(2,1,d))
	end function

	function Eps(a, b, c, d)
	HelType Eps
	integer a, b, c, d
	Eps = .25D0*(
     &    Eps_(a, b, c, d) + Eps_(c, d, a, b) -
     &    Eps_(a, c, b, d) - Eps_(b, d, a, c) +
     &    Eps_(a, d, b, c) + Eps_(b, c, a, d) )
	end function

	function SxS(l1,l2, r1,r2)
	HelType SxS, l1,l2, r1,r2
	SxS = l1*r1 + l2*r2
	end function

	function SxV1(l1,l2, a)
	HelType SxV1, l1,l2
	integer a
	SxV1 = l1*Vec(1,1,a) + l2*Vec(2,1,a)
	end function

	function SxV2(l1,l2, a)
	HelType SxV2, l1,l2
	integer a
	SxV2 = l2*Vec(2,2,a) + l1*Vec(1,2,a)
	end function

	function SxB1(l1,l2, a)
	HelType SxB1, l1,l2
	integer a
	SxB1 = l1*Vec(2,2,a) - l2*Vec(2,1,a)
	end function

	function SxB2(l1,l2, a)
	HelType SxB2, l1,l2
	integer a
	SxB2 = l2*Vec(1,1,a) - l1*Vec(1,2,a)
	end function

	function VxS1(a, r1,r2)
	HelType VxS1, r1,r2
	integer a
	VxS1 = Vec(1,1,a)*r1 + Vec(1,2,a)*r2
	end function

	function VxS2(a, r1,r2)
	HelType VxS2, r1,r2
	integer a
	VxS2 = Vec(2,1,a)*r1 + Vec(2,2,a)*r2
	end function

	function BxS1(a, r1,r2)
	HelType BxS1, r1,r2
	integer a
	BxS1 = Vec(2,2,a)*r1 - Vec(1,2,a)*r2
	end function

	function BxS2(a, r1,r2)
	HelType BxS2, r1,r2
	integer a
	BxS2 = Vec(1,1,a)*r2 - Vec(2,1,a)*r1
	end function

	function SxVxB1(l1,l2, a, b)
	HelType SxVxB1, l1,l2
	integer a, b
	SxVxB1 = SxB1(SxV1(l1,l2, a),SxV2(l1,l2, a), b)
	end function

	function SxVxB2(l1,l2, a, b)
	HelType SxVxB2, l1,l2
	integer a, b
	SxVxB2 = SxB2(SxV1(l1,l2, a),SxV2(l1,l2, a), b)
	end function

	function SxBxV1(l1,l2, a, b)
	HelType SxBxV1, l1,l2
	integer a, b
	SxBxV1 = SxV1(SxB1(l1,l2, a),SxB2(l1,l2, a), b)
	end function

	function SxBxV2(l1,l2, a, b)
	HelType SxBxV2, l1,l2
	integer a, b
	SxBxV2 = SxV2(SxB1(l1,l2, a),SxB2(l1,l2, a), b)
	end function

	function BxVxS1(b, a, r1,r2)
	HelType BxVxS1, r1,r2
	integer a, b
	BxVxS1 = BxS1(b, VxS1(a, r1,r2),VxS2(a, r1,r2))
	end function

	function BxVxS2(b, a, r1,r2)
	HelType BxVxS2, r1,r2
	integer b, a
	BxVxS2 = BxS2(b, VxS1(a, r1,r2),VxS2(a, r1,r2))
	end function

	function VxBxS1(b, a, r1,r2)
	HelType VxBxS1, r1,r2
	integer b, a
	VxBxS1 = VxS1(b, BxS1(a, r1,r2),BxS2(a, r1,r2))
	end function

	function VxBxS2(b, a, r1,r2)
	HelType VxBxS2, r1,r2
	integer b, a
	VxBxS2 = VxS2(b, BxS1(a, r1,r2),BxS2(a, r1,r2))
	end function

	function SxVxBxV1(l1,l2, a, b, c)
	HelType SxVxBxV1, l1,l2
	integer a, b, c
	SxVxBxV1 = SxBxV1(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c)
	end function

	function SxVxBxV2(l1,l2, a, b, c)
	HelType SxVxBxV2, l1,l2
	integer a, b, c
	SxVxBxV2 = SxBxV2(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c)
	end function

	function SxBxVxB1(l1,l2, a, b, c)
	HelType SxBxVxB1, l1,l2
	integer a, b, c
	SxBxVxB1 = SxVxB1(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c)
	end function

	function SxBxVxB2(l1,l2, a, b, c)
	HelType SxBxVxB2, l1,l2
	integer a, b, c
	SxBxVxB2 = SxVxB2(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c)
	end function

	function VxBxVxS1(c, b, a, r1,r2)
	HelType VxBxVxS1, r1,r2
	integer c, b, a
	VxBxVxS1 = VxBxS1(c, b, VxS1(a, r1,r2),VxS2(a, r1,r2))
	end function

	function VxBxVxS2(c, b, a, r1,r2)
	HelType VxBxVxS2, r1,r2
	integer c, b, a
	VxBxVxS2 = VxBxS2(c, b, VxS1(a, r1,r2),VxS2(a, r1,r2))
	end function

	function BxVxBxS1(c, b, a, r1,r2)
	HelType BxVxBxS1, r1,r2
	integer c, b, a
	BxVxBxS1 = BxVxS1(c, b, BxS1(a, r1,r2),BxS2(a, r1,r2))
	end function

	function BxVxBxS2(c, b, a, r1,r2)
	HelType BxVxBxS2, r1,r2
	integer c, b, a
	BxVxBxS2 = BxVxS2(c, b, BxS1(a, r1,r2),BxS2(a, r1,r2))
	end function

	function SxVxBxVxB1(l1,l2, a, b, c, d)
	HelType SxVxBxVxB1, l1,l2
	integer a, b, c, d
	SxVxBxVxB1 = SxBxVxB1(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c, d)
	end function

	function SxVxBxVxB2(l1,l2, a, b, c, d)
	HelType SxVxBxVxB2, l1,l2
	integer a, b, c, d
	SxVxBxVxB2 = SxBxVxB2(SxV1(l1,l2, a),SxV2(l1,l2, a), b, c, d)
	end function

	function SxBxVxBxV1(l1,l2, a, b, c, d)
	HelType SxBxVxBxV1, l1,l2
	integer a, b, c, d
	SxBxVxBxV1 = SxVxBxV1(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c, d)
	end function

	function SxBxVxBxV2(l1,l2, a, b, c, d)
	HelType SxBxVxBxV2, l1,l2
	integer a, b, c, d
	SxBxVxBxV2 = SxVxBxV2(SxB1(l1,l2, a),SxB2(l1,l2, a), b, c, d)
	end function

	function BxVxBxVxS1(d, c, b, a, r1,r2)
	HelType BxVxBxVxS1, r1,r2
	integer d, c, b, a
	BxVxBxVxS1 = BxVxBxS1(d, c, b, VxS1(a, r1,r2),VxS2(a, r1,r2))
	end function

	function BxVxBxVxS2(d, c, b, a, r1,r2)
	HelType BxVxBxVxS2, r1,r2
	integer d, c, b, a
	BxVxBxVxS2 = BxVxBxS2(d, c, b, VxS1(a, r1,r2),VxS2(a, r1,r2))
	end function

	function VxBxVxBxS1(d, c, b, a, r1,r2)
	HelType VxBxVxBxS1, r1,r2
	integer d, c, b, a
	VxBxVxBxS1 = VxBxVxS1(d, c, b, BxS1(a, r1,r2),BxS2(a, r1,r2))
	end function

	function VxBxVxBxS2(d, c, b, a, r1,r2)
	HelType VxBxVxBxS2, r1,r2
	integer d, c, b, a
	VxBxVxBxS2 = VxBxVxS2(d, c, b, BxS1(a, r1,r2),BxS2(a, r1,r2))
	end function

	function ChainV0(iL,eL, eR,iR)
	HelType ChainV0
	integer iL,eL, eR,iR
	ChainV0 = SxS(
     &    SpiLB(iL,eL),
     &    SpiRV(eR,iR) )
	end function

	function ChainB0(iL,eL, eR,iR)
	HelType ChainB0
	integer iL,eL, eR,iR
	ChainB0 = SxS(
     &    SpiLV(iL,eL),
     &    SpiRB(eR,iR) )
	end function

	function ChainV1(iL,eL, a, eR,iR)
	HelType ChainV1
	integer iL,eL, a, eR,iR
	ChainV1 = SxS(
     &    SxV1(SpiLB(iL,eL), a),
     &    SxV2(SpiLB(iL,eL), a),
     &    SpiRB(eR,iR) )
	end function

	function ChainB1(iL,eL, a, eR,iR)
	HelType ChainB1
	integer iL,eL, a, eR,iR
	ChainB1 = SxS(
     &    SxB1(SpiLV(iL,eL), a),
     &    SxB2(SpiLV(iL,eL), a),
     &    SpiRV(eR,iR) )
	end function

	function ChainV2(iL,eL, a, b, eR,iR)
	HelType ChainV2
	integer iL,eL, a, b, eR,iR
	ChainV2 = SxS(
     &    SxV1(SpiLB(iL,eL), a),
     &    SxV2(SpiLB(iL,eL), a),
     &    BxS1(b, SpiRV(eR,iR)),
     &    BxS2(b, SpiRV(eR,iR)) )
	end function

	function ChainB2(iL,eL, a, b, eR,iR)
	HelType ChainB2
	integer iL,eL, a, b, eR,iR
	ChainB2 = SxS(
     &    SxB1(SpiLV(iL,eL), a),
     &    SxB2(SpiLV(iL,eL), a),
     &    VxS1(b, SpiRB(eR,iR)),
     &    VxS2(b, SpiRB(eR,iR)) )
	end function

	function ChainV3(iL,eL, a, b, c, eR,iR)
	HelType ChainV3
	integer iL,eL, a, b, c, eR,iR
	ChainV3 = SxS(
     &    SxVxB1(SpiLB(iL,eL), a, b),
     &    SxVxB2(SpiLB(iL,eL), a, b),
     &    VxS1(c, SpiRB(eR,iR)),
     &    VxS2(c, SpiRB(eR,iR)) )
	end function

	function ChainB3(iL,eL, a, b, c, eR,iR)
	HelType ChainB3
	integer iL,eL, a, b, c, eR,iR
	ChainB3 = SxS(
     &    SxBxV1(SpiLV(iL,eL), a, b),
     &    SxBxV2(SpiLV(iL,eL), a, b),
     &    BxS1(c, SpiRV(eR,iR)),
     &    BxS2(c, SpiRV(eR,iR)) )
	end function

	function ChainV4(iL,eL, a, b, c, d, eR,iR)
	HelType ChainV4
	integer iL,eL, a, b, c, d, eR,iR
	ChainV4 = SxS(
     &    SxVxB1(SpiLB(iL,eL), a, b),
     &    SxVxB2(SpiLB(iL,eL), a, b),
     &    VxBxS1(c, d, SpiRV(eR,iR)),
     &    VxBxS2(c, d, SpiRV(eR,iR)) )
	end function

	function ChainB4(iL,eL, a, b, c, d, eR,iR)
	HelType ChainB4
	integer iL,eL, a, b, c, d, eR,iR
	ChainB4 = SxS(
     &    SxBxV1(SpiLV(iL,eL), a, b),
     &    SxBxV2(SpiLV(iL,eL), a, b),
     &    BxVxS1(c, d, SpiRB(eR,iR)),
     &    BxVxS2(c, d, SpiRB(eR,iR)) )
	end function

	function ChainV5(iL,eL, a, b, c, d, e, eR,iR)
	HelType ChainV5
	integer iL,eL, a, b, c, d, e, eR,iR
	ChainV5 = SxS(
     &    SxVxBxV1(SpiLB(iL,eL), a, b, c),
     &    SxVxBxV2(SpiLB(iL,eL), a, b, c),
     &    BxVxS1(d, e, SpiRB(eR,iR)),
     &    BxVxS2(d, e, SpiRB(eR,iR)) )
	end function

	function ChainB5(iL,eL, a, b, c, d, e, eR,iR)
	HelType ChainB5
	integer iL,eL, a, b, c, d, e, eR,iR
	ChainB5 = SxS(
     &    SxBxVxB1(SpiLV(iL,eL), a, b, c),
     &    SxBxVxB2(SpiLV(iL,eL), a, b, c),
     &    VxBxS1(d, e, SpiRV(eR,iR)),
     &    VxBxS2(d, e, SpiRV(eR,iR)) )
	end function

	function ChainV6(iL,eL, a, b, c, d, e, f, eR,iR)
	HelType ChainV6
	integer iL,eL, a, b, c, d, e, f, eR,iR
	ChainV6 = SxS(
     &    SxVxBxV1(SpiLB(iL,eL), a, b, c),
     &    SxVxBxV2(SpiLB(iL,eL), a, b, c),
     &    BxVxBxS1(d, e, f, SpiRV(eR,iR)),
     &    BxVxBxS2(d, e, f, SpiRV(eR,iR)) )
	end function

	function ChainB6(iL,eL, a, b, c, d, e, f, eR,iR)
	HelType ChainB6
	integer iL,eL, a, b, c, d, e, f, eR,iR
	ChainB6 = SxS(
     &    SxBxVxB1(SpiLV(iL,eL), a, b, c),
     &    SxBxVxB2(SpiLV(iL,eL), a, b, c),
     &    VxBxVxS1(d, e, f, SpiRB(eR,iR)),
     &    VxBxVxS2(d, e, f, SpiRB(eR,iR)) )
	end function

	function ChainV7(iL,eL, a, b, c, d, e, f, g, eR,iR)
	HelType ChainV7
	integer iL,eL, a, b, c, d, e, f, g, eR,iR
	ChainV7 = SxS(
     &    SxVxBxVxB1(SpiLB(iL,eL), a, b, c, d),
     &    SxVxBxVxB2(SpiLB(iL,eL), a, b, c, d),
     &    VxBxVxS1(e, f, g, SpiRB(eR,iR)),
     &    VxBxVxS2(e, f, g, SpiRB(eR,iR)) )
	end function

	function ChainB7(iL,eL, a, b, c, d, e, f, g, eR,iR)
	HelType ChainB7
	integer iL,eL, a, b, c, d, e, f, g, eR,iR
	ChainB7 = SxS(
     &    SxBxVxBxV1(SpiLV(iL,eL), a, b, c, d),
     &    SxBxVxBxV2(SpiLV(iL,eL), a, b, c, d),
     &    BxVxBxS1(e, f, g, SpiRV(eR,iR)),
     &    BxVxBxS2(e, f, g, SpiRV(eR,iR)) )
	end function

	function ChainV8(iL,eL, a, b, c, d, e, f, g, h, eR,iR)
	HelType ChainV8
	integer iL,eL, a, b, c, d, e, f, g, h, eR,iR
	ChainV8 = SxS(
     &    SxVxBxVxB1(SpiLB(iL,eL), a, b, c, d),
     &    SxVxBxVxB2(SpiLB(iL,eL), a, b, c, d),
     &    VxBxVxBxS1(e, f, g, h, SpiRV(eR,iR)),
     &    VxBxVxBxS2(e, f, g, h, SpiRV(eR,iR)) )
	end function

	function ChainB8(iL,eL, a, b, c, d, e, f, g, h, eR,iR)
	HelType ChainB8
	integer iL,eL, a, b, c, d, e, f, g, h, eR,iR
	ChainB8 = SxS(
     &    SxBxVxBxV1(SpiLV(iL,eL), a, b, c, d),
     &    SxBxVxBxV2(SpiLV(iL,eL), a, b, c, d),
     &    BxVxBxVxS1(e, f, g, h, SpiRB(eR,iR)),
     &    BxVxBxVxS2(e, f, g, h, SpiRB(eR,iR)) )
	end function

#endif

