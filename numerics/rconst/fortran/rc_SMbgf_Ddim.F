* Note: this file was generated automatically by RenConst.m.
* Changes by hand will be overwritten the next time RenConst.m is run.

	subroutine calc_renconst
	implicit character (a-s,u-z)
	implicit double complex (t)

#include "model.h"
#include "looptools.h"
#include "renconst.h"

	integer IndexDelta
	external IndexDelta

	integer g1, Gen1

        dMZsq1 = Alfa/(48.D0*CW2**2*Pi*SW2)*
     -   (-8*CW2**3*MZ2 + 3*CW2*dble(A0(MH2)) + 
     -     6*CW2*(9*CW2**2 - 2*CW2*SW2 + SW2**2)*dble(A0(MW2)) + 
     -     3*CW2*dble(A0(MZ2)) + 12*MW2*dble(B0(MZ2,MH2,MZ2)) + 
     -     24*CW2*(MW2 - 4*CW2**2*MZ2)*dble(B0(MZ2,MW2,MW2)) + 
     -     36*CW2*dble(B00(MZ2,0.D0,0.D0)) - 
     -     12*CW2*dble(B00(MZ2,MH2,MZ2)) - 
     -     12*CW2*(9*CW2**2 - 2*CW2*SW2 + SW2**2)*
     -      dble(B00(MZ2,MW2,MW2)) - 
     -     18*CW2*MZ2*dble(B1(MZ2,0.D0,0.D0)))

	do Gen1 = 1, 3
        dMZsq1 = dMZsq1 - 
     -   Alfa/(24.D0*CW2*Pi*SW2)*
     -    (3*(1 - 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(A0(MLE2(Gen1))) + 
     -      (1 + 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(A0(MQD2(Gen1))) + 
     -      (1 - 8*CW2 + 16*CW2**2 + 16*SW2**2)*
     -       dble(A0(MQU2(Gen1))) - 
     -      6*(1 - 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(B00(MZ2,MLE2(Gen1),MLE2(Gen1))) - 
     -      2*(1 + 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(B00(MZ2,MQD2(Gen1),MQD2(Gen1))) + 
     -      (-2 + 16*CW2 - 32*CW2**2 - 32*SW2**2)*
     -       dble(B00(MZ2,MQU2(Gen1),MQU2(Gen1))) + 
     -      3*MZ2*(1 - 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(B1(MZ2,MLE2(Gen1),MLE2(Gen1))) + 
     -      MZ2*(1 + 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(B1(MZ2,MQD2(Gen1),MQD2(Gen1))) + 
     -      MZ2*(1 - 8*CW2 + 16*CW2**2 + 16*SW2**2)*
     -       dble(B1(MZ2,MQU2(Gen1),MQU2(Gen1))) + 
     -      3*(-1 + 2*CW2 + 2*SW2)**2*
     -       (dble(B0(MZ2,MLE2(Gen1),MLE2(Gen1)))*MLE2(Gen1)) + 
     -      (1 + 2*CW2 + 2*SW2)**2*
     -       (dble(B0(MZ2,MQD2(Gen1),MQD2(Gen1)))*MQD2(Gen1)) + 
     -      (-1 + 4*CW2 + 4*SW2)**2*
     -       (dble(B0(MZ2,MQU2(Gen1),MQU2(Gen1)))*MQU2(Gen1)))
	enddo

#ifdef DEBUG
	print *, 'dMZsq1 =', dMZsq1
#endif

        dMWsq1 = Alfa/(48.D0*CW2*Pi*SW2)*
     -   (8*CW2*MW2*(-3 + 2*CW2 + 2*SW2) + 3*CW2*dble(A0(MH2)) + 
     -     30*CW2*dble(A0(MW2)) + 
     -     3*CW2*(1 + 8*CW2)*dble(A0(MZ2)) - 
     -     48*CW2*MW2*SW2*dble(B0(MW2,0.D0,MW2)) + 
     -     12*CW2*MW2*dble(B0(MW2,MH2,MW2)) - 
     -     12*MW2*(7*CW2**2 - SW2**2 + CW2*(-1 + 2*SW2))*
     -      dble(B0(MW2,MW2,MZ2)) - 
     -     96*CW2*SW2*dble(B00(MW2,0.D0,MW2)) - 
     -     12*CW2*dble(B00(MW2,MH2,MW2)) - 
     -     12*CW2*(1 + 8*CW2)*dble(B00(MW2,MZ2,MW2)))

	do Gen1 = 1, 3
        dMWsq1 = dMWsq1 - 
     -   Alfa/(4.D0*Pi*SW2)*(dble(A0(MLE2(Gen1))) + 
     -      3*dble(A0(MQD2(Gen1))) - 
     -      2*dble(B00(MW2,0.D0,MLE2(Gen1))) - 
     -      6*dble(B00(MW2,MQU2(Gen1),MQD2(Gen1))) + 
     -      MW2*dble(B1(MW2,0.D0,MLE2(Gen1))) + 
     -      3*MW2*dble(B1(MW2,MQU2(Gen1),MQD2(Gen1))) + 
     -      3*dble(B0(MW2,MQD2(Gen1),MQU2(Gen1)))*MQU2(Gen1))
	enddo

#ifdef DEBUG
	print *, 'dMWsq1 =', dMWsq1
#endif

        dMHsq1 = Alfa/(32.D0*CW2**2*MW2*Pi*SW2)*
     -   (-4*MW2*((2 + 6*CW2**2)*MW2 + CW2*MZ2) + 
     -     dble(B0(MH2,MZ2,MZ2))*
     -      (12*MW2**2 - 4*CW2*MW2*MH2 + CW2**2*MH2**2) + 
     -     3*CW2**2*(MH2*dble(A0(MH2))) + 
     -     2*CW2**2*((MH2 + 6*MW2)*dble(A0(MW2))) + 
     -     9*CW2**2*(MH2**2*dble(B0(MH2,MH2,MH2))) + 
     -     CW2*(dble(A0(MZ2))*(6*MW2 + CW2*MH2)) + 
     -     2*CW2**2*(dble(B0(MH2,MW2,MW2))*
     -        (MH2**2 + 12*MW2**2 - 4*MW2*MH2)))

	do Gen1 = 1, 3
        dMHsq1 = dMHsq1 - 
     -   Alfa/(4.D0*MW2*Pi*SW2)*
     -    (dble(A0(MLE2(Gen1)))*MLE2(Gen1) + 
     -      MH2*dble(B1(MH2,MLE2(Gen1),MLE2(Gen1)))*MLE2(Gen1) + 
     -      2*dble(B0(MH2,MLE2(Gen1),MLE2(Gen1)))*
     -       MLE2(Gen1)**2 + 3*dble(A0(MQD2(Gen1)))*MQD2(Gen1) + 
     -      3*MH2*dble(B1(MH2,MQD2(Gen1),MQD2(Gen1)))*
     -       MQD2(Gen1) + 6*dble(B0(MH2,MQD2(Gen1),MQD2(Gen1)))*
     -       MQD2(Gen1)**2 + 3*dble(A0(MQU2(Gen1)))*MQU2(Gen1) + 
     -      3*MH2*dble(B1(MH2,MQU2(Gen1),MQU2(Gen1)))*
     -       MQU2(Gen1) + 6*dble(B0(MH2,MQU2(Gen1),MQU2(Gen1)))*
     -       MQU2(Gen1)**2)
	enddo

#ifdef DEBUG
	print *, 'dMHsq1 =', dMHsq1
#endif

        dZAA1 = Alfa/(6.D0*Pi)*
     -   (1 + 12*dble(B0(0.D0,MW2,MW2)) + 
     -     18*dble(DB00(0.D0,MW2,MW2)))

	do Gen1 = 1, 3
        dZAA1 = dZAA1 + Alfa/(3.D0*Pi)*
     -    (3*dble(B1(0.D0,MLE2(Gen1),MLE2(Gen1))) + 
     -      dble(B1(0.D0,MQD2(Gen1),MQD2(Gen1))) + 
     -      4*dble(B1(0.D0,MQU2(Gen1),MQU2(Gen1))) - 
     -      6*dble(DB00(0.D0,MLE2(Gen1),MLE2(Gen1))) - 
     -      2*dble(DB00(0.D0,MQD2(Gen1),MQD2(Gen1))) - 
     -      8*dble(DB00(0.D0,MQU2(Gen1),MQU2(Gen1))))
	enddo

#ifdef DEBUG
	print *, 'dZAA1 =', dZAA1
#endif

	dZZA1 = 0

#ifdef DEBUG
	print *, 'dZZA1 =', dZZA1
#endif

        dTad1 = 1/(64.D0*CW2*MW*Pi**2)*
     -   (EL*(4*MW2*(2*CW2*MW2 + MZ2) + 
     -        dble(A0(MZ2))*(-6*MW2 - CW2*MH2) - 
     -        3*CW2*(MH2*dble(A0(MH2))) - 
     -        2*CW2*((MH2 + 6*MW2)*dble(A0(MW2)))))/SW

	do Gen1 = 1, 3
        dTad1 = dTad1 + 1/(8.D0*MW*Pi**2)*
     -    (EL*(dble(A0(MLE2(Gen1)))*MLE2(Gen1) + 
     -         3*dble(A0(MQD2(Gen1)))*MQD2(Gen1) + 
     -         3*dble(A0(MQU2(Gen1)))*MQU2(Gen1)))/SW
	enddo

#ifdef DEBUG
	print *, 'dTad1 =', dTad1
#endif

	dSWsq1 = CW2*(-(1/MW2*dMWsq1) + 1/MZ2*dMZsq1)

#ifdef DEBUG
	print *, 'dSWsq1 =', dSWsq1
#endif

	dCWsq1 = -dSWsq1

#ifdef DEBUG
	print *, 'dCWsq1 =', dCWsq1
#endif

	dSW1 = dSWsq1/(2.D0*SW)

#ifdef DEBUG
	print *, 'dSW1 =', dSW1
#endif

	dZe1 = -(1/(2.D0*CW)*(dZZA1*SW + CW*dZAA1))

#ifdef DEBUG
	print *, 'dZe1 =', dZe1
#endif

	dZAZ1 = 2/CW*dCWsq1/SW

#ifdef DEBUG
	print *, 'dZAZ1 =', dZAZ1
#endif

	dZZZ1 = dZAA1 + 1/CW2*dCWsq1 - 1/SW2*dCWsq1

#ifdef DEBUG
	print *, 'dZZZ1 =', dZZZ1
#endif

	dZW1 = dZAA1 - 1/SW2*dCWsq1

#ifdef DEBUG
	print *, 'dZW1 =', dZW1
#endif

	dZH1 = dZW1 + 1/MW2*dMWsq1

#ifdef DEBUG
	print *, 'dZH1 =', dZH1
#endif

        dWFZ1 = 1/(24.D0*CW2**2*Pi*SW2)*
     -   (48*Alfa*CW2**3*dble(B0(MZ2,MW2,MW2)) + 
     -     9*Alfa*CW2*dble(B1(MZ2,0.D0,0.D0)) - 
     -     6*Alfa*MW2*dble(DB0(MZ2,MH2,MZ2)) + 
     -     12*Alfa*CW2*(-MW2 + 4*CW2**2*MZ2)*
     -      dble(DB0(MZ2,MW2,MW2)) - 
     -     18*Alfa*CW2*dble(DB00(MZ2,0.D0,0.D0)) + 
     -     6*Alfa*CW2*dble(DB00(MZ2,MH2,MZ2)) + 
     -     6*Alfa*CW2*(9*CW2**2 - 2*CW2*SW2 + SW2**2)*
     -      dble(DB00(MZ2,MW2,MW2)) + 
     -     9*Alfa*CW2*MZ2*dble(DB1(MZ2,0.D0,0.D0)) + 
     -     4*CW2**2*(Alfa*CW2 - 6*Pi*SW2*dZZZ1))

	do Gen1 = 1, 3
        dWFZ1 = dWFZ1 + Alfa/(24.D0*CW2*Pi*SW2)*
     -    (3*(1 - 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(B1(MZ2,MLE2(Gen1),MLE2(Gen1))) + 
     -      (1 + 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(B1(MZ2,MQD2(Gen1),MQD2(Gen1))) + 
     -      (1 - 8*CW2 + 16*CW2**2 + 16*SW2**2)*
     -       dble(B1(MZ2,MQU2(Gen1),MQU2(Gen1))) - 
     -      6*(1 - 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(DB00(MZ2,MLE2(Gen1),MLE2(Gen1))) - 
     -      2*(1 + 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(DB00(MZ2,MQD2(Gen1),MQD2(Gen1))) + 
     -      (-2 + 16*CW2 - 32*CW2**2 - 32*SW2**2)*
     -       dble(DB00(MZ2,MQU2(Gen1),MQU2(Gen1))) + 
     -      3*MZ2*(1 - 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(DB1(MZ2,MLE2(Gen1),MLE2(Gen1))) + 
     -      MZ2*(1 + 4*CW2 + 4*CW2**2 + 4*SW2**2)*
     -       dble(DB1(MZ2,MQD2(Gen1),MQD2(Gen1))) + 
     -      MZ2*(1 - 8*CW2 + 16*CW2**2 + 16*SW2**2)*
     -       dble(DB1(MZ2,MQU2(Gen1),MQU2(Gen1))) + 
     -      3*(-1 + 2*CW2 + 2*SW2)**2*
     -       (dble(DB0(MZ2,MLE2(Gen1),MLE2(Gen1)))*MLE2(Gen1)) + 
     -      (1 + 2*CW2 + 2*SW2)**2*
     -       (dble(DB0(MZ2,MQD2(Gen1),MQD2(Gen1)))*MQD2(Gen1)) + 
     -      (-1 + 4*CW2 + 4*SW2)**2*
     -       (dble(DB0(MZ2,MQU2(Gen1),MQU2(Gen1)))*MQU2(Gen1)))
	enddo

#ifdef DEBUG
	print *, 'dWFZ1 =', dWFZ1
#endif

        dWFW1 = 1/(12.D0*CW2*Pi*SW2)*
     -   (24*Alfa*CW2*SW2*dble(B0(MW2,0.D0,MW2)) + 
     -     24*Alfa*CW2**2*dble(B0(MW2,MW2,MZ2)) + 
     -     12*Alfa*CW2*MW2*SW2*dble(DB0(MW2,0.D0,MW2)) - 
     -     3*Alfa*CW2*MW2*dble(DB0(MW2,MH2,MW2)) + 
     -     3*Alfa*MW2*(7*CW2**2 - SW2**2 + CW2*(-1 + 2*SW2))*
     -      dble(DB0(MW2,MW2,MZ2)) + 
     -     24*Alfa*CW2*SW2*dble(DB00(MW2,0.D0,MW2)) + 
     -     3*Alfa*CW2*dble(DB00(MW2,MH2,MW2)) + 
     -     3*Alfa*CW2*(1 + 8*CW2)*dble(DB00(MW2,MZ2,MW2)) + 
     -     2*CW2*(Alfa - 6*Pi*SW2*dZW1))

	do Gen1 = 1, 3
        dWFW1 = dWFW1 + Alfa/(4.D0*Pi*SW2)*
     -    (dble(B1(MW2,0.D0,MLE2(Gen1))) + 
     -      3*dble(B1(MW2,MQU2(Gen1),MQD2(Gen1))) - 
     -      2*dble(DB00(MW2,0.D0,MLE2(Gen1))) - 
     -      6*dble(DB00(MW2,MQU2(Gen1),MQD2(Gen1))) + 
     -      MW2*dble(DB1(MW2,0.D0,MLE2(Gen1))) + 
     -      3*MW2*dble(DB1(MW2,MQU2(Gen1),MQD2(Gen1))) + 
     -      3*dble(DB0(MW2,MQD2(Gen1),MQU2(Gen1)))*MQU2(Gen1))
	enddo

#ifdef DEBUG
	print *, 'dWFW1 =', dWFW1
#endif

	dZchi1 = dZH1

#ifdef DEBUG
	print *, 'dZchi1 =', dZchi1
#endif

	dZphi1 = dZH1

#ifdef DEBUG
	print *, 'dZphi1 =', dZphi1
#endif


	do g1 = 1, 3
	dMf1(1,g1) = 0

#ifdef DEBUG
	print *, 'dMf1(1,',g1,') =', dMf1(1,g1)
#endif

        dZfL1(1,g1,g1) = 
     -  Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -   (MW2 + 2*CW2*MW2 + 2*MW2*dble(B1(0.D0,0.D0,MZ2)) + 
     -     2*CW2*(dble(B1(0.D0,MLE2(g1),MW2))*(2*MW2 + MLE2(g1))))

#ifdef DEBUG
	print *, 'dZfL1(1,',g1,',',g1,') =', dZfL1(1,g1,g1)
#endif

	dZfR1(1,g1,g1) = 0

#ifdef DEBUG
	print *, 'dZfR1(1,',g1,',',g1,') =', dZfR1(1,g1,g1)
#endif

        dMf1(2,g1) = -(Alfa/(32.D0*CW2*MW2*Pi*SW2)*
     -     ((MW2*(1 + 2*CW2 + 4*SW2 - 8*SW2**2) + 
     -         2*dble(B1(MLE2(g1),MLE2(g1),MZ2))*
     -          (MW2 - 4*MW2*SW2 + 8*MW2*SW2**2 + CW2*MLE2(g1))+
     -           dble(B0(MLE2(g1),MZ2,MLE2(g1)))*
     -          (16*MW2*SW2*(-1 + 2*SW2) + 2*CW2*MLE2(g1)) - 
     -         2*CW2*(dble(B0(MLE2(g1),MH2,MLE2(g1)))*
     -            MLE2(g1)) + 
     -         2*CW2*(dble(B1(MLE2(g1),MLE2(g1),MH2))*
     -            MLE2(g1)) + 
     -         2*CW2*(dble(B1(MLE2(g1),0.D0,MW2))*
     -            (2*MW2 + MLE2(g1))))*MLE(g1)))

#ifdef DEBUG
	print *, 'dMf1(2,',g1,') =', dMf1(2,g1)
#endif

        dZfL1(2,g1,g1) = 
     -  Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -   (MW2*(2*CW2 + (1 - 2*SW2)**2) + 
     -     4*CW2*MW2*dble(B1(MLE2(g1),0.D0,MW2)) + 
     -     dble(B1(MLE2(g1),MLE2(g1),MZ2))*
     -      (2*MW2*(1 - 2*SW2)**2 + CW2*MLE2(g1)) + 
     -     CW2*(dble(B1(MLE2(g1),MLE2(g1),MH2))*MLE2(g1)) - 
     -     2*CW2*(dble(DB0(MLE2(g1),MH2,MLE2(g1)))*
     -        MLE2(g1)**2) + 
     -     2*CW2*(dble(DB1(MLE2(g1),MLE2(g1),MH2))*
     -        MLE2(g1)**2) + 
     -     2*CW2*(dble(DB1(MLE2(g1),0.D0,MW2))*MLE2(g1)*
     -        (2*MW2 + MLE2(g1))) + 
     -     2*dble(DB1(MLE2(g1),MLE2(g1),MZ2))*
     -      (MW2 - 4*MW2*SW2 + 8*MW2*SW2**2 + CW2*MLE2(g1))*
     -      MLE2(g1) + 2*dble(DB0(MLE2(g1),MZ2,MLE2(g1)))*
     -      (8*MW2*SW2*(-1 + 2*SW2) + CW2*MLE2(g1))*MLE2(g1))

#ifdef DEBUG
	print *, 'dZfL1(2,',g1,',',g1,') =', dZfL1(2,g1,g1)
#endif

        dZfR1(2,g1,g1) = 
     -  Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -   (4*MW2*SW2**2 + dble(B1(MLE2(g1),MLE2(g1),MZ2))*
     -      (8*MW2*SW2**2 + CW2*MLE2(g1)) + 
     -     2*CW2*(dble(B1(MLE2(g1),0.D0,MW2))*MLE2(g1)) + 
     -     CW2*(dble(B1(MLE2(g1),MLE2(g1),MH2))*MLE2(g1)) - 
     -     2*CW2*(dble(DB0(MLE2(g1),MH2,MLE2(g1)))*
     -        MLE2(g1)**2) + 
     -     2*CW2*(dble(DB1(MLE2(g1),MLE2(g1),MH2))*
     -        MLE2(g1)**2) + 
     -     2*CW2*(dble(DB1(MLE2(g1),0.D0,MW2))*MLE2(g1)*
     -        (2*MW2 + MLE2(g1))) + 
     -     2*dble(DB1(MLE2(g1),MLE2(g1),MZ2))*
     -      (MW2 - 4*MW2*SW2 + 8*MW2*SW2**2 + CW2*MLE2(g1))*
     -      MLE2(g1) + 2*dble(DB0(MLE2(g1),MZ2,MLE2(g1)))*
     -      (8*MW2*SW2*(-1 + 2*SW2) + CW2*MLE2(g1))*MLE2(g1))

#ifdef DEBUG
	print *, 'dZfR1(2,',g1,',',g1,') =', dZfR1(2,g1,g1)
#endif

        dMf1(3,g1) = -(Alfa/(288.D0*CW2*MW2*Pi*SW2)*
     -     ((MW2*(9 + 18*CW2 + 24*SW2 - 32*SW2**2) + 
     -         36*CW2*(dble(B0(MQU2(g1),MW2,MQD2(g1)))*
     -            MQD2(g1)) + 
     -         dble(B0(MQU2(g1),MZ2,MQU2(g1)))*
     -          (32*MW2*SW2*(-3 + 4*SW2) + 18*CW2*MQU2(g1)) + 
     -         dble(B1(MQU2(g1),MQU2(g1),MZ2))*
     -          (2*MW2*(9 - 24*SW2 + 32*SW2**2) + 
     -            18*CW2*MQU2(g1)) - 
     -         18*CW2*(dble(B0(MQU2(g1),MH2,MQU2(g1)))*
     -            MQU2(g1)) + 
     -         18*CW2*(dble(B1(MQU2(g1),MQU2(g1),MH2))*
     -            MQU2(g1)) + 
     -         18*CW2*(dble(B1(MQU2(g1),MQD2(g1),MW2))*
     -            (2*MW2 + MQD2(g1) + MQU2(g1))))*MQU(g1)))

#ifdef DEBUG
	print *, 'dMf1(3,',g1,') =', dMf1(3,g1)
#endif

        dZfL1(3,g1,g1) = 
     -  Alfa/(144.D0*CW2*MW2*Pi*SW2)*
     -   (MW2*(18*CW2 + (3 - 4*SW2)**2) + 
     -     18*CW2*(dble(B1(MQU2(g1),MQD2(g1),MW2))*
     -        (2*MW2 + MQD2(g1))) + 
     -     dble(B1(MQU2(g1),MQU2(g1),MZ2))*
     -      (2*MW2*(3 - 4*SW2)**2 + 9*CW2*MQU2(g1)) + 
     -     9*CW2*(dble(B1(MQU2(g1),MQU2(g1),MH2))*MQU2(g1)) + 
     -     36*CW2*(dble(DB0(MQU2(g1),MW2,MQD2(g1)))*MQD2(g1)*
     -        MQU2(g1)) + dble(DB0(MQU2(g1),MZ2,MQU2(g1)))*
     -      (32*MW2*SW2*(-3 + 4*SW2)*MQU2(g1) + 
     -        18*CW2*MQU2(g1)**2) - 
     -     18*CW2*(dble(DB0(MQU2(g1),MH2,MQU2(g1)))*
     -        MQU2(g1)**2) + 
     -     18*CW2*(dble(DB1(MQU2(g1),MQU2(g1),MH2))*
     -        MQU2(g1)**2) + 
     -     18*CW2*(dble(DB1(MQU2(g1),MQD2(g1),MW2))*MQU2(g1)*
     -        (2*MW2 + MQD2(g1) + MQU2(g1))) + 
     -     2*dble(DB1(MQU2(g1),MQU2(g1),MZ2))*
     -      (MW2*(9 - 24*SW2 + 32*SW2**2) + 9*CW2*MQU2(g1))*
     -      MQU2(g1))

#ifdef DEBUG
	print *, 'dZfL1(3,',g1,',',g1,') =', dZfL1(3,g1,g1)
#endif

        dZfR1(3,g1,g1) = 
     -  Alfa/(144.D0*CW2*MW2*Pi*SW2)*
     -   (16*MW2*SW2**2 + dble(B1(MQU2(g1),MQU2(g1),MZ2))*
     -      (32*MW2*SW2**2 + 9*CW2*MQU2(g1)) + 
     -     18*CW2*(dble(B1(MQU2(g1),MQD2(g1),MW2))*MQU2(g1)) + 
     -     9*CW2*(dble(B1(MQU2(g1),MQU2(g1),MH2))*MQU2(g1)) + 
     -     36*CW2*(dble(DB0(MQU2(g1),MW2,MQD2(g1)))*MQD2(g1)*
     -        MQU2(g1)) + dble(DB0(MQU2(g1),MZ2,MQU2(g1)))*
     -      (32*MW2*SW2*(-3 + 4*SW2)*MQU2(g1) + 
     -        18*CW2*MQU2(g1)**2) - 
     -     18*CW2*(dble(DB0(MQU2(g1),MH2,MQU2(g1)))*
     -        MQU2(g1)**2) + 
     -     18*CW2*(dble(DB1(MQU2(g1),MQU2(g1),MH2))*
     -        MQU2(g1)**2) + 
     -     18*CW2*(dble(DB1(MQU2(g1),MQD2(g1),MW2))*MQU2(g1)*
     -        (2*MW2 + MQD2(g1) + MQU2(g1))) + 
     -     2*dble(DB1(MQU2(g1),MQU2(g1),MZ2))*
     -      (MW2*(9 - 24*SW2 + 32*SW2**2) + 9*CW2*MQU2(g1))*
     -      MQU2(g1))

#ifdef DEBUG
	print *, 'dZfR1(3,',g1,',',g1,') =', dZfR1(3,g1,g1)
#endif

        dMf1(4,g1) = -(Alfa/(288.D0*CW2*MW2*Pi*SW2)*
     -     ((MW2*(9 + 18*CW2 + 12*SW2 - 8*SW2**2) + 
     -         dble(B0(MQD2(g1),MZ2,MQD2(g1)))*
     -          (16*MW2*SW2*(-3 + 2*SW2) + 18*CW2*MQD2(g1)) + 
     -         dble(B1(MQD2(g1),MQD2(g1),MZ2))*
     -          (MW2*(18 - 24*SW2 + 16*SW2**2) + 18*CW2*MQD2(g1))
     -           - 18*CW2*(dble(B0(MQD2(g1),MH2,MQD2(g1)))*
     -            MQD2(g1)) + 
     -         18*CW2*(dble(B1(MQD2(g1),MQD2(g1),MH2))*
     -            MQD2(g1)) + 
     -         36*CW2*(dble(B0(MQD2(g1),MW2,MQU2(g1)))*
     -            MQU2(g1)) + 
     -         18*CW2*(dble(B1(MQD2(g1),MQU2(g1),MW2))*
     -            (2*MW2 + MQD2(g1) + MQU2(g1))))*MQD(g1)))

#ifdef DEBUG
	print *, 'dMf1(4,',g1,') =', dMf1(4,g1)
#endif

        dZfL1(4,g1,g1) = 
     -  Alfa/(144.D0*CW2*MW2*Pi*SW2)*
     -   (MW2*(18*CW2 + (3 - 2*SW2)**2) + 
     -     dble(B1(MQD2(g1),MQD2(g1),MZ2))*
     -      (2*MW2*(3 - 2*SW2)**2 + 9*CW2*MQD2(g1)) + 
     -     9*CW2*(dble(B1(MQD2(g1),MQD2(g1),MH2))*MQD2(g1)) + 
     -     dble(DB0(MQD2(g1),MZ2,MQD2(g1)))*
     -      (16*MW2*SW2*(-3 + 2*SW2)*MQD2(g1) + 
     -        18*CW2*MQD2(g1)**2) - 
     -     18*CW2*(dble(DB0(MQD2(g1),MH2,MQD2(g1)))*
     -        MQD2(g1)**2) + 
     -     18*CW2*(dble(DB1(MQD2(g1),MQD2(g1),MH2))*
     -        MQD2(g1)**2) + 
     -     36*CW2*(dble(DB0(MQD2(g1),MW2,MQU2(g1)))*MQD2(g1)*
     -        MQU2(g1)) + 18*CW2*
     -      (dble(B1(MQD2(g1),MQU2(g1),MW2))*(2*MW2 + MQU2(g1)))+
     -       18*CW2*(dble(DB1(MQD2(g1),MQU2(g1),MW2))*MQD2(g1)*
     -        (2*MW2 + MQD2(g1) + MQU2(g1))) + 
     -     2*dble(DB1(MQD2(g1),MQD2(g1),MZ2))*
     -      (MW2*(9 - 12*SW2 + 8*SW2**2) + 9*CW2*MQD2(g1))*
     -      MQD2(g1))

#ifdef DEBUG
	print *, 'dZfL1(4,',g1,',',g1,') =', dZfL1(4,g1,g1)
#endif

        dZfR1(4,g1,g1) = 
     -  Alfa/(144.D0*CW2*MW2*Pi*SW2)*
     -   (4*MW2*SW2**2 + dble(B1(MQD2(g1),MQD2(g1),MZ2))*
     -      (8*MW2*SW2**2 + 9*CW2*MQD2(g1)) + 
     -     9*CW2*(dble(B1(MQD2(g1),MQD2(g1),MH2))*MQD2(g1)) + 
     -     18*CW2*(dble(B1(MQD2(g1),MQU2(g1),MW2))*MQD2(g1)) + 
     -     dble(DB0(MQD2(g1),MZ2,MQD2(g1)))*
     -      (16*MW2*SW2*(-3 + 2*SW2)*MQD2(g1) + 
     -        18*CW2*MQD2(g1)**2) - 
     -     18*CW2*(dble(DB0(MQD2(g1),MH2,MQD2(g1)))*
     -        MQD2(g1)**2) + 
     -     18*CW2*(dble(DB1(MQD2(g1),MQD2(g1),MH2))*
     -        MQD2(g1)**2) + 
     -     36*CW2*(dble(DB0(MQD2(g1),MW2,MQU2(g1)))*MQD2(g1)*
     -        MQU2(g1)) + 18*CW2*
     -      (dble(DB1(MQD2(g1),MQU2(g1),MW2))*MQD2(g1)*
     -        (2*MW2 + MQD2(g1) + MQU2(g1))) + 
     -     2*dble(DB1(MQD2(g1),MQD2(g1),MZ2))*
     -      (MW2*(9 - 12*SW2 + 8*SW2**2) + 9*CW2*MQD2(g1))*
     -      MQD2(g1))

#ifdef DEBUG
	print *, 'dZfR1(4,',g1,',',g1,') =', dZfR1(4,g1,g1)
#endif

	enddo

	end

