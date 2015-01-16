* 1to2.F
* the kinematics for a 1 -> 2 process
* this file is part of FormCalc
* last modified 7 Feb 05 th


* some definitions for main.F
#define LEGS 3
#define LEGS_IN 1
#define TYPES TYPE1, TYPE2, TYPE3

#define HEADER "1to2.h"


************************************************************************
** DoIntegrated computes the integrated cross section.

	subroutine DoIntegrated
	implicit none

#include HEADER

	double precision result(NCOMP), error(NCOMP)

	print *, "sqrt(s) =", sqrtS
	call flush(6)

	call IntegratedCS(result, error)

1	format(TAG, G19.10, NCOMP G24.15, NCOMP G15.6)
	print 1, sqrtS, result, error
	call flush(6)
	end


************************************************************************
** DoDifferential computes the differential cross section.

	subroutine DoDifferential
	implicit none

#include HEADER

	double precision result(NCOMP), theta

	do theta = THETACUT, pi - (THETACUT), THETASTEP
	  print *, "sqrt(s) =", sqrtS, "  theta =", theta
	  call flush(6)

	  call DifferentialCS(result, cos(theta))

1	  format(TAG, G19.10, NCOMP G24.15)
	  print 1, theta, result
	  call flush(6)
	enddo
	end


************************************************************************
** SetProcess sets up process-dependent constants.

        subroutine SetProcess(hel, avgfac, sqrtSfrom, sqrtSto)
        implicit none
        integer hel
        double precision avgfac, sqrtSfrom, sqrtSto

#include HEADER

	helicities = hel
	preflux = 1/(2*(2*pi)**2*4)*avgfac

	if( sqrtSto .ge. sqrtSfrom ) then
* the factor 2 pi takes care of the azimuthal integration:
	  preflux = 2*pi*preflux
* make loop over sqrtS terminate immediately:
	  sqrtSto = -1
	endif
	end


************************************************************************
** SetEnergy initializes everything that depends only on sqrtS.

	subroutine SetEnergy(newsqrtS, *)
	implicit none
	double precision newsqrtS

#include HEADER
#include "util.h"

1	continue
	newsqrtS = MASS1
	if( newsqrtS .le. 0 ) newsqrtS = 1
	call ModelVarIni(newsqrtS, *999)
	if( abs(MASS1 - newsqrtS) .gt. 1D-14 ) goto 1
	if( newsqrtS .lt. MASS2 + MASS3 ) then
	  print *, "Below threshold"
	  return 1
	endif

	sqrtS = newsqrtS
	Pout = ThreeMom(sqrtS, dble(MASS2), dble(MASS3))

	flux = preflux*Pout/sqrtS**2

	call VecSet(1, dble(MASS1), 0D0, 0D0, 0D0, 1D0)

	call PREFIX CalcRenConst
	reset = .TRUE.
	return

999	return 1
	end


************************************************************************
** DifferentialCS actually invokes the generated function SquaredME.
** If necessary, soft-photon bremsstrahlung corrections are added.

	subroutine DifferentialCS(result, costh)
	implicit none
	double precision result(*), costh

#include HEADER
#include "util.h"
#include "renconst.h"

#ifdef BREMSSTRAHLUNG
	double precision SoftPhotonFactor
	external SoftPhotonFactor
#endif

	double precision sinth
	integer c

	sinth = sqrt(1 - costh**2)

	call VecSet(2, dble(MASS2), Pout, sinth, 0D0, costh)
	call VecSet(3, dble(MASS3), Pout, -sinth, 0D0, -costh)

	call SquaredME(result, helicities, reset)

	do c = 1, NCOMP
	  result(c) = result(c)*flux
	enddo

#ifdef BREMSSTRAHLUNG
	result(2) = result(2) + SoftPhotonFactor()*result(1)
#endif

#ifdef WF_RENORMALIZATION
	result(2) = result(2) + (WF_RENORMALIZATION)*result(1)
#endif
	end


************************************************************************
** IntegratedCS integrates DifferentialCS over the scattering angle.

	subroutine IntegratedCS(result, error)
	implicit none
	double precision result(*), error(*)

#include HEADER

	double precision costhcut
	integer c, neval, fail
	external DifferentialCS

	costhcut = cos(dble(THETACUT))
	call Patterson(NCOMP, -costhcut, costhcut, DifferentialCS,
     &    dble(ABSACCURACY), dble(RELACCURACY),
     &    neval, fail, result, error)

	if( fail .ne. 0 )
     &    print *, "Failed to reach the desired accuracy."

c	print *, neval, " evaluations used"
	end

