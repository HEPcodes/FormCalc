* cuba.F
* chooser for the Cuba routines
* this file is part of FormCalc
* last modified 17 Feb 05 th


#define PATTERSON 0
#define VEGAS 1
#define SUAVE 2
#define DIVONNE 3
#define CUHRE 4

#define FLAGS_LAST 4


	subroutine Cuba(ndim, integrand, result, error)
	implicit none
	integer ndim
	external integrand
	double precision result(*), error(*)

	integer nregions, neval, fail, c
	double precision prob(NCOMP)

	character*9 name(0:4)
	data name(PATTERSON) /"Patterson"/
	data name(VEGAS) /"Vegas"/
	data name(SUAVE) /"Suave"/
	data name(DIVONNE) /"Divonne"/
	data name(CUHRE) /"Cuhre"/

	if( ndim .eq. 1 ) then

	  call Patterson(NCOMP, 0D0, 1D0, integrand,
     &      dble(ABSACCURACY), dble(RELACCURACY),
     &      neval, fail, result, error)
	  nregions = 1
	  print *, "Patterson integration results:"

	else

#if METHOD == VEGAS

	  call vegas(ndim, NCOMP, integrand,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      VERBOSE, MINEVAL, MAXEVAL,
     &      NSTART, NINCREASE,
     &      neval, fail, result, error, prob)
	  nregions = 1
	  print *, "Vegas integration results:"

#elif METHOD == SUAVE

	  call suave(ndim, NCOMP, integrand,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      VERBOSE + FLAGS_LAST, MINEVAL, MAXEVAL,
     &      NNEW, dble(FLATNESS),
     &      nregions, neval, fail, result, error, prob)
	  print *, "Suave integration results:"

#elif METHOD == DIVONNE

	  call divonne(ndim, NCOMP, integrand,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      VERBOSE, MINEVAL, MAXEVAL,
     &      KEY1, KEY2, KEY3, MAXPASS,
     &      dble(BORDER), dble(MAXCHISQ), dble(MINDEVIATION),
     &      0, NDIM, 0, 0, 0,
     &      nregions, neval, fail, result, error, prob)
	  print *, "Divonne integration results:"

#else

	  call cuhre(ndim, NCOMP, integrand,
     &      dble(RELACCURACY), dble(ABSACCURACY),
     &      VERBOSE + FLAGS_LAST, MINEVAL, MAXEVAL,
     &      KEY,
     &      nregions, neval, fail, result, error, prob)
	  print *, "Cuhre integration results:"

#endif
	endif

	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail

	if( fail .gt. 0 ) then
	  print *, "Failed to reach the desired accuracy."
	else if( fail .lt. 0 ) then
	  print *, "Integration error."
	endif

100	format(I2, G24.15, " +- ", G24.15, "    p = ", F5.3)
	print 100, (c, result(c), error(c), prob(c), c = 1, NCOMP)
	end

