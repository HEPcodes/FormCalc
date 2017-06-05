* CompAZ.F
* compute the spectrum of a backscattered laser photon at a
* photon collider
* Code adapted from http://hep.fuw.edu.pl/u/zarnecki/compaz/compaz.f
* last modified 10 Jun 16 th

#include "types.h"


	subroutine CompAZnew(res, x, E, icase)
	implicit none
	RealType res, x, E
	integer icase

* Input parameters:
* -----------------
*
*   x             Fraction of beam energy carried by photon  ( 0 < x < 1 )
*
*   E             Electron beam energy [GeV]
*
*   icase         Choice of output distribution:
*
*   icase = 0     Full spectra - sum of three contributions
*           1     Main contribution from direct Compton process
*           2     Contribution from scattering with 2 photons
*           3     Contribution from scattering on secondary electrons
*
*          10     Polarization of the full spectra
*          11     Polarization of the 1st component (direct Compton)
*          12     Polarization of the 2nd component (two photon)
*          13     Polarization of the 3rd component (secondary electron)
*
*          20     Effective z parameter
*          21     maximum x value for direct Compton
*          22     maximum x value for scattering with 2 photons

* Parameters of the distribution, partly fitted to the spectra from V. Telnov
	RealType Br2, gpow, polL, pole, pole2

* - Laser parameter
#define z0fun(E) (18.08800D0*E)

* - Nonlinear corrections
#define xifun(E) (.06594662D0 + .7060851D0*E)

* - Total normalization
#define totfun(E) (.7115863D0 - .6776124D0*E)

* - 2nd component fraction
#define c2fun(xi) (.5540019D0*(1 - exp(-37.38912D0*xi**2)))

* - 3rd component fraction
#define c3fun(E) (.7257064D0 + 1.517959D0*E)

* - Corrections to Compton spectra
	parameter (Br2 = 3.315865D0/8, gpow = 1.387423D0/2)

* - Electron/laser polarisation
	parameter (polL = 1, pole = .85D0, pole2 = 0)

	RealType Ee, z0, xi, z
	RealType xmax, xmax2, g
	RealType tot, tot2, tot3
	RealType gflux1, gflux2, gflux3
	RealType pol1, pol2, pol3
	RealType g1, g2
	RealType ymin, yspan, ymid, y, w

	logical DoPolar

	integer n

	RealType p, pe, pL
	RealType gflux, gpol, norm, CompAZc, CompAZp
	RealType gfun, pfun

	integer nodes
	parameter (nodes = 12)
	RealType gaussdata(2,nodes)

	data gaussdata /
     &    .99518721999702136017999740970073681D0,
     &    .012341229799987199546805667070037292D0,
     &    .97472855597130949819839199300816906D0,
     &    .028531388628933663181307815951878286D0,
     &    .93827455200273275852364900170872145D0,
     &    .044277438817419806168602748211338229D0,
     &    .88641552700440103421315434198219676D0,
     &    .059298584915436780746367758500108585D0,
     &    .82000198597390292195394987266974521D0,
     &    .073346481411080305734033615253116518D0,
     &    .74012419157855436424382810309997843D0,
     &    .086190161531953275917185202983742667D0,
     &    .64809365193697556925249578691074763D0,
     &    .097618652104113888269880664464247154D0,
     &    .54542147138883953565837561721837237D0,
     &    .10744427011596563478257734244660622D0,
     &    .43379350762604513848708423191334971D0,
     &    .11550566805372560135334448390678356D0,
     &    .31504267969616337438679329131981024D0,
     &    .12167047292780339120446315347626243D0,
     &    .19111886747361630915863982075706963D0,
     &    .12583745634682829612137538251118369D0,
     &    .064056892862605626085043082624745039D0,
     &    .12793819534675215697405616522469537D0/

* Photon flux, as given in hep-ph/9508397
	gflux(x, z, p) = 2 - x + x/(1 - x)*
     &    (1 + (2/z*x/(1 - x) - 1)*(2/z + p*(2 - x)) - 2/z)

* Photon polarized flux
	gpol(x, z, pe, pL) =
     &    pe*x*(1/(1 - x) + (2*x/(1 - x)/z - 1)**2) +
     &    pL*(1 - x + 1/(1 - x))*(2*x/(1 - x)/z - 1)

* Flux normalisation (corrected!  formula in paper is wrong)
	norm(z, p) = .5D0*(1 + 16/z + (p*z**2 - 1)/(1 + z)**2) +
     &    (1 - (4/z + p)*(1 + 2/z))*log(1 + z) + 2*p

* 'Ideal' spectrum of compton backscattered photons
	CompAZc(x, z, p) = gflux(x, z, p)/norm(z, p)

* Spectrum of compton backscattered photons
	CompAZp(x, z, pe, pL) = gpol(x, z, pe, pL)/gflux(x, z, pe*pL)

	gfun(y) = log(y*z + 1)/(log(z + 1)*y) *
     &    CompAZc(x/y, y*z, pole2*polL)*
     &    CompAZc(1 - y, z, pole*polL)

	pfun(y) = CompAZp(x/y, y*z, pole2, polL)


* Check input flag
	if( .NOT. ((icase .ge. 0  .and. icase .le. 3 ) .or.
     &             (icase .ge. 10 .and. icase .le. 13) .or.
     &             (icase .ge. 20 .and. icase .le. 22)) )
     &    stop "CompAZ: Invalid icase flag"

* Check input energy
	if( E .lt. 50 .or. E .gt. 500 )
     &    stop "CompAZ: Beam energy outside range 50..500"

* Polarisation calculation flag

	DoPolar = icase .ge. 10 .and. icase .le. 13

* Set parameter values for given energy
* -------------------------------------

* Beam energy in TeV
	Ee = e/1000D0

* Nominal value of z = 4*Ee*El/me^2
	z0 = z0fun(Ee)

* correction nonlinear effects
	xi = xifun(Ee)

* effective z
	z = z0/(1 + xi)

* Maximum x values for different contributions
	xmax = z/(1 + z)
	xmax2 = 2*z/(1 + 2*z)

	res = z
	if( icase .eq. 20 ) return
	res = xmax
	if( icase .eq. 21 ) return
	res = xmax2
	if( icase .eq. 22 ) return

* Leading Compton term
* --------------------

	gflux1 = 0
	pol1 = 0
	tot = totfun(Ee)
	if( x .gt. 0 .and. x .le. xmax .and. tot .gt. 0 ) then
	  g = dim(z/x, z + 1)
	  gflux1 = tot*CompAZc(x, z, pole*polL)*exp(-Br2*g)
	  if( DoPolar ) pol1 = CompAZp(x, z, pole, polL)
	endif

	res = gflux1
	if( icase .eq. 1 ) return
	res = pol1
	if( icase .eq. 11) return

* Scattering on 2 photons
* -----------------------

	gflux2 = 0
	pol2 = 0
	tot2 = tot*c2fun(xi)
	if( x .gt. 0 .and. x .lt. xmax2 .and. tot2 .gt. 0 ) then
	  g = dim(2*z/x, 2*z + 1)
	  gflux2 = tot2*CompAZc(x, 2*z, pole*polL)*exp(-Br2*g)*g**gpow
	  if( DoPolar ) pol2 = CompAZp(max(x, xmax), 2*z, pole, polL)
	endif

	res = gflux2
	if( icase .eq. 2 ) return
	res = pol2
	if( icase .eq. 12 ) return

* Electron rescattering
* ---------------------

	gflux3 = 0
	pol3 = 0
	tot3 = tot*c3fun(Ee)
	if( x .gt. 0 .and. x .le. xmax .and. tot3 .gt. 0 ) then
	  ymin = max(.5D0*x*(1 + sqrt(4/(x*z) + 1)), 1/(z + 1))
	  yspan = .5D0*(1 - ymin)
	  if( yspan .gt. 0 ) then
	    ymid = ymin + yspan
	    do n = 1, nodes
	      y = yspan*gaussdata(1,n)
	      w = yspan*gaussdata(2,n)
	      g1 = gfun(ymid - y)
	      g2 = gfun(ymid + y)
	      gflux3 = gflux3 + w*(g1 + g2)
	      if( DoPolar ) pol3 = pol3 +
     &          w*(g1*pfun(ymid - y) + g2*pfun(ymid + y))
	    enddo
	    gflux3 = tot3*gflux3
	    pol3 = tot3*pol3
	  endif
	endif

	res = gflux3
	if( icase .eq. 3 ) return
	if( icase .eq. 13 ) then
	  if( res .ne. 0 ) res = pol3/res
	  return
	endif

	res = gflux1 + gflux2 + gflux3
	if( icase .eq. 0 ) return
	if( icase .eq. 10 ) then
	  if( res .ne. 0 ) res = (gflux1*pol1 + gflux2*pol2 + pol3)/res
	  return
	endif
	res = 0
	end

