* parton.h
* setup for a single partonic process
* is included by partonic.h for every partonic process
* this file is part of FormCalc
* last modified 18 Jun 15

#undef CPall
#if defined CPini || defined CPxs
#define CPall
#endif

#ifndef PARTON1
#define PARTON1 -1
#endif

#ifndef PARTON2
#define PARTON2 -1
#endif

CPdecl	external SQUAREDME_FUNC

CPall	pid = PID
CPall	parton1 = PARTON1
CPall	parton2 = PARTON2
CPall	type = Generic(ARG_Ptyp,JOIN_OCT)
CPall	hel = iand(helmask, Generic(ARG_Phel,JOIN_HEL))

CPini	avgfac(pid) = Re(COLOURFACTOR)/(IDENTICALFACTOR*spin_df(hel))
#ifdef UNPOLARIZED_DIRAC_FERMIONS
CPini	avgfac(pid) = avgfac(pid)*Generic(ARG_Pfac,JOIN_MUL)
#endif

#ifdef BREMSSTRAHLUNG
CPini	call SetArray(charge(1,pid), Charge(ARG_RE,JOIN_SEQ))
#endif

CPini	PartonicIni(SQUAREDME_FUNC)

CPxs	PartonicXS(SQUAREDME_FUNC)


#undef PID
#undef PARTON1
#undef PARTON2

