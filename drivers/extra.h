#if 0
	extra.h
	extra user preprocessor defs
	(included by both Fortran and C)
	this file is part of FormCalc
	last modified 19 Nov 13 th
#endif


#if 0
Possibly some wave-function renormalization
(e.g. if calculating in the background-field method)

Wave-function renormalization is added as in
  Mloop += (WF_RENORMALIZATION)*Mtree

For example, the background-field formulation of the SM
(SMbgf.mod) needs wf.ren. for external gauge bosons, as in

#define WF_RENORMALIZATION (nW*dWFW1 + nZ*dWFZ1)
#endif

