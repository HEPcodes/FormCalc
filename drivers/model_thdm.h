* model_thdm.h
* declarations for model_thdm.F
* this file is part of FormCalc
* last modified 2 Oct 10 th


#include "model_sm.h"

	double precision Lambda5
	double precision Yuk1, Yuk2, Yuk3
	logical thdm_digest

	common /thdm_para/ Lambda5
	common /thdm_para/ Yuk1, Yuk2, Yuk3
	common /thdm_para/ thdm_digest
