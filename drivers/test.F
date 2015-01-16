#ifdef __AVX__
#define SIMD2 2
#elif defined __SSE3__
#define SIMD2 1
#else
#define SIMD2 0
#endif

	program test
#if defined SIMD && SIMD == SIMD2
#else
	print *, SIMD2
#endif
	end
