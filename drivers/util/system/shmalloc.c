/*
	shmalloc.c
		(De)allocation of IPC shared memory
		for use in the DifferentialCS integrand function
		this file is part of FormCalc
		last modified 4 Nov 13 th

Example usage:

	integer*8 i
* base must be long enough to hold two pointers,
* i.e. 16 bytes on a 64-bit platform
	RealType base(2)

* allocate an array of 100 RealTypes:
	call shmalloc(i, base, 100, 8)

	call DifferentialCS
* in DifferentialCS:
*	access base(i), base(i+1), ..., base(i+99)

* allocated shared memory is not automatically freed at exit
* (can be freed manually with ipcrm), so make sure you do:
	call shmfree(base)

*/


#include <string.h>
#include <sys/shm.h>
#include <assert.h>

#if NOUNDERSCORE
#define shmalloc_ shmalloc
#define shmfree_ shmfree
#endif


typedef long long int memindex;

typedef struct {
  void *addr;
  int id;
} shminfo;

static inline memindex PtrDiff(const void *a, const void *b) {
  return (char *)a - (char *)b;
}


void shmalloc_(memindex *index, shminfo *base, const int *n, const int *size)
{
  base->id = shmget(IPC_PRIVATE, (*n + 1)**size - 1, IPC_CREAT | 0600);
  assert(base->id != -1);
  base->addr = shmat(base->id, NULL, 0);
  assert(base->addr != (void *)-1);
  *index = PtrDiff(base->addr + *size - 1, base)/(long)*size;
}


void shmfree_(shminfo *base)
{
  shmdt(base->addr);
  shmctl(base->id, IPC_RMID, NULL);
}

