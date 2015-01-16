* IdentityMatrix.F
* set the n-by-n matrix A to the identity matrix
* this file is part of FormCalc
* last modified 14 Sep 04 th


	subroutine IdentityMatrix(n, A, lda)
	implicit none
	integer n, lda
	double complex A(lda,*)

	integer i, j

	do i = 1, n
	  do j = 1, n
	    A(j,i) = 0
	  enddo
	  A(i,i) = 1
	enddo
	end

