* ReverseEigen.F
* recover the original matrix A from the eigenvalues and
* the transformation matrix U by forming A = U^+ diag(d) U.
* this file is part of FormCalc
* last modified 20 Oct 05 th


	subroutine ReverseEigen(n, A, lda, d, U, ldu)
	implicit none
	integer n, lda, ldu
	double complex A(lda,*), U(ldu,*)
	double precision d(*)

	integer i, j, k
	double complex s

	do j = 1, n
	  do i = 1, n
	    s = 0
	    do k = 1, n
	      s = s + d(k)*DCONJG(U(k,i))*U(k,j)
	    enddo
	    A(i,j) = s
	  enddo
	enddo
	end

