* Jacobi.f
* computes the variables for an elementary Jacobi rotation
* code adapted from the "Handbook" routines for complex M12
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of FormCalc
* last modified 8 Dec 05 th


	subroutine Jacobi(delta, s, tau, diff, M12, M12sq)
	implicit none
	double precision delta, tau, diff, M12sq
	double complex s, M12

	double precision h, invc

	h = 2/(diff + sign(sqrt(diff**2 + 4*M12sq), diff))

* delta = DCONJG(M12)*tan(phi) is the shift in the diagonal elements,
* i.e. M11 -> M11 + delta, M22 -> M22 - delta
	delta = h*M12sq

* invc = 1/cos(phi)
	invc = sqrt(delta*h + 1)

* s = sin(phi)
	s = h/invc*M12

* tau = DCONJG(sin(phi))*tan(phi/2) = 1 - cos(phi)
	tau = h*delta/(invc*(invc + 1))
	end

