* multigauss.F
* multi-dimensional adaptive Gauss-Kronrod integration
* based on code by S. Berge
* this file is part of FormCalc
* last modified 12 Sep 01 th

* NCOMP is the number of components of the integrand vector.
*#define NCOMP 2


************************************************************************
** multigauss integrates a vector using an adaptive Gauss-Kronrod
** quadrature. It attempts to evaluate the integral with a relative
** error less than accuracy, but breaks off if the number of intervals
** into which it has to subdivide the integration region becomes larger
** than nmax.

	options /recursive

	subroutine multigauss(result, accuracy, para, intv, nmax, nodes)
	implicit none
	double precision result(NCOMP), accuracy
	double precision para(-1:NCOMP, 0:*)
	integer intv(3), nmax, nodes

#define X(n) para(-1, n)
#define ERR(n) para(0, n)
#define VAL(n) para(NCOMP, n)

	double precision err
	integer i, k, c, n

	do c = 0, NCOMP
	  para(c, 0) = 0
	enddo
	n = intv(1)
	do i = 1, n
	  call kronrod(para, i, nodes)
	enddo

* as long as the sum of errors of all intervals is larger than
* accuracy, subdivide the interval with the largest error and repeat.

	do while(VAL(0) .ne. 0 .and.
     +           ERR(0) .gt. abs(accuracy*VAL(0)))
	  if(n .ge. nmax) then
	    print *,
     +      "Warning: multigauss failed to reach the desired accuracy."
	    print *, "Remaining relative error: ", ERR(0)/VAL(0)
	    goto 999
	  endif

* find interval with largest error
	  err = ERR(1)
	  i = 1
	  do k = 1, n
	    if(ERR(k) .gt. err) then
	      err = ERR(k)
	      i = k
	    endif
	  enddo

* make room for one more subdivision
	  do k = n, i, -1
	    do c = -1, NCOMP
	      para(c, k + 1) = para(c, k)
	    enddo
	  enddo
	  n = n + 1

	  do c = 0, NCOMP
	    para(c, 0) = para(c, 0) - para(c, i)
	  enddo
	  X(i) = X(i - 1) + .5D0*(X(i) - X(i - 1))
	  call kronrod(para, i, nodes)
	  call kronrod(para, i + 1, nodes)
	enddo

999	continue
	do c = 1, NCOMP
	  result(c) = para(c, 0)
	enddo

	intv(3) = intv(3) + n - intv(1)
	intv(2) = max(intv(2), n)

* join intervals with sufficient accuracy
	err = VAL(0)*accuracy**2
	i = 1
	do while(i .lt. n)
	  if(ERR(i) + ERR(i + 1) .lt. err) then
	    n = n - 1
	    do k = i, n
	      X(k) = X(k + 1)
	    enddo
	  else
	    i = i + 1
	  endif
	enddo
	intv(1) = n
	end 


************************************************************************
** kronrod integrates the interval X(i - 1)...X(i) using a
** nodes/(2 nodes + 1) Gauss-Kronrod rule. See KronrodPoints.m for the
** computation of the array krdata.

	options /recursive

	subroutine kronrod(para, i, nodes)
	implicit none
	double precision para(-1:NCOMP, 0:*)
	integer i, nodes

	double precision res(NCOMP, 6), range, mid, x
	integer c, k, p

	double precision krdata(105)
	data krdata /
     +    0.34644298189013636168107712823159978D0,
     +    0.33998104358485626480266575910324469D0,
     +    0.64028621749630998240468902315749202D0,
     +    0.32694918960145162955845946561731919D0,
     +    0.26679834045228444803277062841785566D0,
     +    0.65214515486254614262693605077800059D0,
     +    0.86113631159405257522394648889280951D0,
     +    0.97656025073757311153450535936991963D0,
     +    0.17005360533572272680273885329620659D0,
     +    0.062977373665473014765492488552818676D0,
     +    0.34785484513745385737306394922199941D0,
     +    0.24107258017346476191063599297275917D0,
     +    0.23861918608319690863050172168071194D0,
     +    0.46311821247530461215675836401917663D0,
     +    0.23377086411699440662283572598899837D0,
     +    0.2132096522719622791628941635168893D0,
     +    0.46791393457269104738987034398955099D0,
     +    0.66120938646626451366139959501990535D0,
     +    0.82137334086502794004564983424395025D0,
     +    0.18107199432313761518699209331551194D0,
     +    0.13732060463444692308714987253378181D0,
     +    0.36076157304813860756983351383771611D0,
     +    0.93246951420315202781230155449399461D0,
     +    0.98870320261267885750464595171218508D0,
     +    0.083694440446906626132845603482411108D0,
     +    0.030396154119819768851964544676027885D0,
     +    0.17132449237917034504029614217273289D0,
     +    0.18444640574469164352897095570564293D0,
     +    0.18343464249564980493947614236018398D0,
     +    0.36070109792813195719254862229689144D0,
     +    0.18140002506803464306174852517255044D0,
     +    0.17207060855521131185729488020385709D0,
     +    0.36268378337836198296515044927719561D0,
     +    0.52553240991632898581773904918924635D0,
     +    0.67235407094515867715631073809283105D0,
     +    0.15665260616818840049024808848696874D0,
     +    0.1362631092551722152623387452545062D0,
     +    0.31370664587788728733796220198660131D0,
     +    0.79666647741362673959155393647583044D0,
     +    0.89412090684745642194836101753825132D0,
     +    0.11164637082683961322210815893394148D0,
     +    0.082482298931358330688625193445607895D0,
     +    0.22238103445337447054435599442624088D0,
     +    0.96028985649753623168356086856947299D0,
     +    0.9933798758817161559358880690196708D0,
     +    0.049439395002139308500363969446996895D0,
     +    0.01782238332071035515278696120274979D0,
     +    0.10122853629037625915253135430996219D0,
     +    0.1494455540029169056649364683898212D0,
     +    0.14887433898163121088482600112971998D0,
     +    0.29439286270146019813112660310386557D0,
     +    0.14773910490133849137484151597206805D0,
     +    0.14277593857706008079709427313871706D0,
     +    0.29552422471475287017389299465133833D0,
     +    0.43339539412924719079926594316578416D0,
     +    0.56275713466860468333900009927269414D0,
     +    0.13470921731147332592805400177170683D0,
     +    0.12349197626206585107795810983107416D0,
     +    0.26926671930999635509122692156946935D0,
     +    0.67940956829902440623432736511487358D0,
     +    0.78081772658641689706371757834504238D0,
     +    0.10938715880229764189921059032580496D0,
     +    0.093125454583697605535065465083366344D0,
     +    0.21908636251598204399553493422816319D0,
     +    0.86506336668898451073209668842349305D0,
     +    0.93015749135570822600120718005950835D0,
     +    0.075039674810919952767043140916190009D0,
     +    0.054755896574351996031381300244580176D0,
     +    0.14945134915058059314577633965769733D0,
     +    0.97390652851717172007796401208445205D0,
     +    0.99565716302580808073552728068900285D0,
     +    0.032558162307964727478818972459389761D0,
     +    0.011694638867371874278064396062192048D0,
     +    0.066671344308688137593568809893331793D0,
     +    0.12555689390547433530429613286007933D0,
     +    0.12523340851146891547244136946385313D0,
     +    0.2485057483204692762677909603627181D0,
     +    0.12458416453615607343731247320922899D0,
     +    0.12162630352394838324609975809130963D0,
     +    0.24914704581340278500056243604295121D0,
     +    0.36783149899818019375269153664371756D0,
     +    0.48133945047815709293594361501883226D0,
     +    0.11671205350175682629358074530573005D0,
     +    0.11002260497764407263590739874224961D0,
     +    0.23349253653835480876084989892487806D0,
     +    0.58731795428661744729670241894053428D0,
     +    0.68405989547005589394492910034115448D0,
     +    0.10164973227906027771568877049122751D0,
     +    0.091549468295049210528171939739614172D0,
     +    0.20316742672306592174906445580979838D0,
     +    0.76990267419430468703689383321281808D0,
     +    0.84355812416115324479214188505983926D0,
     +    0.079920275333601701493392609529783355D0,
     +    0.067250907050839930304940940047316074D0,
     +    0.16007832854334622633465252954335907D0,
     +    0.90411725637047485667846586611909619D0,
     +    0.95053779594312129654906019513161863D0,
     +    0.053697017607756251228889163320458185D0,
     +    0.038915230469299477115089632285862896D0,
     +    0.10693932599531843096025471819399622D0,
     +    0.98156063424671925069054909014928082D0,
     +    0.99693392252959542691235023725838527D0,
     +    0.023036084038982232591084580367968981D0,
     +    0.0082577114331683957576939224392115776D0,
     +    0.04717533638651182719461596148501706D0 /

	p = nodes*(5*nodes - 6)/8 - 6
	range = .5D0*(X(i) - X(i - 1))
	mid = X(i - 1) + range
	call sample(res(1, 2), mid)
	do c = 1, NCOMP
	  res(c, 1) = 0
	  res(c, 2) = res(c, 2)*krdata(p)
	enddo
	do k = 1, nodes/2
	  x = range*krdata(p + 1)
	  call sample(res(1, 3), mid + x)
	  call sample(res(1, 4), mid - x)
	  x = range*krdata(p + 2)
	  call sample(res(1, 5), mid + x)
	  call sample(res(1, 6), mid - x)
	  do c = 1, NCOMP
	    x = res(c, 3) + res(c, 4)
	    res(c, 2) = res(c, 2) + x*krdata(p + 3)
     +        + (res(c, 5) + res(c, 6))*krdata(p + 4)
	  enddo
	  res(NCOMP, 1) = res(NCOMP, 1) + x*krdata(p + 5)
	  p = p + 5
	enddo

* value of integral:
	do c = 1, NCOMP
	  para(c, i) = range*res(c, 2)
	  para(c, 0) = para(c, 0) + para(c, i)
	enddo

* error estimate:
	ERR(i) = abs(range*(res(NCOMP, 2) - res(NCOMP, 1)))
* the original Kronrod error estimate includes the following line,
* but for most cases that's just too conservative:
*	ERR(i) = (200*ERR(i))**1.5D0
	ERR(0) = ERR(0) + ERR(i)

	end

