* this is pave.frm, the Form part of FormCalc
* last modified 24 Aug 98 th
* Note: this program can also run with Form 1, but will possibly
* leave spinor chains not fully simplified

.sort

s r2, CW, SW, C4, S4, C2, S2, EL, a2;

id r2^2 = 2;
id r2^-2 = 1/2;
id CW^-4 = C4;
id SW^-4 = S4;
id CW^-2 = C2;
id SW^-2 = S2;
id EL^4 = 16 * Pi^2 * a2;

#call mandelstam{}

.sort

i I1, I2, I3, I4, C;
s M1, M2, M3, M4, X1, X2, unknown;
v P1, P2, P3;
f ga, ga5, fc;
cf DEN, eq, neq, repl, hide, pivot, fme;
cf d0, c0, b0, a0, pv, A0, B0m, B1m, B00m, B11m, pave3, pave4;
set bcd: b0, c0, d0;


* >>> Dirac algebra for open fermion chains:

#procedure diracsimplify()
#if 'VERSION_' > 1
  repeat id,disorder, ga(C?, I1?) * ga(C?, I2?) =
    2 * d_(I1, I2) - ga(C, I2) * ga(C, I1);
#endif
  id ga(C?, P1?) * ga(C?, P1?) = P1.P1;
  id ga(C?, I1?) * ga(C?, I1?) = d_(I1, I1);
  repeat;
    id ga(C?, P2?) * ga(C?, I1?) * pivot(C?, P1?, P2?) =
      (2 * d_(P2, I1) - ga(C, I1) * ga(C, P2)) * pivot(C, P1, P2);
  endrepeat;
  id ga(C?, P1?) * Spinor(P1?, M1?) = M1 * Spinor(P1, M1);
  repeat;
    id ga(C?, I2?) * ga(C?, P1?) * pivot(C?, P1?, P2?) =
      (2 * d_(P1, I2) - ga(C, P1) * ga(C, I2)) * pivot(C, P1, P2);
  endrepeat;
  id Spinor(P1?, M1?) * ga(C?, P1?) = -M1 * Spinor(P1, M1);
  id Spinor(P1?, M1?) * ga5(C?) * ga(C?, P1?) =
    M1 * Spinor(P1, M1) * ga5(C);
#endprocedure

id gi_(C?) = 1;
id Spinor(?) * gi_(C?) * Spinor(??) = Spinor(.) * Spinor(..);

if( count(g_, 1) > 0 );
  id g_(C?, I1?) = fc(C, ga(C, I1));
  id g5_(C?) = fc(C, ga5(C));
  repeat;
    id fc(C?, X1?) * fc(C?, X2?) = fc(C, X1 * X2);
  endrepeat;
  id Spinor(P1?, M1?) * fc(C?, X1?) * Spinor(P2?, M2?) =
    pivot(C, P1, P2) * Spinor(P1, M1) * X1 * Spinor(P2, M2);
  id fc(C?, X1?) = X1;
  repeat;
#call diracsimplify{}
  endrepeat;
endif;

.sort


* >>> cancel q^2's in the numerator

id q1.q1 * b0(P1?, M1?, M2?) =
  a0(M2) + M1 * b0(P1, M1, M2);

repeat;
  if( match(q1.q1) );
    id,once, q1.q1 * c0(P1?, P2?, M1?, M2?, M3?) =
      repl(q1 - P1) * b0(P2 - P1, M2, M3) +
        M1 * c0(P1, P2, M1, M2, M3);
    id,once, q1.q1 * d0(P1?, P2?, P3?, M1?, M2?, M3?, M4?) =
      repl(q1 - P1) * c0(P2 - P1, P3 - P1, M2, M3, M4) +
        M1 * d0(P1, P2, P3, M1, M2, M3, M4);
    repeat;
      id q1.q1 * repl(P1?) = hide(P1, P1) * repl(P1);
      id d_(q1, I1?) * repl(P1?) = hide(P1, I1) * repl(P1);
      id q1.P2? * repl(P1?) = hide(P1, P2) * repl(P1);
      id e_(q1, I1?, I2?, I3?) * repl(P1?) =
        hide(P1, I1, I2, I3) * repl(P1);
      id g_(C?, q1) * repl(P1?) = hide(0, C, P1) * repl(P1);
    endrepeat;
    id repl(?) = 1;
    id hide(0, C?, I1?) = g_(C, I1);
    id hide(I1?, I2?, I3?, I4?) = e_(I1, I2, I3, I4);
    id hide(I1?, I2?) = d_(I1, I2);
  endif;
endrepeat;

.sort


* >>> tensor reduction:

repeat;
  id,once, q1.q1 * pv?bcd(?) = pv(I4, I4, .);
  sum I4;
endrepeat;
repeat;
  id,once, q1.P1? * pv?bcd(?) = d_(P1, I4) * pv(I4, .);
  sum I4;
endrepeat;
repeat;
  id,once, d_(q1, I4?) * pv?bcd(?) = pv(I4, .);
  sum I4;
endrepeat;
repeat;
  id,once, ga(C?, q1) * pv?bcd(?) = ga(C, I4) * pv(I4, .);
  sum I4;
endrepeat;
repeat;
  id,once, e_(q1, I1?, I2?, I3?) * pv?bcd(?) =
    pv(I4, .) * e_(I4, I1, I2, I3);
  sum I4;
endrepeat;

#if 'VERSION_' > 1
b a0, b0, c0, d0;
.sort
keep brackets;
#else
.sort
#endif

#if 'DIM' == 4
* add local terms for dimred/CDR

id d0(I1?, I2?, I3?, I4?, P1?, P2?, P3?, M1?, M2?, M3?, M4?) =
  d0(I1, I2, I3, I4, P1, P2, P3, M1, M2, M3, M4) -
  eq(e_(I1, I2, I3, I4)) *
    5/144 * (d_(I1, I2) * d_(I3, I4) +
             d_(I1, I3) * d_(I2, I4) +
             d_(I1, I4) * d_(I2, I3)) +
  1/8 * (eq(I1, I2) * neq(I3, I4) + eq(I1, I3) * neq(I2, I4) +
         eq(I1, I4) * neq(I2, I3) + eq(I2, I3) * neq(I1, I4) +
         eq(I2, I4) * neq(I1, I3) + eq(I3, I4) * neq(I1, I2));

id eq(I1?, I1?) = 1;
id neq(I1?, I1?) = 0;
id neq(I1?, I2?) = d_(I1, I2);

id c0(I1?, I2?, I3?, P1?, P2?, M1?, M2?, M3?) =
  c0(I1, I2, I3, P1, P2, M1, M2, M3) +
  eq(e_(I1, I2, I3)) *
    1/36 * (d_(I1, I2) * (d_(P1, I3) + d_(P2, I3)) +
            d_(I1, I3) * (d_(P1, I2) + d_(P2, I2)) +
            d_(I2, I3) * (d_(P1, I1) + d_(P2, I1)));

id c0(I1?, I1?, P1?, P2?, M1?, M2?, M3?) =
  c0(I1, I1, P1, P2, M1, M2, M3) - 1/2;

id eq(0) = 1;
id eq(?) = 0;

id b0(I1?, I1?, P1?, M1?, M2?) = a0(M2) + M1 * b0(P1, M1, M2);
#endif


id d0(I1?, I2?, I3?, I4?, P1?, P2?, P3?, M1?, M2?, M3?, M4?) =
  pave4(0, 0, 0, 0, P1, P2, P3, M1, M2, M3, M4) *
     (d_(I1, I2) * d_(I3, I4) +
      d_(I1, I3) * d_(I2, I4) +
      d_(I1, I4) * d_(I2, I3))
#do i = 1, 3
#do j = 1, 3
  + pave4(0, 0, 'i', 'j', P1, P2, P3, M1, M2, M3, M4) *
     (d_(I1, I2) * d_(P'i', I3) * d_(P'j', I4) +
      d_(I2, I3) * d_(P'i', I1) * d_(P'j', I4) +
      d_(I1, I3) * d_(P'i', I2) * d_(P'j', I4) +
      d_(I1, I4) * d_(P'i', I2) * d_(P'j', I3) +
      d_(I2, I4) * d_(P'i', I1) * d_(P'j', I3) +
      d_(I3, I4) * d_(P'i', I1) * d_(P'j', I2))
#do k = 1, 3
#do l = 1, 3
  + pave4('i', 'j', 'k', 'l', P1, P2, P3, M1, M2, M3, M4) *
      d_(P'i', I1) * d_(P'j', I2) * d_(P'k', I3) * d_(P'l', I4)
#enddo
#enddo
#enddo
#enddo
  ;

id d0(I1?, I2?, I3?, P1?, P2?, P3?, M1?, M2?, M3?, M4?) =
#do i = 1, 3
  + pave4(0, 0, 'i', P1, P2, P3, M1, M2, M3, M4) *
     (d_(I1, I2) * d_(P'i', I3) +
      d_(I2, I3) * d_(P'i', I1) +
      d_(I1, I3) * d_(P'i', I2))
#do j = 1, 3
#do k = 1, 3
  + pave4('i', 'j', 'k', P1, P2, P3, M1, M2, M3, M4) *
      d_(P'i', I1) * d_(P'j', I2) * d_(P'k', I3)
#enddo
#enddo
#enddo
  ;

id d0(I1?, I2?, P1?, P2?, P3?, M1?, M2?, M3?, M4?) =
  pave4(0, 0, P1, P2, P3, M1, M2, M3, M4) * d_(I1, I2)
#do i = 1, 3
#do j = 1, 3
  + pave4('i', 'j', P1, P2, P3, M1, M2, M3, M4) *
      d_(P'i', I1) * d_(P'j', I2)
#enddo
#enddo
  ;

id d0(I1?, P1?, P2?, P3?, M1?, M2?, M3?, M4?) =
#do i = 1, 3
  + pave4('i', P1, P2, P3, M1, M2, M3, M4) * d_(P'i', I1)
#enddo
  ;

id c0(I1?, I2?, I3?, P1?, P2?, M1?, M2?, M3?) =
#do i = 1, 2
  + pave3(0, 0, 'i', P1, P2, M1, M2, M3) *
     (d_(I1, I2) * d_(P'i', I3) +
      d_(I2, I3) * d_(P'i', I1) +
      d_(I1, I3) * d_(P'i', I2))
#do j = 1, 2
#do k = 1, 2
  + pave3('i', 'j', 'k', P1, P2, M1, M2, M3) *
      d_(P'i', I1) * d_(P'j', I2) * d_(P'k', I3)
#enddo
#enddo
#enddo
  ;

id c0(I1?, I2?, P1?, P2?, M1?, M2?, M3?) =
  pave3(0, 0, P1, P2, M1, M2, M3) * d_(I1, I2)
#do i = 1, 2
#do j = 1, 2
  + pave3('i', 'j', P1, P2, M1, M2, M3) * d_(P'i', I1) * d_(P'j', I2)
#enddo
#enddo
  ;

id c0(I1?, P1?, P2?, M1?, M2?, M3?) =
#do i = 1, 2
  + pave3('i', P1, P2, M1, M2, M3) * d_(P'i', I1)
#enddo
  ;

id b0(I1?, I2?, P1?, M1?, M2?) =
  B00m(P1, M1, M2) * d_(I1, I2) +
  B11m(P1, M1, M2) * d_(P1, I1) * d_(P1, I2);

id b0(I1?, P1?, M1?, M2?) = B1m(P1, M1, M2) * d_(P1, I1);

symm pave4:11 1, 2, 3, 4;
symm pave4:10 1, 2, 3;
symm pave4:9 1, 2;
symm pave3:8 1, 2, 3;
symm pave3:7 1, 2;

id a0(0) = 0;
id a0(M1?) = A0(M1);
id b0(P1?, M1?, M2?) = B0m(P1, M1, M2);
symm B0m 2, 3;
id c0(P1?, P2?, M1?, M2?, M3?) = pave3(0, P1, P2, M1, M2, M3);
id d0(P1?, P2?, P3?, M1?, M2?, M3?, M4?) =
  pave4(0, P1, P2, P3, M1, M2, M3, M4);

id a0(?) = unknown;
id b0(?) = unknown;
id c0(?) = unknown;
id d0(?) = unknown;

.sort


* >>> Dirac algebra for open fermion chains:

if( count(ga, 1, ga5, 1, Spinor, 1) > 0 );
  repeat;
#call diracsimplify{}
#if 'DIM' == 4
* apply the Chisholm identity
    repeat;
      id,once, ga(C?, I1?) * ga(C?, I2?) * ga(C?, I3?) =
        ga5(C) * ga(C, I4) * e_(I1, I2, I3, I4) +
        d_(I1, I2) * ga(C, I3) - d_(I1, I3) * ga(C, I2) +
        d_(I2, I3) * ga(C, I1);
      id ga5(C?) * ga5(C?) = 1;
      sum I4;
    endrepeat;
    contract,0;
#endif
  endrepeat;
  id pivot(?) = 1;
  id ga(C?, I1?) = fc(C, ga(I1));
  id ga5(C?) = fc(C, ga5);
  repeat;
    id fc(C?, ?) * fc(C?, ??) = fc(C, ., ..);
  endrepeat;
  id Spinor(P1?, M1?) * fc(C?, ?) * Spinor(P2?, M2?) =
    fme(Spinor(P1, M1), ., Spinor(P2, M2));
  id Spinor(P1?, M1?) * Spinor(P2?, M2?) =
    fme(Spinor(P1, M1), Spinor(P2, M2));
  id fc(C?, ?) = fme(.);
endif;

.sort


* >>> cleaning up:

#call momcons{}

#if 'DIM' == 0
repeat;
  id,once, D * A0(M1?) = 4 * A0(M1) - 2 * M1;
  id,once, D * B0m(?) = 4 * B0m(.) - 2;
  id,once, D * B11m(?) = 4 * B11m(.) - 2/3;
  id,once, D * B1m(?) = 4 * B1m(.) + 1;
  id,once, D * B00m(P1?, M1?, M2?) =
    4 * B00m(P1, M1, M2) + 1/6 * P1.P1 - 1/2 * M1 - 1/2 * M2;
  id,once, D * pave3(0, 0, P1?, P2?, M1?, M2?, M3?) =
    4 * pave3(0, 0, P1, P2, M1, M2, M3) - 1/2;
#if 'VERSION_' > 1
  id,once, D * pave3(0, 0, I1?, P1?, P2?, M1?, M2?, M3?) =
    4 * pave3(0, 0, I1, P1, P2, M1, M2, M3) + 1/6;
#else
  id,once, D * pave3(I1?, 0, 0, P1?, P2?, M1?, M2?, M3?) =
    4 * pave3(I1, 0, 0, P1, P2, M1, M2, M3) + 1/6;
#endif
  id,once, D * pave4(0, 0, 0, 0, ?) = 4 * pave4(0, 0, 0, 0, .) - 1/12;
endrepeat;
id D = 4;
#endif

#call mandelstam{}

b A0, B0m, B1m, B00m, B11m, pave3, pave4, DEN, i_;

print +s amp;

