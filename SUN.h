* SUN.h
* simplification of SU(N) structures
* this file is part of FormCalc
* last modified 26 Apr 00 th

* The algorithm implemented here is an extension of the one given in
*   J.A.M. Vermaseren, The use of computer algebra in QCD,
*   in: H. Latal, W. Schweiger, Proceedings Schladming 1996,
*       Springer Verlag, ISBN 3-540-62478-3.

* The idea is to transform all SU(N) objects to generators, SUNT.
* In the output, only two types of objects can appear:
* - chains of SUNTs (with external colour indices), or
* - traces of SUNTs.
* A chain of SUNTs is denoted by SUNT(a, b, ..., i, j), where
* a, b, ... are gluon indices and i and j are colour indices.
* SUNT(i, j) is the special case of the identity in colour space.
* A trace over SUNTs is marked by both colour indices being zero,
* i.e. SUNT(a, b, ..., 0, 0).


#if 'VERSION_' > 1
b SUNT, SUNTSum, SUNF, SumOver;
.sort
keep brackets;
#else
.sort
#endif

d 'SUNN';
i COL1, COL2, COL3, COL4;
i GLU1, GLU2, GLU3, GLU4, GLU5;
s X1, X2;
cf bead;


id SUNT(?) = SUNT(.) * sun(.);
id SUNTSum(?) = SUNTSum(.) * sun(.);
id SUNF(?) = SUNF(.) * sun(.);
repeat;
  id sun(COL1?, ?) * SumOver(COL1?, COL2?, X1?) = sun(COL1, .);
  id sun(COL1?, ?) = sun(.);
endrepeat;
id sun(?) = 1;


if(count(SUNF, 1) > 0);

  repeat;
    id,once, SUNF(?, GLU1?, GLU2?, GLU3?, GLU4?) =
      SUNF(., GLU1, GLU2, GLU5) * SUNF(GLU5, GLU3, GLU4) * SumOver(GLU5);
    sum GLU5;
  endrepeat;

* f^{abc} = 2 i Tr(T^c T^b T^a - T^a T^b T^c)

  id SUNF(GLU1?, GLU2?, GLU3?) =
    2*i_*(SUNT(GLU3, GLU2, GLU1, 0, 0) - SUNT(GLU1, GLU2, GLU3, 0, 0));

endif;


repeat;
  id,once, SUNT(?, 0, 0) = SUNT(., COL1, COL1) * SumOver(COL1);
  sum COL1;
endrepeat;

repeat;
  id,once, SUNT(?, GLU1?, GLU2?, COL1?, COL2?) =
    SUNT(., GLU1, COL1, COL3) * SUNT(GLU2, COL3, COL2) * SumOver(COL3);
  sum COL3;
endrepeat;


* T^a_{ij} T^a_{kl} =
*   1/2 (delta_{il} delta_{jk} - 1/N delta_{ij} delta_{kl})

id SUNT(GLU1?, COL1?, COL2?) * SUNT(GLU1?, COL3?, COL4?) *
    SumOver(GLU1?, ?) =
  1/2 * SUNT(COL1, COL4) * SUNT(COL2, COL3) -
  1/2/'SUNN' * SUNT(COL1, COL2) * SUNT(COL3, COL4);

id SUNTSum(COL1?, COL2?, COL3?, COL4?) =
  1/2 * SUNT(COL1, COL4) * SUNT(COL2, COL3) -
  1/2/'SUNN' * SUNT(COL1, COL2) * SUNT(COL3, COL4);


* cleaning up, step 1: get rid of the deltas

repeat;
  id,once, SUNT(COL1?, COL2?) * SumOver(COL1?, ?) = d_(COL1, COL2);
  id,once, SUNT(COL2?, COL1?) * SumOver(COL1?, ?) = d_(COL1, COL2);
endrepeat;

id SUNT(COL1?, COL1?) = 1;
id SUNT(GLU1?, COL1?, COL1?) * SumOver(COL1?, ?) = 0;

symm SUNT:2, 1, 2;

* cleaning up, step 2: bead up the SUNTs into chains

repeat;
  id,once, SUNT(?, GLU1?, COL1?, COL2?) = bead(., GLU1, COL1, COL2);
  repeat;
    id bead(?, COL1?, COL2?) * SUNT(??, COL2?, COL3?) *
        SumOver(COL2?, ???) =
      bead(., .., COL1, COL3);
    id SUNT(?, COL1?, COL2?) * bead(??, COL2?, COL3?) *
        SumOver(COL2?, ???) =
      bead(., .., COL1, COL3);
  endrepeat;
  id bead(?, COL1?, COL1?) * SumOver(COL1?, ??) = bead(., 0, 0);
* special case of Tr(T^a T^b) = 1/2 delta_{ab}
  id bead(GLU1?, GLU1?, 0, 0) = 1/2;
  id bead(?) = sun(SUNT(.));
endrepeat;

id SUNT(?) = sun(SUNT(.));

repeat;
  id sun(X1?) * sun(X2?) = sun(X1 * X2);
endrepeat;

.sort

id sun(X1?) = Mat(sun(X1));
repeat;
  id Mat(X1?) * Mat(X2?) = Mat(X1 * X2);
endrepeat;

