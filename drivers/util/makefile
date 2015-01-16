LIB = $(CURDIR).a

common = cuba/common
vegas = cuba/vegas
suave = cuba/suave
divonne = cuba/divonne
cuhre = cuba/cuhre

CFLAGS += -I$(common)

VEGAS_C = $(vegas)/Vegas.c
VEGAS_F = $(vegas)/vegas-f.c
VEGAS_H = $(vegas)/decl.h $(common)/stddecl.h
VEGAS_SRCS = $(common)/debug.c $(common)/Random.c \
  $(common)/ChiSquare.c $(common)/Erf.c $(vegas)/util.c \
  $(vegas)/common.c $(vegas)/Grid.c $(vegas)/Integrate.c

SUAVE_C = $(suave)/Suave.c
SUAVE_F = $(suave)/suave-f.c
SUAVE_H = $(suave)/decl.h $(common)/stddecl.h
SUAVE_SRCS = $(common)/debug.c $(common)/Random.c \
  $(common)/ChiSquare.c $(common)/Erf.c $(suave)/util.c \
  $(suave)/common.c $(suave)/Grid.c $(suave)/Fluct.c \
  $(suave)/Sample.c $(suave)/Integrate.c

DIVONNE_C = $(divonne)/Divonne.c
DIVONNE_F = $(divonne)/divonne-f.c
DIVONNE_H = $(divonne)/decl.h $(common)/stddecl.h
DIVONNE_SRCS = $(common)/debug.c $(common)/Random.c \
  $(common)/ChiSquare.c $(common)/Erf.c $(divonne)/util.c \
  $(divonne)/common.c $(divonne)/KorobovCoeff.c \
  $(divonne)/Rule.c $(divonne)/Sample.c $(divonne)/FindMinimum.c \
  $(divonne)/Explore.c $(divonne)/Split.c $(divonne)/Integrate.c

CUHRE_C = $(cuhre)/Cuhre.c
CUHRE_F = $(cuhre)/cuhre-f.c
CUHRE_H = $(cuhre)/decl.h $(common)/stddecl.h
CUHRE_SRCS = $(common)/debug.c \
  $(common)/ChiSquare.c $(common)/Erf.c $(cuhre)/util.c \
  $(cuhre)/common.c $(cuhre)/Rule.c $(cuhre)/Integrate.c

CUBA = \
  $(VEGAS_C) $(VEGAS_F) \
  $(SUAVE_C) $(SUAVE_F) \
  $(DIVONNE_C) $(DIVONNE_F) \
  $(CUHRE_C) $(CUHRE_F)


OBJS = \
  system/logfile.o \
  kin/ThreeMom.o \
  kin/VecSet.o \
  kin/Split.o \
  kin/SInvariant.o \
  kin/TInvariant.o \
  kin/Pair.o \
  kin/Eps.o \
  kin/SxS.o \
  kin/SeS.o \
  kin/VxS.o \
  kin/VeS.o \
  kin/BxS.o \
  kin/BeS.o \
  diag/Jacobi.o \
  diag/Eigen2x2.o \
  diag/NoMixing.o \
  diag/Eigensystem.o \
  diag/SingularValues.o \
  diag/PseudoEigensystem.o \
  diag/ApplyEigen.o \
  diag/ApplySVD.o \
  diag/ReverseEigen.o \
  diag/ReverseSVD.o \
  diag/IdentityMatrix.o \
  univariate/Patterson.o \
  univariate/Gauss.o \
  $(CUBA:.c=.o) \
  compaz/CompAZ.o


$(LIB): $(LIB)($(OBJS))


$(LIB)(univariate/Gauss.o univariate/Patterson.o): univariate/config.h

$(LIB)(diag/Eigensystem.o diag/SingularValues.o diag/PseudoEigensystem.o \
  diag/ApplyEigen.o diag/ApplySVD.o): diag/config.h


$(LIB)($(VEGAS_C:.c=.o)): $(VEGAS_C) $(VEGAS_H) $(VEGAS_SRCS)
$(LIB)($(VEGAS_F:.c=.o)): $(VEGAS_F) $(VEGAS_H)

$(LIB)($(SUAVE_C:.c=.o)): $(SUAVE_C) $(SUAVE_H) $(SUAVE_SRCS)
$(LIB)($(SUAVE_F:.c=.o)): $(SUAVE_F) $(SUAVE_H)

$(LIB)($(DIVONNE_C:.c=.o)): $(DIVONNE_C) $(DIVONNE_H) $(DIVONNE_SRCS)
$(LIB)($(DIVONNE_F:.c=.o)): $(DIVONNE_F) $(DIVONNE_H)

$(LIB)($(CUHRE_C:.c=.o)): $(CUHRE_C) $(CUHRE_H) $(CUHRE_SRCS)
$(LIB)($(CUHRE_F:.c=.o)): $(CUHRE_F) $(CUHRE_H)


clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) $(LIB)

