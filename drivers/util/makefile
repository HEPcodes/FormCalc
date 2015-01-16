parent := $(DIR)/..
system := $(DIR)/system
kin := $(DIR)/kin
cuttools := $(DIR)/cuttools
samurai := $(DIR)/samurai
diag := $(DIR)/diag
univariate := $(DIR)/univariate
compaz := $(DIR)/compaz
cuba := $(DIR)/cuba

common := $(cuba)/common
vegas := $(cuba)/vegas
suave := $(cuba)/suave
divonne := $(cuba)/divonne
cuhre := $(cuba)/cuhre

CFLAGS += -I$(common) -I$(@D) -DHAVE_CONFIG_H

COMMON_DEPS := $(common)/stddecl.h \
  $(common)/ChiSquare.c $(common)/Erf.c \
  $(common)/CSample.c $(common)/Fork.c

VEGAS_C := $(vegas)/Vegas.c
VEGAS_DEPS := $(COMMON_DEPS) $(common)/Random.c \
  $(vegas)/decl.h $(vegas)/common.c \
  $(vegas)/Grid.c $(vegas)/Integrate.c

SUAVE_C := $(suave)/Suave.c
SUAVE_DEPS := $(COMMON_DEPS) $(common)/Random.c \
  $(suave)/decl.h $(suave)/common.c \
  $(suave)/Grid.c $(suave)/Fluct.c $(suave)/Sample.c \
  $(suave)/Integrate.c

DIVONNE_C := $(divonne)/Divonne.c
DIVONNE_DEPS := $(COMMON_DEPS) $(common)/Random.c \
  $(divonne)/decl.h $(divonne)/common.c \
  $(divonne)/KorobovCoeff.c $(divonne)/Rule.c $(divonne)/Sample.c \
  $(divonne)/FindMinimum.c $(divonne)/Iterate.c $(divonne)/Explore.c \
  $(divonne)/Split.c $(divonne)/Integrate.c

CUHRE_C := $(cuhre)/Cuhre.c
CUHRE_DEPS := $(COMMON_DEPS) \
  $(cuhre)/decl.h $(cuhre)/common.c \
  $(cuhre)/Rule.c $(cuhre)/Integrate.c

CUBA := $(VEGAS_C) $(SUAVE_C) $(DIVONNE_C) $(CUHRE_C) \
  $(common)/WorkerIni.c


OBJS := \
  $(CUBA:.c=.o) \
  $(system)/logfile.o \
  $(system)/setfpu.o \
  $(system)/Error.o \
  $(system)/fortranflush.o \
  $(system)/parallel.o \
  $(system)/shmalloc.o \
  $(kin)/VecSet.o \
  $(kin)/VecCopy.o \
  $(kin)/Split.o \
  $(kin)/PrintVec.o \
  $(cuttools)/CTMomDecode.o \
  $(cuttools)/CTmpnum.o \
  $(cuttools)/CTBcut.o \
  $(cuttools)/CTCcut.o \
  $(cuttools)/CTDcut.o \
  $(cuttools)/CTEcut.o \
  $(cuttools)/CTFcut.o \
  $(samurai)/SAMomDecode.o \
  $(samurai)/SABcut.o \
  $(samurai)/SACcut.o \
  $(samurai)/SADcut.o \
  $(samurai)/SAEcut.o \
  $(samurai)/SAFcut.o \
  $(diag)/Eigen2x2.o \
  $(diag)/NoMixing.o \
  $(diag)/HEigensystem.o \
  $(diag)/SEigensystem.o \
  $(diag)/CEigensystem.o \
  $(diag)/TakagiFactor.o \
  $(diag)/SVD.o \
  $(univariate)/Patterson.o \
  $(univariate)/Gauss.o \
  $(compaz)/CompAZ.o


$(LIB): $(LIB)($(OBJS))


$(LIB)( $(univariate)/Gauss.o $(univariate)/Patterson.o ): \
  $(univariate)/config.h $(parent)/F/types.h

$(LIB)( $(diag)/HEigensystem.o $(diag)/SEigensystem.o \
  $(diag)/CEigensystem.o $(diag)/TakagiFactor.o $(diag)/SVD.o ): \
  $(diag)/diag-f.h $(parent)/F/types.h

$(LIB)( $(samurai)/SABcut.o $(samurai)/SACcut.o $(samurai)/SADcut.o \
  $(samurai)/SAEcut.o $(samurai)/SAFcut.o ): \
  $(samurai)/SA.h $(parent)/F/types.h

$(LIB)( $(cuttools)/CTBcut.o $(cuttools)/CTCcut.o $(cuttools)/CTDcut.o \
  $(cuttools)/CTEcut.o $(cuttools)/CTFcut.o ): \
  $(cuttools)/CT.h $(parent)/F/types.h

$(LIB)( $(samurai)/SAMomDecode.o $(cuttools)/CTMomDecode.o ): \
  $(parent)/F/types.h $(parent)/F/util.h

$(LIB)( $(kin)/VecCopy.o ): \
  $(parent)/distrib.h simd.h $(parent)/F/types.h $(parent)/F/util.h

$(LIB)( $(system)/parallel.o ): $(parent)/distrib.h

$(LIB)($(VEGAS_C:.c=.o)): $(VEGAS_C) $(VEGAS_DEPS)
$(LIB)($(SUAVE_C:.c=.o)): $(SUAVE_C) $(SUAVE_DEPS)
$(LIB)($(DIVONNE_C:.c=.o)): $(DIVONNE_C) $(DIVONNE_DEPS)
$(LIB)($(CUHRE_C:.c=.o)): $(CUHRE_C) $(CUHRE_DEPS)
$(LIB)($(common)/WorkerIni.o): $(common)/stddecl.h

LIBS += $(LIB)

