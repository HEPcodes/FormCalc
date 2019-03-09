common := $(DIR)common
vegas := $(DIR)vegas
suave := $(DIR)suave
divonne := $(DIR)divonne
cuhre := $(DIR)cuhre

COMMON_DEPS := $(common)/stddecl.h \
  $(common)/ChiSquare.c $(common)/Erf.c \
  $(common)/CSample.c $(common)/Parallel.c $(common)/sock.h

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
  $(common)/Fork.c $(common)/Global.c $(common)/Data.c
CUBA_OBJS = $(CUBA:.c=.o)

OBJS += $(CUBA_OBJS)

$(UTILLIB)($(CUBA_OBJS)): CFLAGS += -I$(common) -I$(@D) -DHAVE_CONFIG_H

$(UTILLIB)($(VEGAS_C:.c=.o)): $(VEGAS_C) $(VEGAS_DEPS) 
$(UTILLIB)($(SUAVE_C:.c=.o)): $(SUAVE_C) $(SUAVE_DEPS) 
$(UTILLIB)($(DIVONNE_C:.c=.o)): $(DIVONNE_C) $(DIVONNE_DEPS) 
$(UTILLIB)($(CUHRE_C:.c=.o)): $(CUHRE_C) $(CUHRE_DEPS) 
$(UTILLIB)($(common)/Fork.o): $(common)/Fork.c $(common)/stddecl.h $(common)/sock.h
$(UTILLIB)($(common)/Global.o): $(common)/Global.c $(common)/stddecl.h
$(UTILLIB)($(common)/Data.o): $(common)/Data.c $(common)/stddecl.h

