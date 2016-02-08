OBJS_KIN := \
  $(DIR)VecSet.o \
  $(DIR)VecCopy.o \
  $(DIR)Split.o \
  $(DIR)PrintVec.o

OBJS += $(OBJS_KIN)

$(UTILLIB)($(OBJS_KIN)): $(SRC)/F/types.h

$(UTILLIB)($(DIR)VecCopy.o $(DIR)VecSet.o): \
  $(SRC)/distrib.h simd.h $(SRC)/F/util.h

$(UTILLIB)($(DIR)Split.o): $(SRC)/F/util.h

