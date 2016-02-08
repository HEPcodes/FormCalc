OBJS_NINJA := \
  $(DIR)NJBcut.o \
  $(DIR)NJCcut.o \
  $(DIR)NJDcut.o \
  $(DIR)NJEcut.o \
  $(DIR)NJFcut.o

OBJS += $(DIR)NJMomDecode.o $(DIR)None.o $(OBJS_NINJA)

$(UTILLIB)($(OBJS_NINJA)): $(DIR)NJ.h $(SRC)/F/types.h

$(UTILLIB)($(DIR)NJMomDecode.o): $(SRC)/F/types.h $(SRC)/F/util.h

