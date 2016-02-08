OBJS_SAMURAI := \
  $(DIR)SABcut.o \
  $(DIR)SACcut.o \
  $(DIR)SADcut.o \
  $(DIR)SAEcut.o \
  $(DIR)SAFcut.o

OBJS += $(DIR)SAMomDecode.o $(OBJS_SAMURAI)

$(UTILLIB)($(OBJS_SAMURAI)): $(DIR)SA.h $(SRC)/F/types.h

$(UTILLIB)($(DIR)SAMomDecode.o): $(SRC)/F/types.h $(SRC)/F/util.h

