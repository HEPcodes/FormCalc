OBJS_CUTTOOLS := \
  $(DIR)CTBcut.o \
  $(DIR)CTCcut.o \
  $(DIR)CTDcut.o \
  $(DIR)CTEcut.o \
  $(DIR)CTFcut.o

OBJS += $(DIR)CTMomDecode.o $(DIR)CTmpnum.o $(OBJS_CUTTOOLS)

$(UTILLIB)($(OBJS_CUTTOOLS)): $(DIR)CT.h $(SRC)/F/types.h

$(UTILLIB)($(DIR)CTMomDecode.o): $(SRC)/F/types.h $(SRC)/F/util.h

