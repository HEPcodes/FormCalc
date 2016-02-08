OBJS_DIAG := \
  $(DIR)Eigen2x2.o \
  $(DIR)NoMixing.o \
  $(DIR)HEigensystem.o \
  $(DIR)SEigensystem.o \
  $(DIR)CEigensystem.o \
  $(DIR)TakagiFactor.o \
  $(DIR)SVD.o

OBJS += $(OBJS_DIAG)

$(UTILLIB)($(OBJS_DIAG)): $(DIR)diag-f.h $(SRC)/F/types.h

