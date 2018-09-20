#Files
EXEC := forward.out
SRC := $(wildcard *.c)
OBJ := $(patsubst %.c,%.o,$(SRC))
#Options
CC := icc
CFLAGS := -openmp -O2
TACC_GRVY_DIR := -L$(TACC_GRVY_LIB)
TACC_HDF5_DIR := -L$(TACC_HDF5_LIB)
LDLIBS := -lgrvy -openmp 
#-lhdf5 -lz
TACC_GRVY_COMP := -I$(TACC_GRVY_INC)
TACC_HDF5_COMP := -I$(TACC_HDF5_INC) 
#Rules
$(EXEC) : $(OBJ)
	$(CC) $(TACC_GRVY_DIR)  $(LDLIBS)  -o $@ $^
%.o : %.c
	$(CC) $(CFLAGS) $(TACC_GRVY_COMP)  -c $<
clean:
	$(RM) $(EXEC)
	$(RM) $(OBJ) 
