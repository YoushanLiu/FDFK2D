EXEC = ~/bin/gauss_smooth

OBJ_DIR = .
SRC_DIR = ./


FC = ifort
#FC = gfortran


ifeq ($(FC),ifort)
	FFLAGS = -O3 -fpp -qopenmp -xHost -unroll -ftz -fma -fpp -fpe0 -assume buffered_io -assume byterecl -align sequence -diag-disable 6477 -implicitnone -gen-interfaces -mcmodel=large -shared-intel -arch AVX
else ifeq ($(FC),gfortran)
	FFLAGS = -O3 -Ofast -cpp -fopenmp -g -funroll-loops -malign-double -ffast-math -mtune=native
else
	FFLAGS = -O3 -fast -mp -m64 -Mpreprocess -Mnobounds -Mautoinline -Mcache_align -Mfma -Mprefetch -Munroll -Mdaz -mcmodel=medium -Mmpi=mpich2 -flushz -Mmavx
endif

LIBS =


OBJEXT = o
LINKOPT = 

UNAME := $(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
   OBJ = obj
   LINKOPT = /F150994944
   FFLAGS = /nologo /QxHost /fpp /O3 /Qparallel /Qopt-matmul /Qopenmp /Qmkl:parallel
endif


ifeq ($(findstring gfortran, $(FC)), gfortran)
	# gfortran
	FFLAGS += -J$(OBJ_DIR)
else
	# ifort or pgf90
	FFLAGS += -module $(OBJ_DIR)
endif

FFLAGS += -I$(OBJ_DIR)



OBJS = \
	gauss_smooth.$(OBJEXT) \
	$(EMPTY_MACRO)


OBJS := $(addprefix $(OBJ_DIR)/,$(OBJS))



all: $(OBJS)
	mkdir -p ~/bin
	-rm -rf $(EXEC)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -rf $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*_genmod.f90

cleanall: clean
	rm -rf $(INST_DIR)/$(BIN)


.PHONY: all


$(OBJ_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/%.o: %.for
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

