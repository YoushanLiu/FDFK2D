# Makefile for Teleseismic Wavefield modeling to be release

EXEC = ~/bin/FDFK2D

OBJ_DIR = ../obj
SRC_DIR = ./


FC = ifort
FC = gfortran
FC = pgfortran


ifeq ($(FC),ifort)
	FFLAGS = -O3 -qopenmp -xHost -unroll -ftz -fma -fpp -fpe0 -assume buffered_io -assume byterecl -align sequence -diag-disable 6477 -implicitnone -gen-interfaces -mcmodel=large -shared-intel -arch AVX
	LIBS = -lfftw3 -llapack -lblas -L../lib_intel -lHybrid -lPML -lFK3D
else ifeq ($(FC),gfortran)
	FFLAGS = -O3 -Ofast -fopenmp -g -funroll-loops -malign-double -ffast-math -mtune=native
	LIBS = -lfftw3 -llapack -lblas -L../lib_gnu -lHybrid -lPML -lFK3D
else
	FFLAGS = -O3 -fast -mp -m64 -Mpreprocess -Mnobounds -Mautoinline -Mcache_align -Mfma -Mprefetch -Munroll -Mdaz -mcmodel=medium -Mflushz
	LIBS = -lfftw3 -llapack -lblas -L../lib_pgi -lHybrid -lPML -lFK3D
endif

#LIBS = -lfftw3 -llapack -lblas -lHybrid -lPML -lFK3D
#LIBS = -lfftw3 -llapack -lblas -L../release_lib/ libHybrid.a libPML.a libFK3D.a
#LIBS = -lfftw3 -llapack -lblas -L../release_lib -lHybrid -lPML -lFK3D


OBJEXT = o
LINKOPT = 

UNAME := $(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
   OBJ = obj
   LINKOPT = /F150994944
   FFLAGS = /nologo /QxHost /fpp /O3 /Qparallel /Qopt-matmul /Qopenmp /Qmkl:parallel
endif


ifeq ($(FC),gfortran)
	# gfortran
	FFLAGS += -J$(OBJ_DIR)
else
	# ifort or pgf90
	FFLAGS += -module $(OBJ_DIR)
endif


FFLAGS += -I$(OBJ_DIR)



OBJS = \
	constants.$(OBJEXT) \
	calc_backazimuth.$(OBJEXT)  \
	define_rotation_matrix.$(OBJEXT) \
	coord_transform.$(OBJEXT) \
	FDFK_par.$(OBJEXT) \
	Input.$(OBJEXT) \
	Submain.$(OBJEXT) \
	utils.$(OBJEXT) \
	FDFK2D.$(OBJEXT)\
	$(EMPTY_MACRO)


OBJS := $(addprefix $(OBJ_DIR)/,$(OBJS))



all: build


lib:
	@if [ "$(FC)" = "ifort" ]; then \
		echo "You are using $(FC) compiler"; \
		cp ../lib_intel/*.so /usr/lib64; \
	elif [ "$(FC)" = "gfortran" ]; then \
		echo "You are using $(FC) compiler"; \
		cp ../lib_gnu/*.so /usr/lib64; \
	else \
		echo "You are using $(FC) compiler"; \
		cp ../lib_pgi/*.so /usr/lib64; \
	fi


build: $(OBJS)
	-rm -rf $(EXEC)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -rf $(OBJ_DIR)/* *.mod *_genmod.f90

cleanall: clean
	rm -rf $(INST_DIR)/$(BIN)


$(OBJ_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/%.o: %.for
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@
	
	
.PHONY: all



# dependencies
$(OBJ_DIR)/calc_backazimuth.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT)
$(OBJ_DIR)/define_rotation_matrix.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT)
$(OBJ_DIR)/coord_transform.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT) $(OBJ_DIR)/calc_backazimuth.$(OBJEXT) $(OBJ_DIR)/FDFK_par.$(OBJEXT)
$(OBJ_DIR)/FDFK2D.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT) $(OBJ_DIR)/FDFK_par.$(OBJEXT)
$(OBJ_DIR)/Input.$(OBJEXT): $(OBJ_DIR)/FDFK_par.$(OBJEXT)

