# Makefile for Teleseismic Wavefield modeling to be release

EXEC = ~/bin/FDFK2D

OBJ_DIR = ../obj
SRC_DIR = ./


FC = ifort
FFLAGS = -O3 -xHost -fpp -fpe0 -ftz -mcmodel=medium -qopt-matmul -qopenmp #-mkl=parallel
LIBS = -lfftw3 -llapack -lblas -lhybrid -lPML -lFK3D


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



all: lib build install


lib:
	cp ../lib/*.so /usr/lib64

build:
$(OBJ_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/%.o: %.for
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJ_DIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

install: $(OBJS)
	mkdir -p ~/bin
	-rm -rf $(EXEC)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -rf $(OBJ_DIR)/* *.mod *_genmod.f90

cleanall: clean
	rm -rf $(INST_DIR)/$(BIN)


.PHONY: all





# dependencies
$(OBJ_DIR)/calc_backazimuth.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT)
$(OBJ_DIR)/define_rotation_matrix.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT)
$(OBJ_DIR)/coord_transform.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT) $(OBJ_DIR)/calc_backazimuth.$(OBJEXT) $(OBJ_DIR)/FDFK_par.$(OBJEXT)
$(OBJ_DIR)/FDFK2D.$(OBJEXT): $(OBJ_DIR)/constants.$(OBJEXT) $(OBJ_DIR)/FDFK_par.$(OBJEXT)
$(OBJ_DIR)/Input.$(OBJEXT): $(OBJ_DIR)/FDFK_par.$(OBJEXT)

