# ================= 3rd Party Libraries ====================================

OPENMESH_DIR = /Users/karrywong/Desktop/research/OpenMesh-7.1/
OPENMESH_INCLUDE_DIR = $(OPENMESH_DIR)/src/
OPENMESH_LIB_DIR = $(OPENMESH_DIR)/buildDebug/Build/lib/
OPENMESH_LIBS = -lOpenMeshCore -lOpenMeshTools

EIGEN_INCLUDE_DIR = /Users/karrywong/Desktop/research/eigen3

#extra directory paths if needed
##SUITESPARSE_DIR = /homes/home03/w/wong0925/openmesh/SuiteSparse
##SUITESPARSE_INCLUDE_DIR = $(SUITESPARSE_DIR)/include/
##SUITESPARSE_LIB_DIR = $(SUITESPARSE_DIR)/lib/
##SUITESPARSE_LIBS = -lcholmod -lamd -lcamd -lcolamd -lccolamd  -lsuitesparseconfig -lmetis -lopenblas -lumfpack

##HOME_DIR = /homes/home03/w/wong0925/openmesh/OpenMatch
##HOME_LIB_DIR = $(HOME_DIR)/lib/
#HOME_LIBS = -lgomp

# ================= Project Directories ====================================

BIN_DIR = ../bin

# ================= Project Name ===========================================

EXT=.exe
NAME=Map2Sphere
NAMEFUL = $(BIN_DIR)/$(NAME)$(EXT)

# ================= Compilers and Flags ====================================

CC = gcc 
CFLAGS = -c -O
CPP = g++
CPPFLAGS = -c -O -std=c++11
INCLUDE_DIRS = -I$(OPENMESH_INCLUDE_DIR) -I$(EIGEN_INCLUDE_DIR)
LIB_DIRS = -L$(OPENMESH_LIB_DIR) 
LIBS = $(OPENMESH_LIBS) 
LOAD_LIB_PATH = -Wl,-rpath, $(OPENMESH_LIB_DIR) 
LD_FLAGS = -O

.cpp.o :
	$(CPP) $(CPPFLAGS) $(INCLUDE_DIRS) $<

.c.o :
	$(C) $(CFLAGS) $(INCLUDE_DIRS) $<

.h.o :
	$(H) $(HFLAGS) $(INCLUDE_DIRS) $<

# ================= Compile source code ====================================

OBJECTS = \
$(NAME).o \
CheckMesh.o\

# ================= Generate Executable ====================================

$(NAMEFUL) : $(OBJECTS)
	$(CPP) -o $(NAMEFUL) $(LD_FLAGS) $(OBJECTS) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

all: $(OBJECTS)
	$(CPP) -o $(NAMEFUL) $(LD_FLAGS) $(OBJECTS) $(LIB_DIRS) $(LIBS) $(LOAD_LIB_PATH)

clean:
	touch junk.o; rm -f *.o $(NAMEFUL)
