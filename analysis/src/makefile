#
# C++ compiler and options
#
CXX = g++
CXXFLAGS = -O3 -Wall -Wshadow -std=c++11
LIB = -lm

#
# Directories for source and object files and executables
#
SDIR = ./
LDIR = ../lammps
ODIR = ../obj
OBDIR = $(ODIR)/bin
OLDIR = $(ODIR)/lib
EDIR = ./

#
# Header files
#
#SRCH_NAME

#
# Source files
#
SRC_NAME = Get_Bead_Position_LAMMPS.cpp
#SRC_NAME = ReadLAMMPS.cpp

#
# Library (LAMMPS) files
LH_NAME = Bead.h Polymer.h LAMMPS.h DataManager.h
LSRC_NAME = Bead.cpp Polymer.cpp LAMMPS.cpp DataManager.cpp

#
# Executable
#
EXE_NAME = Get_Bead_Position_LAMMPS
#EXE_NAME = ReadLAMMPS

#
# Set correct file directories for source, object and executable files
#
SH = $(patsubst %, $(SDIR)/%, $(SH_NAME))
SRC = $(patsubst %, $(SDIR)/%, $(SRC_NAME))
LH = $(patsubst %, $(LDIR)/%, $(LH_NAME))
LSRC = $(patsubst %, $(LDIR)/%, $(LSRC_NAME))
OBJ_NAME = $(SRC_NAME:.cpp=.o)
LOBJ_NAME = $(LSRC_NAME:.cpp=.o)
OBJ = $(patsubst %, $(OBDIR)/%, $(OBJ_NAME))
LOBJ = $(patsubst %, $(OLDIR)/%, $(LOBJ_NAME))
EXE = $(patsubst %, $(EDIR)/%, $(EXE_NAME))

#
# Compile
#
all: $(EXE)

odir:
	mkdir -p $(ODIR)

obdir:
	mkdir -p $(OBDIR)

oldir:
	mkdir -p $(OLDIR)

edir:
	mkdir -p $(EDIR)

$(OBDIR)/%.o: $(SDIR)/%.cpp obdir
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OLDIR)/%.o: $(LDIR)/%.cpp oldir
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXE): $(LOBJ) $(OBJ) edir
	$(CXX) $(CXXFLAGS) -o $(EXE) $(LOBJ) $(OBJ) $(LIB)

#
# Clean out object files and the executable
#

clean:
	rm -r $(EDIR)
	rm -r $(ODIR)
