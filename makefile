#
# C++ compiler and options
#
CXX = g++
CXXFLAGS = -O3 -Wall -Wshadow -std=c++11
LIB = -lm

#
# Directories for source and object files and executables
#
SDIR = lammps
ODIR = obj
EDIR = exe

#
# Header files
#
HEAD_NAME = Bead.h Polymer.h LAMMPS.h

#
# Source files
#
SRC_NAME = Bead.cpp Polymer.cpp LAMMPS.cpp Gen_46Chr_1BallEach.cpp
#SRC_NAME = Bead.cpp Polymer.cpp LAMMPS.cpp Gen_46Chr_fromBall.cpp

#
# Executable
#
EXE_NAME = Gen_46Chr_1BallEach
#EXE_NAME = Gen_46Chr_fromBall

#
# Set correct file directories for source, object and executable files
#
HEAD = $(patsubst %, $(SDIR)/%, $(HEAD_NAME))
SRC = $(patsubst %, $(SDIR)/%, $(SRC_NAME))
OBJ_NAME = $(SRC_NAME:.cpp=.o)
OBJ = $(patsubst %, $(ODIR)/%, $(OBJ_NAME))
EXE = $(patsubst %, $(EDIR)/%, $(EXE_NAME))

#
# Compile
#
all: $(EXE)

odir:
	mkdir -p $(ODIR)

edir:
	mkdir -p $(EDIR)

$(ODIR)/%.o: $(SDIR)/%.cpp $(INC) odir
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXE): $(OBJ) edir
	$(CXX) $(CXXFLAGS) -o $(EXE) $(OBJ) $(LIB)

#
# Clean out object files and the executable
#

clean:
	rm -r $(EDIR)
	rm -r $(ODIR)