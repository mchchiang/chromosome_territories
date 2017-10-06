// LAMMPS.h

#ifndef LAMMPS_H
#define LAMMPS_H

#include <map>
#include <vector>
#include <sstream>
#include <memory>
#include <string>
#include "Polymer.h"
#include "Bead.h"

using std::map;
using std::string;
using std::vector;
using std::stringstream;

class LAMMPS {
  
private:
  vector< shared_ptr<Polymer> > polymers {};
  vector< shared_ptr<Bead> > beads {};
  
  // Box size
  double lx {};
  double ly {};
  double lz {};

  int typesOfBeads {1};
  int typesOfBonds {0};
  int typesOfAngles {0};

  const int preci {15}; // precision for printing floating point numbers
  
  // Internal functions
  void writePositionAndVelocity(const shared_ptr<Bead>& bead,
								map< shared_ptr<Bead>, int >& beadIndexMap,
								stringstream& positionWriter,
								stringstream& velocityWriter,
								int& beadIndexCount);
  void writeBondAndAngle(const shared_ptr<Bead>& bead,
						 map< shared_ptr<Bead>, int >& beadIndexMap,
						 map< shared_ptr<Bead::Bond>, int >& bondIndexMap,
						 map< shared_ptr<Bead::Angle>, int >& angleIndexMap,
						 stringstream& bondWriter,
						 stringstream& angleWriter,
						 int& bondIndexCount, int& angleIndexCount);
  void writeHeader(stringstream& writer, int nBeads, int nBonds, int nAngles);
  void writePosition(stringstream& writer, 
		     const shared_ptr<Bead>& bead, int beadIndex);
  void writeVelocity(stringstream& writer, 
		     const shared_ptr<Bead>& bead, int beadIndex);
  void writeBond(stringstream& writer, int bondIndex, int bondType, 
				 int bead1Index, int bead2Index);
  void writeAngle(stringstream& writer, int angleIndex, int angleType,
				  int bead1Index, int bead2Index, int bead3Index);

public:

  // Constructors
  LAMMPS();
  LAMMPS(double x, double y, double OBz);

  // Accessor methods
  shared_ptr<Polymer> getPolymer(int id);
  shared_ptr<Bead> getBead(int id);

  void setLx(double lx);
  double getLx();
  void setLy(double ly);
  double getLy();
  void setLz(double lz);
  double getLz();

  int getNumOfBeads();
  int getNumOfBonds();
  int getNumOfAngles();
  void setTypesOfBeads(int type);
  int getTypesOfBeads();
  void setTypesOfBonds(int type);
  int getTypesOfBonds();
  void setTypesOfAngles(int type);
  int getTypesOfAngles();
  
  void addBead(int id, shared_ptr<Bead> bead);
  void removeBead(int id);

  void addPolymer(int id, shared_ptr<Polymer> polymer);
  void removePolymer(int id);

  bool importData(string inFile, string mapFile);
  bool exportData(string outFile, string mapFile);

};

#endif
