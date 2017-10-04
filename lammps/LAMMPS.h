// LAMMPS.h

#ifndef LAMMPS_H
#define LAMMPS_H

#include <vector>
#include <sstream>
#include <memory>
#include <string>
#include "Polymer.h"
#include "Bead.h"
#include "Bond.h"
#include "Angle.h"

using std::string;
using std::vector;
using std::stringstream;

class LAMMPS {
  
private:
  vector< shared_ptr<Polymer> > polymers;
  vector< shared_ptr<Bead> > beads;
  
  // Box size
  double lx;
  double ly;
  double lz;

  int typesOfBeads {1};
  int typesOfBonds {0};
  int typesOfAngles {0};

  const int preci {16}; // precision for printing numbers
  
  // Internal functions
  void writePosition(stringstream& writer, 
		     const shared_ptr<Bead>& bead, int beadIndex);
  void writeVelocity(stringstream& writer, 
		     const shared_ptr<Bead>& bead, int beadIndex);
  void writeBond(stringstream& writer,
		 const shared_ptr<Bond>& bond, int bondIndex);
  void writeAngle(stringstream& writer,
		  const shared_ptr<Angle>& angle, int angleIndex);

public:
  
  // Accessor methods
  shared_ptr<Polymer> getPolymer(int id);
  shared_ptr<Bead> getBeadr(int id);

  double getLx();
  double getLy();
  double getLz();

  int getNumOfBeads();
  int getNumOfBonds();
  int getNumOfAngles();
  int getTypesOfBeads();
  int getTypesOfBonds();
  int getTypesOfAngles();

  int addPolymer(shared_ptr<Polymer> polymer);
  void removePolymer(int id);

  bool importData(string inFile, string mapFile);
  bool exportData(string outFile, string mapFile);

};

#endif
