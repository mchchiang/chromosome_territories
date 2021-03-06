// LAMMPS.hpp

#ifndef LAMMPS_HPP
#define LAMMPS_HPP

#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>
#include "Polymer.hpp"
#include "Bead.hpp"

using std::map;
using std::string;
using std::vector;
using std::ifstream;
using std::stringstream;

class LAMMPS {
  
private:
  map<int,shared_ptr<Polymer> > polymers {};
  map<int,shared_ptr<Bead> > beads {};
  
  // Box size
  double lx {};
  double ly {};
  double lz {};

  int typesOfBeads {1};
  int typesOfBonds {0};
  int typesOfAngles {0};

  const int preci {15}; // precision for printing floating point numbers
  
  // Internal functions
  bool readHeader (ifstream& reader, int& numOfBeads, 
		   int& numOfBonds, int& numOfAngles);

  bool readBoxSize(ifstream& reader);

  bool readPosition(ifstream& reader, int& numOfBeads,
		    map<int,shared_ptr<Bead> >& beadIndexMap);
  bool readVelocity(ifstream& reader, int& numOfBeads,
		    map<int,shared_ptr<Bead> >& beadIndexMap);

  bool readBond(ifstream& reader, int& numOfBonds,
		map<int,shared_ptr<Bead> >& beadIndexMap);
  
  bool readAngle(ifstream& reader, int& numOfAngles,
		 map<int,shared_ptr<Bead> >& beadIndexMap);

  void readInputMap(ifstream& reader,
		    map<int,shared_ptr<Bead> >& beadIndexMap);

  void writePositionAndVelocity(const shared_ptr<Bead>& bead,
				map<shared_ptr<Bead>,int >& beadIndexMap,
				stringstream& positionWriter,
				stringstream& velocityWriter,
				int& beadIndexCount);
  void writeBondAndAngle(const shared_ptr<Bead>& bead,
			 map< shared_ptr<Bead>,int >& beadIndexMap,
			 map< shared_ptr<Bond>,int >& bondIndexMap,
			 map< shared_ptr<Angle>,int >& angleIndexMap,
			 stringstream& bondWriter,
			 stringstream& angleWriter,
			 int& bondIndexCount, int& angleIndexCount);
  void writeHeader(stringstream& writer, 
		   int nBeads, int nBonds, int nAngles);
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
  LAMMPS(double x, double y, double z);

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

  void removeBead(int id);
  void removeAllBeads();
  void removePolymer(int id);
  void removeAllPolymers();
  void clear();

  bool changeBeadID(int oldID, int newID);
  bool changePolymerID(int oldID, int newID);

  bool importData(string inFile, string mapFile);
  bool exportData(string outFile, string mapFile);
  
  shared_ptr<Polymer> createPolymer(int id, int numOfBeads,
				    int beadType = 1, 
				    int bondType = 1, 
				    int angleType = 1);
  shared_ptr<Polymer> createRandomWalkPolymer(int id, int numOfBeads,
					      int beadType,
					      double x0, double y0, double z0,
					      double rx, double ry, double rz);
  shared_ptr<Bead> createBead(int id, int beadType = 1);

private:
  void addBead(int id, shared_ptr<Bead> bead);
  void addPolymer(int id, shared_ptr<Polymer> polymer);
};

#endif
