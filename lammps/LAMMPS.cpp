// LAMMPS.cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include "Bead.h"
#include "Polymer.h"
#include "LAMMPS.h"

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::vector;
using std::string;
using std::map;
using std::ofstream;
using std::stringstream;

// Constructors
LAMMPS::LAMMPS() {}

LAMMPS::LAMMPS(double x, double y, double z) : lx {x}, ly {y}, lz {z} {}


// Accessor methods
shared_ptr<Polymer> LAMMPS::getPolymer(int id){
  return polymers[id];
}

shared_ptr<Bead> LAMMPS::getBead(int id){
  return beads[id];
}

void LAMMPS::setLx(double x){
  lx = x;
}

double LAMMPS::getLx(){
  return lx;
}

void LAMMPS::setLy(double y){
  ly = y;
}

double LAMMPS::getLy(){
  return ly;
}

void LAMMPS::setLz(double z){
  lz = z;
}

double LAMMPS::getLz(){
  return lz;
}

int LAMMPS::getNumOfBeads(){
  int total {};
  for (auto const& polymer : polymers){
    total += polymer->getNumOfBeads();
  }
  total += beads.size();
  return total;
}

int getNumOfBonds(){return 0;}
int getNumOfAngles(){return 0;}

void LAMMPS::setTypesOfBeads(int type){
  typesOfBeads = type;
}

int LAMMPS::getTypesOfBeads(){
  return typesOfBeads;
}

void LAMMPS::setTypesOfBonds(int type){
  typesOfBonds = type;
}

int LAMMPS::getTypesOfBonds(){
  return typesOfBonds;
}

void LAMMPS::setTypesOfAngles(int type){
  typesOfAngles = type;
}

int LAMMPS::getTypesOfAngles(){
  return typesOfAngles;
}

void LAMMPS::addBead(int id, shared_ptr<Bead> bead){
  beads.insert(beads.begin()+id, bead);
}

void LAMMPS::removeBead(int id){
  beads[id]->removeAllBonds();
  beads[id]->removeAllAngles();
  beads.erase(beads.begin()+id);
}

void LAMMPS::addPolymer(int id, shared_ptr<Polymer> polymer){
  polymers.insert(polymers.begin()+id, polymer);
}

void LAMMPS::removePolymer(int id){
  polymers[id]->removeAllBeads();
  polymers.erase(polymers.begin()+id);
}

bool LAMMPS::importData(string inFile, string mapFile){
  return true;
}

bool LAMMPS::exportData(string outFile, string mapFile){
  int beadIndexCount {1};
  int bondIndexCount {1};
  int angleIndexCount {1};

  map< shared_ptr<Bead>, int> beadIndexMap {};
  map< shared_ptr<Bead::Bond>, int> bondIndexMap {};
  map< shared_ptr<Bead::Angle>, int> angleIndexMap {};  

  cout << "Start writting LAMMPS file..." << endl;

  ofstream writer;
  writer.open(outFile);

  stringstream headerWriter;
  stringstream positionWriter;
  stringstream velocityWriter;
  stringstream bondWriter;
  stringstream angleWriter;

  positionWriter << "\nAtoms\n" << endl;
  positionWriter << std::defaultfloat;
  velocityWriter << "\nVelocities\n" << endl;
  velocityWriter << std::defaultfloat;
  bondWriter << "\nBonds\n" << endl;
  bondWriter << std::defaultfloat;
  angleWriter << "\nAngles\n" << endl;
  angleWriter << std::defaultfloat;

  /*positionWriter << std::fixed;
  velocityWriter << std::fixed;
  bondWriter << std::fixed;
  angleWriter << std::fixed;*/

  // Write positions and velocities
  for (auto const& p : polymers){
    for (auto const& b : p->getBeads()){
	  writePositionAndVelocity(b, beadIndexMap, 
							   positionWriter, velocityWriter, 
							   beadIndexCount);
	}
  }

  for (auto const& b : beads){
	writePositionAndVelocity(b, beadIndexMap, 
							   positionWriter, velocityWriter, 
							   beadIndexCount);
  }

  // Write bonds and angles
  for (auto const& p : polymers){
	for (auto const& b : p->getBeads()){
	  writeBondAndAngle(b, beadIndexMap, bondIndexMap, angleIndexMap,
						 bondWriter, angleWriter,
						 bondIndexCount, angleIndexCount);
    }
  }
 
  for (auto const& b : beads){
	writeBondAndAngle(b, beadIndexMap, bondIndexMap, angleIndexMap,
						bondWriter, angleWriter,
						bondIndexCount, angleIndexCount);
  }

  writeHeader(headerWriter,
			  beadIndexCount-1, bondIndexCount-1, angleIndexCount-1);
  
  headerWriter << endl;
  positionWriter << endl;
  velocityWriter << endl;
  bondWriter << endl;
  angleWriter << endl;

  writer << headerWriter.str();
  writer << positionWriter.str();
  writer << velocityWriter.str();
  if (bondIndexCount > 0)
	writer << bondWriter.str();
  if (angleIndexCount > 0)
	writer << angleWriter.str();
  
  writer.close();

  cout << "Finish writting LAMMPS file..." << endl;

  return true;
}

void LAMMPS::writePositionAndVelocity(const shared_ptr<Bead>& bead,
									  map< shared_ptr<Bead>, int>& beadIndexMap,
									  stringstream& positionWriter,
									  stringstream& velocityWriter,
									  int& beadIndexCount){
  if (beadIndexMap.count(bead) == 0){
	writePosition(positionWriter, bead, beadIndexCount);
	writeVelocity(velocityWriter, bead, beadIndexCount);
	beadIndexMap[bead] = beadIndexCount;
	beadIndexCount++;
  }
}

void LAMMPS::writeBondAndAngle(const shared_ptr<Bead>& bead,
							   map< shared_ptr<Bead>, int >& beadIndexMap,
							   map< shared_ptr<Bead::Bond>, int >& bondIndexMap,
							   map< shared_ptr<Bead::Angle>, int >& angleIndexMap,
							   stringstream& bondWriter,
							   stringstream& angleWriter,
							   int& bondIndexCount, int& angleIndexCount){
  vector< shared_ptr<Bead::Bond> > bondList = bead->getBonds();
  vector< shared_ptr<Bead::Angle> > angleList = bead->getAngles();
  
  for (auto const& bond : bondList){
	if (bondIndexMap.count(bond) == 0){
	  bondIndexMap[bond] = bondIndexCount;
	  int type = bond->getType();
	  int bead1Index = beadIndexMap[bond->getBead(0)];
	  int bead2Index = beadIndexMap[bond->getBead(1)];
	  writeBond(bondWriter, bondIndexCount, type, bead1Index, bead2Index);
	  bondIndexCount++;
	}
  }
  
  for (auto const& angle : angleList){
	if (angleIndexMap.count(angle) == 0){
	  angleIndexMap[angle] = angleIndexCount;
	  int type = angle->getType();
	  int bead1Index = beadIndexMap[angle->getBead(0)];
	  int bead2Index = beadIndexMap[angle->getBead(1)];
	  int bead3Index = beadIndexMap[angle->getBead(2)];
	  writeAngle(angleWriter, angleIndexCount, type, 
				 bead1Index, bead2Index, bead3Index);
	  angleIndexCount++;
	}
  }
} 

void LAMMPS::writeHeader(stringstream& writer, 
						 int nBeads, int nBonds, int nAngles){
  string header {
    "LAMMPS data file from restart file: timestep = 0,\tprocs = 1"};

  writer << header << endl;
  writer << endl;
  writer << nBeads << " atoms " << endl;
  writer << nBonds << " bonds " << endl;
  writer << nAngles << " angles " << endl;
  writer << "\n";
  writer << typesOfBeads << " atom types " << endl;
  writer << typesOfBonds << " bond types " << endl;
  writer << typesOfAngles << " angle types " << endl;
  writer << "\n";
  writer << -lx/2.0 << " " << (lx-lx/2.0) << " xlo xhi" << endl;
  writer << -ly/2.0 << " " << (ly-ly/2.0) << " ylo yhi" << endl;
  writer << -lz/2.0 << " " << (lz-lz/2.0) << " zlo zhi" << endl;
  
  writer << "\nMasses\n" << endl;
  for (int i {1}; i <= typesOfBeads; i++)
    writer << i << " " << 1 << endl;
}

void LAMMPS::writePosition(stringstream& writer,
						   const shared_ptr<Bead>& bead, int beadIndex){
  writer << beadIndex << " "
         << bead->getLabel() << " "
         << bead->getType() << " ";
  writer << std::scientific;
  writer << setprecision(preci) << bead->getPosition(0) << " "
         << setprecision(preci) << bead->getPosition(1) << " "
         << setprecision(preci) << bead->getPosition(2) << " ";
  writer << std::defaultfloat;
  writer << bead->getBoundaryCount(0) << " "
         << bead->getBoundaryCount(1) << " "
         << bead->getBoundaryCount(2) << endl;
}

void LAMMPS::writeVelocity(stringstream& writer,
						   const shared_ptr<Bead>& bead, int beadIndex){
  writer << beadIndex << " ";
  writer << std::scientific;
  writer << setprecision(preci) << bead->getVelocity(0) << " "
         << setprecision(preci) << bead->getVelocity(1) << " "
         << setprecision(preci) << bead->getVelocity(2) << endl;
  writer << std::defaultfloat;
}

void LAMMPS::writeBond(stringstream& writer, int bondIndex, int bondType,
					   int bead1Index, int bead2Index){
  writer << bondIndex << " "
         << bondType << " "
         << bead1Index << " "
         << bead2Index << endl;
}

void LAMMPS::writeAngle(stringstream& writer, int angleIndex, int angleType,
						int bead1Index, int bead2Index, int bead3Index){
  writer << angleIndex << " "
         << angleType  << " "
         << bead1Index << " "
         << bead2Index << " "
         << bead3Index << endl;
}
