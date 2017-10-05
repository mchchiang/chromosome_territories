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
using std::vector;
using std::string;
using std::map;
using std::ofstream;
using std::stringstream;
using Bond = Bead::Bond;
using Angle = Bead::Angle;

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

int getNumOfBonds(){}
int getNumOfAngles(){}

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

void LAMMPS::addPolymer(int id, shared_ptr<Polymer> polymer){
  polymers.insert(polymers.begin()+id, polymer);
}

void LAMMPS::removePolymer(int id){
  polymers.erase(polymers.begin()+id);
}

bool importData(string inFile, string mapFile){}

bool exportData(string outFile, string mapFile){
  int beadIndexCount {1};
  int bondIndexCount {1};
  int angleIndexCount {1};

  map< shared_ptr<Bead>, int> indexMap;
  map< shared_ptr<Bond>, int> bondIndexMap;
  map< shared_ptr<Angle>, int> angleIndexMap;

  string header {
    "LAMMPS data file from restart file: timestep = 0,\tprocs = 1"};

  cout << "Start writting LAMMPS file..." << endl;

  stringstream headerWriter;
  
  // Write atoms' positions and velocities
  stringstream positionWriter;
  stringstream velocityWriter;
  
  positionWriter << "\nAtoms\n" << endl;
  positionWriter << std::defaultfloat;
  velocityWriter << "\nVelocities\n" << endl;
  velocityWriter << std::defaultfloat;

  
  for (auto const& p : polymers){
    for (auto const& b : p->getBeads()){
      
    }
  }
  
  return true;
}

void LAMMPS::writePosition(stringstream& writer,
			   const shared_ptr<Bead>& bead, int beadIndex){
  
}

