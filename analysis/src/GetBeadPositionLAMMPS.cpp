/* GetBeadPositionLAMMPS.cpp
 * A code to get the position of beads from LAMMPS output file
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include "LAMMPS.hpp"
#include "Polymer.hpp"
#include "Bead.hpp"
#include "Angle.hpp"
#include "Bond.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::shared_ptr;
using std::make_shared;
using std::string;
using std::vector;

int main(int argc, char* argv[]){
  string inputFile (argv[1]);

  ifstream reader;
  istringstream iss;
  reader.open(inputFile);
  if (!reader){
    cout << "Unable to read input file!" << endl;
    return 1;
  }
  
  int num {};
  string line {};

  // Read input and output filenames
  string lammpsFile {}, mapFile {}, outFile {}, selectionMode {};
  reader >> lammpsFile;
  reader >> mapFile;
  reader >> outFile;

  // Read polymers to select
  int numOfSelectedPolymers {}; 
  vector<int> polymerKeys;
  reader >> numOfSelectedPolymers;
  if (numOfSelectedPolymers <= 0){
    cout << "Must select at least one polymer!" << endl;
    return 1;
  }
  polymerKeys.reserve(numOfSelectedPolymers);
  cout << "Polymers to select: ";
  for (int i {}; i < numOfSelectedPolymers; i++){
    reader >> num;
    cout << num << " ";
    polymerKeys.push_back(num);
  }
  cout << endl;

  // Read types of beads to select
  int numOfSelectedTypes {};
  vector<int> types;
  reader >> numOfSelectedTypes;
  cout << "Types to select: ";
  if (numOfSelectedTypes <= 0){
    cout << "Must select at least one type of beads!" << endl;
    return 1;
  }
  types.reserve(numOfSelectedTypes);
  for (int i {}; i < numOfSelectedTypes; i++){
    reader >> num;
    cout << num << " ";
    types.push_back(num);
  }
  cout << endl;

  // Read selection mode
  double xlo {}, xhi {}, ylo {}, yhi {}, zlo {}, zhi {}, rin {}, rout {};
  reader >> selectionMode;
  cout << selectionMode << endl;
  if (selectionMode == "cubic"){
    reader >> xlo >> xhi >> ylo >> yhi >> zlo >> zhi;
  } else if (selectionMode == "shell"){
    reader >> rin >> rout;
  } else {
    cout << "Error: incorrect selection mode" << endl;
    return 1;
  }
  reader.close();
  
  // Get the bead position
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>();
  lammps->importData(lammpsFile, mapFile);
  
  double lx {lammps->getLx()};
  double ly {lammps->getLy()};
  double lz {lammps->getLz()};

  cout << "Selecting beads ... " << endl;
  // For storing the selected beads' positions
  vector< vector<double> >* position 
  {new vector< vector<double> >()};
  
  shared_ptr<Polymer> polymer {};
  polymer = lammps->getPolymer(20);
  shared_ptr<Bead> bead {};
  int numOfBeads {};
  int ix, iy, iz;
  double x, y, z, r;
  for (int& key : polymerKeys){
    polymer = lammps->getPolymer(key);
    if (polymer == nullptr){
      cout << "Can't find polymer!" << endl;
      return 1;
    }
    numOfBeads = polymer->getNumOfBeads();
    for (int j {}; j < numOfBeads; j++){
      bead = polymer->getBead(j);
      for (int& type : types){
	if (type == bead->getType()){
	  x = bead->getPosition(0);
	  y = bead->getPosition(1);
	  z = bead->getPosition(2);
	  ix = bead->getBoundaryCount(0);
	  iy = bead->getBoundaryCount(1);
	  iz = bead->getBoundaryCount(2);

	  // Select beads within the defined geometry
	  if (selectionMode == "cubic"){
	    if (x >= xlo && x <= xhi &&
		y >= ylo && y <= yhi &&
		z >= zlo && z <= zhi)
	      position->push_back({x+ix*lx, y+iy*ly, z+iz*lz});
	  } else { // selectionMode = "shell"
	    r = sqrt(x*x+y*y+z*z);
	    if (r >= rin && r <= rout)
	      position->push_back({x+ix*lx, y+iy*ly, z+iz*lz});
	  }
	  break;
	}
      }
    }
  }

  size_t numOfSelectedBeads {position->size()};
  
  // Output the selected beads positions
  cout << "Outputting bead position ..." << endl;
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  writer << numOfSelectedBeads << endl;
  writer << std::setprecision(10) 
	 << std::scientific;
  for (size_t i {}; i < numOfSelectedBeads; i++){
    writer << position->at(i)[0] << " "
	   << position->at(i)[1] << " "
	   << position->at(i)[2] << endl;
  }
  writer.close();

  // Delete resources
  delete position;
}
