/* GenBlockCopolymer.cpp
 * This is a code that generaetes a random walk block copolymer
 * with the number of beads of each type specified by the user.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <cmath>
#include <random>
#include <algorithm>
#include "Bead.hpp"
#include "Polymer.hpp"
#include "LAMMPS.hpp"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

int main(int argc, char * argv[]){
  if (argc < 9){
    cout << "Not enough arguments! Generation aborted." << endl;
  }

  int argi {};
  int numOfBeads = stoi(string(argv[++argi]), nullptr, 10);
  int numOfType1Beads = stoi(string(argv[++argi]), nullptr, 10);
  double lx = stod(string(argv[++argi]), nullptr);
  double ly = stod(string(argv[++argi]), nullptr);
  double lz = stod(string(argv[++argi]), nullptr);
  double buffer = stod(string(argv[++argi]), nullptr);
  string outFile (argv[++argi]);
  string outMapFile (argv[++argi]);

  // Generate polymer
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>(lx, ly, lz);
  shared_ptr<Polymer> polymer {};
  shared_ptr<Bead> bead {};
  const int label {1};
  
  lammps->setTypesOfBeads(3);
  lammps->setTypesOfBonds(1);
  lammps->setTypesOfAngles(1);

  // Init random generator
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> randInt(1,2);

  // Generate the type of each bead randomly and then reshuffle them
  const int maxType1Count {numOfType1Beads};
  const int maxType2Count {numOfBeads-numOfType1Beads};
  int type1Count {};
  int type2Count {};

  vector<int> beadType (numOfBeads, 0);
  for (int i {}; i < numOfBeads; i++){
    if (type1Count >= maxType1Count){
      beadType[i] = 2;
      type2Count++;
    } else if (type2Count >= maxType2Count){
      beadType[i] = 1;
      type1Count++;
    } else {
      int r {randInt(mt)};
      if (r == 1){
	beadType[i] = 1;
	type1Count++;
      } else if (r == 2){
	beadType[i] = 2;
	type2Count++;
      }
    }
  }
  
  shuffle(beadType.begin(), beadType.end(), mt);

  polymer = lammps->createRandomWalkPolymer(label, numOfBeads, 0, 
					    0.0, 0.0, 0.0, 
					    lx-buffer, ly-buffer, lz-buffer);
  // Set bead type
  for (int i {}; i < numOfBeads; i++){
    bead = lammps->getPolymer(label)->getBead(i);
    bead->setLabel(label);
    bead->setType(beadType[i]);
  }

  // Write the input file
  lammps->exportData(outFile, outMapFile);
}
