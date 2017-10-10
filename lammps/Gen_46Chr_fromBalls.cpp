/* Gen_46Chr_fromBalls.cpp
 * This code creates a 1Mbp/bead fibre for each chromosome of
 * a human cell based on the positions of the spheres that
 * represent the chromosome in a simpler model
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <cmath>
#include "LAMMPS.h"
#include "Polymer.h"
#include "Bead.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;

int main(int argc, char * argv[]){
  // Read input arguments for data file names
  if (argc < 4){
    cout << "Not enough arguments! Generation aborted." << endl;
    return 1;
  }
  
  char* chromoFile {argv[1]};
  char* inFile {argv[2]};
  char* inMapFile {argv[3]};
  char* outFile {argv[4]};

  const int haploidNum {23}; // For human
  const int numOfFibres {haploidNum*2};
  
  // Determine the length of the fibre representing each chromosome
  long fibreLength [numOfFibres];
  const long bpPerBead {1e6};

  // Read the number of bp in each chromosome
  string header {};
  ifstream chromoReader;
  chromoReader.open(chromoFile);
  
  if (!chromoReader){
    cout << "Unable to read chromo file \""
         << chromoFile << "\"\n"
         << "Aborting reading process" << endl;
    return 1;
  }

  // Skip the header of the file                                                 
  getline(chromoReader, header);
  getline(chromoReader, header);

  for (int i {}, j {}, long length {}; i < numOfFibres; i++){
    chromoReader >> j >> length;
    fibreLength[i] = length / bpPerBead;
  }

  // Read the positions of the chromosome spheres
  bool inputDataOK {false};
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>();
  inputDataOK = lammps->importData(inFile, inMapFile);
 
  if (!inputDataOK){
    cout << "Problem with reading LAMMPS data file ... "
	 << "Generation aborted." << endl;
    return 1;
  }

  // Generate random walk polymer
  double x0 {}, y0 {}, z0 {}, x {}, y {}, z {};
  const double sphereRadius {25.0};
  const double lx {
  shared_ptr<Bead> bead {};
  for (int i {}, i < numOfFibres; i++){
    bead = lammps->getPolymer(i)->getBead(0);
    x0 = bead->getPosition(0);
    y0 = bead->getPosition(1);
    z0 = bead->getPosition(2);
    lammps->createRandomWalkPolymer(
  }

}


