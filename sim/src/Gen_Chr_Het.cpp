/* Gen_Chr_Het.cpp
 * This is a code that generaetes a chromosome with heterchromatin
 * (determined by the signal of H3K9me) and laminar contact (LaminB signal)
 * info encoded
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <cmath>
#include "Bead.hpp"
#include "Polymer.hpp"
#include "LAMMPS.hpp"
#include "DataManager.hpp"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

int getBeadType(double fracOfLam, double fracOfHet);

int main(int argc, char * argv[]){
  if (argc != 11){
    cout << "Usage: Gen_Chr_Het [chromoFile] [lamFile] [hetFile] [chrNum] "
	 << "[lx] [ly] [lz] [buffer] [outFile] [outMapFile]" << endl;
    return 1;
  }

  string chromoFile (argv[1]);
  string lamFile (argv[2]);
  string hetFile (argv[3]);
  int chrNum = stoi(string(argv[4]), nullptr, 10); // Between 1 and 46
  double lx = stod(string(argv[5]), nullptr);
  double ly = stod(string(argv[6]), nullptr);
  double lz = stod(string(argv[7]), nullptr);
  double buffer = stod(string(argv[8]), nullptr);
  string outFile (argv[9]);
  string outMapFile (argv[10]);

  const int haploidNum {23}; // For human
  const int numOfFibres {haploidNum*2};

  // Determine the length of the fibre representing each chromosome
  int fibreLength [numOfFibres];
  const int bpPerBead {10000};

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

  // Read the number of bp in each chromosome
  int index {};
  double length {};
  for (int i {}; i < numOfFibres; i++){
    chromoReader >> index >> length;
    fibreLength[i] = static_cast<int>(ceil(length / bpPerBead));
  }

  // Create vectors for storing laminar contact and heterochromatic region
  vector< vector<double> > fracOfLam {};
  vector< vector<double> > fracOfHet {};

  for (int i {}; i < numOfFibres; i++){
    vector<double> vec (fibreLength[i], 0.0);
    fracOfLam.push_back(vec);
    fracOfHet.push_back(vec);
  }

  // Read the bioinformatics files (laminar and HET)
  DataManager dataMan {haploidNum, false, bpPerBead};
  bool readOK {false};
  const string chrPrefix {};
  // Read laminar file
  readOK = dataMan.getFracContent(chrPrefix, lamFile, fracOfLam, 
				  0, -1, 0, 1, 2, -1, 7);
  if (!readOK) return 1;

  // Read H3K9me3 (heterochromatin) file for IMR90 cell line
  readOK = dataMan.getFracContent(chrPrefix, hetFile, fracOfHet, 
				  0, -1, 0, 1, 2, -1, 7);
  if (!readOK) return 1;

  // Generate polymer
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>(lx, ly, lz);
  shared_ptr<Polymer> polymer {};
  shared_ptr<Bead> bead {};
  
  lammps->setTypesOfBeads(4);
  lammps->setTypesOfBonds(1);
  lammps->setTypesOfAngles(1);

  polymer = lammps->createRandomWalkPolymer(chrNum, fibreLength[chrNum-1], 0, 
					    0.0, 0.0, 0.0, 
					    lx-buffer, ly-buffer, lz-buffer);
  for (int i {}; i < fibreLength[chrNum-1]; i++){
    bead = lammps->getPolymer(chrNum)->getBead(i);
    bead->setLabel(chrNum);
    bead->setType(getBeadType(fracOfLam[chrNum-1][i], fracOfHet[chrNum-1][i]));
  }

  // Write the input file
  lammps->exportData(outFile, outMapFile);
}

int getBeadType(double fracOfLam, double fracOfHet){
  const int eu {1};
  const int het {2};

  if (fracOfLam > 0.0 || fracOfHet > 0.0){
    return het;
  }
  return eu;
}
