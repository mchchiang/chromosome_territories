/* Gen_46Chr_fromBalls.cpp
 * This code creates a 1Mbp/bead fibre for each chromosome of
 * a human cell based on the positions of the spheres that
 * represent the chromosome in a simpler model
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include "../lammps/LAMMPS.h"
#include "../lammps/Polymer.h"
#include "../lammps/Bead.h"
#include "../lammps/DataManager.h"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

int getBeadType(double fracOfProm, double fracOfLam, double fracOfPCG);

int main(int argc, char * argv[]){
  // Read input arguments for data file names
  if (argc < 9){
    cout << "Not enough arguments! Generation aborted." << endl;
    return 1;
  }
  
  string chromoFile(argv[1]);
  string promFile (argv[2]);
  string lamFile (argv[3]);
  string pcgFile (argv[4]);
  string inFile (argv[5]);
  string inMapFile (argv[6]);
  string outFile (argv[7]);
  string outMapFile (argv[8]);

  const int haploidNum {23}; // For human
  const int numOfFibres {haploidNum*2};
  
  // Determine the length of the fibre representing each chromosome
  int fibreLength [numOfFibres];
  const int bpPerBead {1000000};

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

  // Create vectors for storing promoter, laminar, PCG scores
  vector< vector<double> > fracOfProm {};
  vector< vector<double> > fracOfLam {};
  vector< vector<double> > fracOfPCG {};
  for (int i {}; i < numOfFibres; i++){
    vector<double> vec (fibreLength[i], 0.0);
    fracOfProm.push_back(vec);
    fracOfLam.push_back(vec);
    fracOfPCG.push_back(vec);
  }
  
  // Read the bioinformatics files (promoter, laminar, PCG)
  DataManager dataMan {haploidNum, false, bpPerBead}; 
  bool readOK {false};

  // Read promoter file
  readOK = dataMan.getFracContent(promFile, fracOfProm, 1, 0, 1, 2, 3, 5, 10);
  if (!readOK) return 1;
  
  // Read laminar file
  readOK = dataMan.getFracContent(lamFile, fracOfLam, 1, -1, 1, 2, 3, -1, 4);
  if (!readOK) return 1;
  
  // Read PCG file
  readOK = dataMan.getFracContent(pcgFile, fracOfPCG, 1, 15.0, 0, 1, 2, 3, 4);
  if (!readOK) return 1;
  
  // Read the positions of the chromosome spheres
  bool inputDataOK {false};
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>();
  inputDataOK = lammps->importData(inFile, inMapFile);
 
  if (!inputDataOK){
    cout << "Unable to read LAMMPS data file\n "
	 << "Generation aborted" << endl;
    return 1;
  }

  lammps->setTypesOfBeads(6);
  lammps->setTypesOfBonds(1);
  lammps->setTypesOfAngles(1);

  // Generate random walk polymer
  double x0 {}, y0 {}, z0 {};
  const double r {20.0};
  const double lx {r*2.0/sqrt(3.0)};
  const double ly {r*2.0/sqrt(3.0)};
  const double lz {r*2.0/sqrt(3.0)};
  shared_ptr<Bead> bead {};
  shared_ptr<Polymer> polymer {};
  int type {};

  // Get each sphere's position and generate a polymer
  // starting from that position
  for (int i {}; i < numOfFibres; i++){
    bead = lammps->getPolymer(i)->getBead(0);
    type = bead->getType();
    x0 = bead->getPosition(0);
    y0 = bead->getPosition(1);
    z0 = bead->getPosition(2);
    polymer = lammps->createRandomWalkPolymer(i, fibreLength[i], type,
					      x0, y0, z0, lx, ly, lz);
    for (int j {}; j < fibreLength[i]; j++){
      bead = polymer->getBead(j);
      bead->setLabel(i+1);
      bead->setType(getBeadType(fracOfProm[i][j], fracOfLam[i][j], 
				fracOfPCG[i][j]));
    }
  }

  // Remove all the laminar beads
  lammps->removeAllBeads();

  // Write input file
  lammps->exportData(outFile, outMapFile);
}

// Determine the type/colour of each bead
int getBeadType(double fracOfProm, double fracOfLam, double fracOfPCG){
  bool hasProm {fracOfProm > 0.0};
  bool hasLam {fracOfLam > 0.0};
  bool hasPCG {fracOfPCG > 0.0};
  
  if (hasProm && (hasLam || hasPCG)){
	if (fracOfProm > 3.0*fracOfLam)
	  hasLam = false;
	else
	  hasProm = false;
	if (fracOfProm > 3.0*fracOfPCG)
	  hasPCG = false;
	else
	  hasProm = false;
  }

  if (hasPCG && hasLam){
	if (fracOfPCG > 2.0*fracOfLam)
	  hasLam = false;
	else
	  hasPCG = false;
  }

  if (fracOfProm <= 0.5) hasProm = false;
  if (fracOfLam > 0.8){
	hasProm = false;
	hasLam = true;
	hasPCG = false;
  }

  int colour {1}; // No dominating type
  if (hasProm) colour = 2; // Promoter
  if (hasLam && !hasPCG) colour = 3; // Laminar
  if (hasPCG && !hasLam) colour = 4; // PCG
  if (hasLam && hasPCG) colour = 5; // Mixed
  return colour;
}
