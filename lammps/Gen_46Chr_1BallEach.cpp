/* Gen_46Chr_1BallEach.cpp
 * This code generates a bead to model each chromosome of 
 * a human cell. It classifies the beads into nine types (1-9)
 * depending on their degree of contact with the nucleus laminar
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include "LAMMPS.h"
#include "Polymer.h"
#include "Bead.h"

using namespace std;

int main(int argc, char * argv[])
{
  if (argc < 2)
	{
	  cout << "Not enough arguments! Generation aborted." << endl;
	  return 1;
	}
  char * outFile {argv[1]};

  const int haploidNum {23};
  const int numOfFibres {haploidNum*2};
  const long chromoLength [numOfFibres]
     {249250621, 243199373, 198022430, 191154276, 180915260,
	  171115067, 159138663, 146364022, 141213431, 135534747,
	  135006516, 133851895, 115169878, 107349540, 102531392,
	  90354753,  81195210,  78077248,  59128983,  63025520,
	  48129895,  51304566,  155270560,
	  249250621, 243199373, 198022430, 191154276, 180915260,
	  171115067, 159138663, 146364022, 141213431, 135534747,
	  135006516, 133851895, 115169878, 107349540, 102531392,
	  90354753,  81195210,  78077248,  59128983,  63025520,
	  48129895,  51304566,  155270560
	  };
  const int numOfBeads {1};

  // Box size of simulation
  double lx {120};
  double ly {120};
  double lz {120};
  
  // Initialise the laminar contact and colour of each bead                    
  vector<double> laminar(numOfFibres, 0.0);
  vector<int> colour(numOfFibres, 0);
  
  // Reading laminar file
  cout << "Reading laminar file " << endl;

  ifstream readLAM;
  string readLAMFile {"LAD.Pk.genome.full.dat"};
  readLAM.open(readLAMFile);

  if (!readLAM){
	cout << "Unable to read the LAM file \"" 
		 << readLAMFile << "\"\n"
		 << "Aborting reading process" << endl;
	return 1;
  }

  // Remove header line
  string header;
  getline(readLAM, header); 

  const string chrPrefix {"chr"};
  double lamStart {};
  double lamEnd {};
  double lamContent {};
  string ch {};
  int bin {};
  int chromo {};

  while (!readLAM.eof()){
	readLAM >> bin >> ch >> lamStart >> lamEnd;
	ch.erase(0, chrPrefix.length());
	
	if (ch != "" && ch != "Y"){ // Ignore chrY data - model female cell
	  if (ch == "X")
		chromo = haploidNum;
	  else 
		chromo = stoi(ch, nullptr, 10);
	  
	  lamContent = (lamEnd - lamStart) / chromoLength[chromo-1];
	  laminar[chromo-1] += lamContent;
	  laminar[chromo+haploidNum-1] += lamContent;
	}
  }

  readLAM.close();

  cout << "Finish reading laminar file" << endl;

  // Assigning colour depending on laminar content
  // Binning contact strength into 10 categories (1-10)
  cout << "Assigning contact strength type to each chromosome" << endl;
  
  int numOfFibreTypes {10};
  double min {0.1};
  double max {0.6};
  double binWidth = (max - min) / numOfFibreTypes;
  int binNum {};

  for (int i {}; i < numOfFibres; i++){
	binNum = static_cast<int>((laminar[i] - min) / binWidth + 1);
	if (binNum < 0) binNum = 0;
	colour[i] = binNum;
  }
  
  // Output total laminar content for each chromosome
  cout << "The laminar content for each chromosome is as follows:" << endl;
  cout << setw(15) << "Chromosome #"
       << setw(15) << "Length [Mbp]"
       << setw(15) << "LAM content"
       << setw(15) << "Bin/Type"
       << endl;
  
  for (int i {}; i < numOfFibres; i++)
    {
      cout << setw(15) << (i+1)
		   << setw(15) << (chromoLength[i] / 1000000.0)
		   << setw(15) << laminar[i]
		   << setw(15) << colour[i]
		   << endl;
    }
  
  // Creating the polymer
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>(lx, ly, lz);
  
  // Create all polymers 
  vector< shared_ptr<Polymer> > polymers(numOfFibres);
  for (int i {}; i < numOfFibres; i++){
	polymers[i] = lammps->createPolymer(i, numOfBeads);
  }
  
  // Generate each bead's position randomly within the laminar
  
  const double buffer {10.0};
  const double maxRadius {lx/2-buffer};
  double pi {M_PI};
  
  double u {}, v {}, w {};
  double r {}, theta {}, phi {};
  double x {}, y {}, z {};
  shared_ptr<Bead> bead {};

  for (int i {}; i < numOfFibres; i++){
	u = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	v = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	w = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	r = maxRadius * pow(u, 1.0/3.0);
	theta = acos(1 - 2 * v);
	phi = 2 * pi * w;
	x = r * sin(theta) * cos(phi);
	y = r * sin(theta) * sin(phi);
	z = r * cos(theta);
	bead = polymers[i]->getBead(0);
	bead->setPosition(0, x);
	bead->setPosition(1, y);
	bead->setPosition(2, z);
	bead->setType(colour[i]);
  }

  lammps->exportData(outFile, "LAMMPS_MAP.dat");
}
