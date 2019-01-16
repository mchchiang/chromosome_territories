/* ClusterDistributionFromWall.cpp
 * A code which determines the number of clusters in a simulation box.
 * The code stores a linked-list of bead indices. The list contains
 * multiple (cyclic) sub-linked-lists in which each list comprises of 
 * the indices of all the beads belonging to the same cluster. The code
 * is based on the algorithm described in:
 * Stoddard J. Comp. Phys. 27, 291, 1977
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;
using std::map;

double length(double x, double y, double z);

int main(int argc, char* argv[]){
  if (argc != 11){
    cout << "Usage: [numOfBeads] [lx] [ly] [lz] [dz] "
	 << "[startTime] [endTime] [timcInc] [posFile] [outFile]" 
	 << endl;
    return 1;
  }
  
  int argi {};
  const int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  const double lx {stod(string(argv[++argi]), nullptr)};
  const double ly {stod(string(argv[++argi]), nullptr)};
  const double lz {stod(string(argv[++argi]), nullptr)};
  const double dz {stod(string(argv[++argi]), nullptr)};
  const int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  const int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  const int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  const string posFile (argv[++argi]);
  const string outFile (argv[++argi]);
  
  const int numOfTypes {3};
  const double zMax {lz/2.0};

  // Set up position reader
  long time {};
  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    cout << "Problem with reading position file!" << endl;
    return 1;
  }
  
  // For creating the distribution
  int numOfBins {static_cast<int>(ceil(lz/dz))};
  vector<double> zeroVec (numOfTypes, 0.0);
  vector<vector<double> >* distribution
  {new vector<vector<double> >(numOfBins, zeroVec)};

  int index;
  double z;
  // Read the beads' positions
  while (reader.nextFrame()){
    time = reader.getTime();
    cout << time << endl;
    if (time >= startTime && time <= endTime){
      
      // Get the radial position of each bead from the CM
      for (int i {}; i < numOfBeads; i++){
	z = zMax-reader.getUnwrappedPosition(i,2);
	index = {static_cast<int>(floor(z/dz))};
	if (index >= numOfBins){
	  cout << "Bead " << i << " (z =  " << z << ") is outside max range "
	       << " z_max = " << zMax << endl;
	} else {
	  (*distribution)[index][reader.getType(i)-1] += 1.0;
	}
      }
      
    } // Close time selection
  } // Close reading loop
  
  // Normalise the distribution by the volume of each concentric shell
  double volume {dz*lx*ly};
  double timeFrames {static_cast<double>((endTime-startTime+1)/timeInc)};
  for (int i {}; i < numOfBins; i++){
    for (int t {}; t < numOfTypes; t++){
      (*distribution)[i][t] /= (volume*timeFrames);
    }
  }

  // Normalise the distribution by the maximum number density for each type
  vector<vector<double> >* normDistrb
  {new vector<vector<double> >(numOfBins, zeroVec)};

  vector<double> max (numOfTypes,0.0);
  double distrbValue;
  // Find the maximum of each type
  for (int i {}; i < numOfBins; i++){
    for (int t {}; t < numOfTypes; t++){
      distrbValue = (*distribution)[i][t];
      if (distrbValue >= max[t]){
	max[t] = distrbValue;
      }
    }
  }
  //Normalise
  for (int i {}; i < numOfBins; i++){
    for (int t {}; t < numOfTypes; t++){
      (*normDistrb)[i][t] = (*distribution)[i][t] / max[t];
    }
  }
  
  // Write the distribution to file
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Unable to open output file!" << endl;
    return 1;
  }

  double left, centre, right;
  for (int i {}; i < numOfBins; i++){
    left = i*dz;
    right = (i+1)*dz;
    centre = (left+right)/2.0;
    writer << left << " " << centre << " " << right << " ";
    for (int t {}; t < 3; t++){
      writer << (*distribution)[i][t] << " " 
	     << (*normDistrb)[i][t] << " ";
    }
    writer << endl;
  }
  writer.close();

  // Clear resources
  delete distribution;
  delete normDistrb;
}

double length(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}
