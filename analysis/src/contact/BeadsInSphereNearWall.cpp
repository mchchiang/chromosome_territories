/* BeadsInSphere.cpp
 *
 * A program that reads the position file and counts                          
 * the average number of heterochromatin baeds within
 * a cutoff radius from a heterochromatin bead
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;

const double pi {M_PI};
double length(double x, double y, double z);
double enclosedVolume(double radius, double distOverWall);

int main(int argc, char* argv[]){

  if (argc < 13){
    cout << "Usage: [numOfBeads] [lx] [ly] [lz] [cutoff] [wallDist] "
	 << "[localDist] [startTime] [endTime] [timeInc] "
	 << "[posFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
  double cutoff {stod(string(argv[++argi]), nullptr)};
  double wallDist {stod(string(argv[++argi]), nullptr)};
  int localDist {stoi(string(argv[++argi]), nullptr, 10)};
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile (argv[++argi]);
  string outFile (argv[++argi]);

  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    return 1;
  }

  // Beads in sphere
  vector<double>* avgZPos {new vector<double>(numOfBeads, 0.0)};
  vector<double>* localBeadsInSphere {new vector<double>(numOfBeads, 0.0)};
  vector<double>* distalBeadsInSphere {new vector<double>(numOfBeads, 0.0)};
  vector<double>* localDensity {new vector<double>(numOfBeads, 0.0)};
  vector<double>* distalDensity {new vector<double>(numOfBeads, 0.0)};
  vector<int>* localCount {new vector<int>(numOfBeads, 0)};
  vector<int>* distalCount {new vector<int>(numOfBeads, 0)};

  int count {};
  long time {};
  while (reader.nextFrame()){
    time = reader.getTime();
    if (time >= startTime && time <= endTime){

      // Reset the count
      std::fill(localCount->begin(), localCount->end(), 0);
      std::fill(distalCount->begin(), distalCount->end(), 0);

      // Compute contact
      double dx, dy, dz, volume;
      for (int i {}; i < numOfBeads; i++){
	(*localCount)[i]++;
	for (int j {}; j < i; j++){
	  dx = reader.getUnwrappedPosition(i, 0) - 
	    reader.getUnwrappedPosition(j, 0);
	  dy = reader.getUnwrappedPosition(i, 1) - 
	    reader.getUnwrappedPosition(j, 1);
	  dz = reader.getUnwrappedPosition(i, 2) - 
	    reader.getUnwrappedPosition(j, 2);
	  if (length(dx, dy, dz) <= cutoff){
	    if (abs(i-j) <= localDist){
	      (*localCount)[i]++;
	      (*localCount)[j]++;
	    } else {
	      (*distalCount)[i]++;
	      (*distalCount)[j]++;
	    }
	  }
	}
      }
      
      double z;
      for (int i {}; i < numOfBeads; i++){
	z = lz/2.0-reader.getUnwrappedPosition(i, 2);
	(*avgZPos)[i] += z;
	(*localBeadsInSphere)[i] += (*localCount)[i];
	(*distalBeadsInSphere)[i] += (*distalCount)[i];
	volume = enclosedVolume(cutoff, cutoff-z);
	(*localDensity)[i] += static_cast<double>((*localCount)[i])/volume;
	(*distalDensity)[i] += static_cast<double>((*distalCount)[i])/volume;
      }
      count++;
    } else if (time > endTime){
      break;
    }
  }

  // Average quantities by time
  for (int i {}; i < numOfBeads; i++){
    (*avgZPos)[i] /= static_cast<double>(count);
    (*localBeadsInSphere)[i] /= static_cast<double>(count);
    (*distalBeadsInSphere)[i] /= static_cast<double>(count);
    (*localDensity)[i] /= static_cast<double>(count);
    (*distalDensity)[i] /= static_cast<double>(count);
  }
  
  // Average quantites by beads (for those near the wall)
  double avgLocalBeadsInSphere {}, avgDistalBeadsInSphere {};
  double avgLocalDensity {}, avgDistalDensity {};
  double avgBeadsInSphere {}, avgDensity {};
  int numOfWallBeads {};
  for (int i {}; i < numOfBeads; i++){
    if ((*avgZPos)[i] <= wallDist){
      avgLocalBeadsInSphere += (*localBeadsInSphere)[i];
      avgDistalBeadsInSphere += (*distalBeadsInSphere)[i];
      avgLocalDensity += (*localDensity)[i];
      avgDistalDensity += (*distalDensity)[i];
      numOfWallBeads++;
    }
  }
  
  if (numOfWallBeads > 0){
    avgLocalBeadsInSphere /= static_cast<double>(numOfWallBeads);
    avgDistalBeadsInSphere /= static_cast<double>(numOfWallBeads);
    avgLocalDensity /= static_cast<double>(numOfWallBeads);
    avgDistalDensity /= static_cast<double>(numOfWallBeads);
    avgBeadsInSphere = avgLocalBeadsInSphere + avgDistalBeadsInSphere;
    avgDensity = avgLocalDensity + avgDistalDensity;
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  writer << std::setprecision(10) << std::fixed;
  writer << avgLocalBeadsInSphere << " "
	 << avgDistalBeadsInSphere << " "
	 << avgBeadsInSphere << " "
	 << avgLocalDensity << " "
	 << avgDistalDensity << " "
	 << avgDensity << endl;

  writer.close();

  // Delete resources
  delete avgZPos;
  delete localCount;
  delete distalCount;
  delete localBeadsInSphere;
  delete distalBeadsInSphere;
  delete localDensity;
  delete distalDensity;
}

double length(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}

double enclosedVolume(double r, double h){
  double sphereVolume {4.0/3.0*pi*r*r*r};
  double sphereCapVolume {};
  if (h > 0.0){
    sphereCapVolume = pi*h*h*(3*r-h)/3.0;
  }
  return sphereVolume - sphereCapVolume;
}

