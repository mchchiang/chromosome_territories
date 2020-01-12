/* Cluster.cpp
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
int getMaxIndex(const vector<int>& list);
vector<double> centreOfMass(const vector<vector<double> >& position);


int main(int argc, char* argv[]){
  if (argc != 14){
    cout << "Usage: [numOfBeads] [numOfHetBeads] [lx] [ly] [lz] [rc] [dr] "
	 << "[rMax] [startTime] [endTime] [timcInc] [posFile] [outFile]" 
	 << endl;
    return 1;
  }
  
  int argi {};
  const int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  const int numOfHetBeads {stoi(string(argv[++argi]), nullptr, 10)};
  const double lx {stod(string(argv[++argi]), nullptr)};
  const double ly {stod(string(argv[++argi]), nullptr)};
  const double lz {stod(string(argv[++argi]), nullptr)};
  const double rc {stod(string(argv[++argi]), nullptr)};
  const double dr {stod(string(argv[++argi]), nullptr)};
  const double rMax {stod(string(argv[++argi]), nullptr)};
  const int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  const int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  const int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  const string posFile (argv[++argi]);
  const string outFile (argv[++argi]);
  
  const int numOfTypes {3};

  // Set up position reader
  long time {};
  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    cout << "Problem with reading position file!" << endl;
    return 1;
  }
  
  // For finding clusters
  int beadsInCluster {};
  vector<int> clusterSize;
  vector<int> clusterStartIndex;
  vector<int>* list {new vector<int>(numOfHetBeads, 0)};
  
  // For storing heterochromatin beads' positions
  vector<double> zeroVec (numOfTypes, 0.0);
  vector<vector<double> >* position 
  {new vector<vector<double> >(numOfHetBeads, zeroVec)};

  // Store the positions of beads in the largest cluster
  vector<vector<double> >* bigClusterPosition;

  // For creating the radial distribution
  int numOfBins {static_cast<int>(ceil(rMax/dr))};
  vector<vector<double> >* distribution
  {new vector<vector<double> >(numOfBins, zeroVec)};

  // Read the beads' positions
  while (reader.nextFrame()){
    time = reader.getTime();
    cout << time << endl;
    if (time >= startTime && time <= endTime){
  
      // Get heterochromatin (K9me3) beads only
      int count {};
      for (int i {}; i < numOfBeads; i++){
	if (reader.getType(i) == 2){
	  for (int j {}; j < 3; j++){
	    (*position)[count][j] = reader.getUnwrappedPosition(i,j);
	  }
	  count++;
	}
      }
      
      // Reset list
      for (int i {}; i < numOfHetBeads; i++){
	(*list)[i] = i;
      }

      // Find clusters
      int j, temp;
      double dx, dy, dz;
      
      for (size_t i {}; i < list->size(); i++){
	if (i == (*list)[i]){
	  j = i;
	  beadsInCluster++;
	  clusterStartIndex.push_back(i);
	  do {
	    for (size_t k {i+1}; k < list->size(); k++){
	      if (k != (*list)[k]) continue;
	      
	      dx = (*position)[j][0] - (*position)[k][0];
	      dy = (*position)[j][1] - (*position)[k][1];
	      dz = (*position)[j][2] - (*position)[k][2];
	      
	      // Swap index j and k if r_jk <= r_c
	      if (length(dx, dy, dz) <= rc){
		temp = (*list)[j];
		(*list)[j] = (*list)[k];
		(*list)[k] = temp;
		beadsInCluster++;
	      }
	    }
	    j = list->at(j);
	  } while (j != i);
	  clusterSize.push_back(beadsInCluster);
	  beadsInCluster = 0;
	}
      }

      // Find the largest cluster
      int maxClusterIndex {getMaxIndex(clusterSize)};
      int maxClusterSize {clusterSize[maxClusterIndex]};
      
      bigClusterPosition = 
	new vector<vector<double> >(maxClusterSize, zeroVec);
      int startIndex {clusterStartIndex[maxClusterIndex]};
      int index {startIndex};
      count = 0;
      do {
	for (int i {}; i < 3; i++){
	  (*bigClusterPosition)[count][i] = (*position)[index][i];
	}
	count++;
	index = (*list)[index];
      } while (index != startIndex);

      // Compute the centre of mass of the largest cluster
      vector<double> cm {centreOfMass(*bigClusterPosition)};      
      cout << time << " " << cm[0] << " " << cm[1] << " " << cm[2] << endl;
      
      delete bigClusterPosition;
      bigClusterPosition = nullptr;

      // Determine the number density of beads of each kind from
      // the centre of mass of the largest heterochromatin cluster
      double r;
      // Get the radial position of each bead from the CM
      for (int i {}; i < numOfBeads; i++){
	dx = reader.getUnwrappedPosition(i,0) - cm[0];
	dy = reader.getUnwrappedPosition(i,1) - cm[1];
	dz = reader.getUnwrappedPosition(i,2) - cm[2];
	r = length(dx,dy,dz);
	index = {static_cast<int>(floor(r/dr))};
	
	if (index >= numOfBins){
	  cout << "Bead " << i << " (r =  " << r << ") is outside max range "
	       << " r_max = " << rMax << endl;
	} else {
	  (*distribution)[index][reader.getType(i)-1] += 1.0;
	}
      }
      
    } // Close time selection
  } // Close reading loop
  
  // Normalise the distribution by the volume of each concentric shell
  double volume;
  double prefactor {4.0/3.0*M_PI};
  double timeFrames {static_cast<double>((endTime-startTime+1)/timeInc)};
  for (int i {}; i < numOfBins; i++){
    volume = prefactor*(pow((i+1)*dr,3)-pow(i*dr,3));
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
  
  // Don't need the heterochromatin data anymore
  delete position; 
  delete list;
  
  // Write the distribution to file
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Unable to open output file!" << endl;
    return 1;
  }

  double left, centre, right;
  for (int i {}; i < numOfBins; i++){
    left = i*dr;
    right = (i+1)*dr;
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

int getMaxIndex(const vector<int>& vec){
  int max {}, maxIndex {};
  for (size_t i {}; i < vec.size(); i++){
    if (vec[i] > max){
      max = vec[i];
      maxIndex = i;
    }
  }
  return maxIndex;
}

vector<double> centreOfMass(const vector<vector<double> >& position){
  vector<double> cm (3, 0.0);
  for (size_t i {}; i < position.size(); i++){
    for (int j {}; j < 3; j++){
      cm[j] += position[i][j];
    }
  }
  for (int i {}; i < 3; i++){
    cm[i] /= static_cast<double>(position.size());
  }
  return cm;
}

