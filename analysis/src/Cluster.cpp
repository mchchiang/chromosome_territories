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

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;
using std::map;

double distance(double x, double y, double z);

int main(int argc, char* argv[]){
  if (argc < 4){
    cout << "Not enough arguments!" << endl;
    return 1;
  }
  double rc = stod(string(argv[1]), nullptr);
  string posFile (argv[2]);
  string outFile (argv[3]);
  
  // Read bead positions
  cout << "Reading beads' positions ..." << endl;
  ifstream reader;
  istringstream iss;
  reader.open(posFile);
  if (!reader){
    cout << "Unable to read position file!" << endl;
    return 1;
  }

  string line;
  double x {}, y {}, z {};
  
  int numOfBeads;
  reader >> numOfBeads;
  
  if (numOfBeads == 0){
    // Write cluster info
    cout << "Outputting data ... " << endl;
    ofstream writer;
    writer.open(outFile);
    if (!writer){
      cout << "Unable to open output file!" << endl;
      return 1;
    }
    writer << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
    writer.close();
    return 0;
  }

  vector< vector<double> >* position
  {new vector< vector<double> >()};
  position->reserve(numOfBeads);

  for (int i {}; i < numOfBeads && !reader.eof(); i++){
    getline(reader, line);
    iss.clear();
    iss.str(line);
    iss >> x >> y >> z;
    position->push_back({x, y, z});
  }
  reader.close();
  
  vector<int>* list {new vector<int>()};
  list->reserve(numOfBeads);
  
  // Initialise the linked-list
  cout << "Initialising the linked-list ..." << endl;
  for (int i {}; i < position->size(); i++){
    list->push_back(i);
  }

  // Sort the list - using the algorithm of the referenced paper
  cout << "Sorting the linked-list ..." << endl;
  int j, temp;
  double dx, dy, dz;
  int clusterCount {};
  int beadsInCluster {};
  vector<int> clusterSize;

  for (int i {}; i < list->size()-1; i++){
    if (i == list->at(i)){
      j = i;
      beadsInCluster++;
      do {
	for (int k {i+1}; k < list->size(); k++){
	  if (k != list->at(k)) continue;
	  
	  dx = position->at(j)[0] - position->at(k)[0];
	  dy = position->at(j)[1] - position->at(k)[1];
	  dz = position->at(j)[2] - position->at(k)[2];
	  
	  // Swap index j and k if r_jk <= r_c
	  if (distance(dx, dy, dz) <= rc){
	    temp = list->at(j);
	    list->at(j) = list->at(k);
	    list->at(k) = temp;
	    beadsInCluster++;
	  }
	}
	j = list->at(j);
      } while (j != i);
      clusterCount++;
      clusterSize.push_back(beadsInCluster);
      beadsInCluster = 0;
    }
  }

  // Get cluster statistics
  cout << "Getting cluster statistics ... " << endl;
  
  // Compute average cluster size
  double avgSize {}, avgSizeSq {}, varSize {}, sigmaSize {}, errorSize {};
  double nCluster {static_cast<double>(clusterCount)};
  for (int size : clusterSize){
    avgSize += size;
    avgSizeSq += size*size;
  }
  cout << "Total number of beads: " << avgSize << endl;

  avgSize /= nCluster;
  avgSizeSq /= nCluster;

  if (clusterCount > 1){
    varSize = nCluster / (nCluster-1) * (avgSizeSq - avgSize*avgSize);
    sigmaSize = sqrt(varSize);
  }
  errorSize = sigmaSize / sqrt(nCluster);
  
  // Write cluster info
  cout << "Outputting data ... " << endl;
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Unable to open output file!" << endl;
    return 1;
  }

  writer << clusterCount <<  " " 
	 << avgSize << " "
	 << sigmaSize << " "
	 << errorSize << endl;
  writer.close();

  // Delete resources
  delete position;
  delete list;
}

double distance(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}
