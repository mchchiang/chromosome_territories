/* OCI.cpp
 *
 * A program that calculates the Open Chromatin Index (OCI)
 * for rolling windows of fixed width in the contact map.
 *
 * The percetange of distal counts is given by
 * 
 * % of distal counts = (distal count)/(local count + distal count) * 100%
 * 
 * where all contacts within the window contribute to the local count
 * and those interactions across the window contribute to the distal count.
 *
 * The OCI is then defined to be
 *
 * OCI = % of distal counts - m
 *
 * where m is the median of the % of distal counts across all windows.
 * 
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <armadillo>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using std::string;
using std::isfinite;
using namespace arma;

int main(int argc, char* argv[]){

  if (argc < 6){
    cout << "Usage: [numOfBeads] [localDist] [mode=full/upper] "
	 << "[matrixFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int localDist {stoi(string(argv[++argi]), nullptr, 10)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string outFile (argv[++argi]);

  bool full {true};
  if (mode != "full") full = false;

  CMap map = ContactMap::createFromMatrixFile(numOfBeads, full, matrixFile);

  int numOfBins {numOfBeads-localDist+1};
  vec* oci {new vec(numOfBins, fill::zeros)};
  
  int startBead, endBead;
  double localCount, distalCount;
  for (int i {}; i < numOfBins; i++){
    startBead = i;
    endBead = i+localDist-1;
    localCount = 0.0;
    distalCount = 0.0;

    // local contacts
    for (int j {startBead}; j <= endBead; j++){
      for (int k {startBead}; k <= j; k++){
	localCount += map->get(j, k);
      }
    }

    // distal contacts
    for (int j {}; j < startBead; j++){
      for (int k {startBead}; k <= endBead; k++){
	distalCount += map->get(j, k);
      }
    }
    for (int j {startBead}; j <= endBead; j++){
      for (int k {endBead+1}; k < numOfBeads; k++){
	distalCount += map->get(j, k);
      }
    }
    
    (*oci)(i) = distalCount / (localCount + distalCount) * 100.0;
  }
  
  // Find median (ignore invalid regions)
  vector<double>* validValues {new vector<double>()};
  for (int i {}; i < numOfBins; i++){
    if (isfinite((*oci)(i))){
      validValues->push_back((*oci)(i));
    }
  }

  // Subtract % distal count by median
  double med {median(vec(*validValues))};
  oci->transform([&med](double val){return isfinite(val) ? val-med : 0.0;});

  // Output
  ofstream writer;
  writer.open(outFile);
  
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  writer << std::setprecision(5) << std::fixed;
  for (int i {}; i < numOfBins; i++){
    writer << i << " " << (*oci)(i) << endl;
  }

  // Delete resources
  delete oci;
  delete validValues;
}
