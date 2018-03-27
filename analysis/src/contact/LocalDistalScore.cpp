/* LocalDistalScore.cpp
 *
 * A program that calculates the local and distal score
 * and the open chromatin index (OCI) for each bead in
 * the contact map. The OCI is defined as:
 * OCI = log2(distal/local)
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using std::string;

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
  
  // Compute local and distal score
  vector<double>* localScore {new vector<double>(numOfBeads, 0.0)};
  vector<double>* distalScore {new vector<double>(numOfBeads, 0.0)};
  vector<int>* localCount {new vector<int>(numOfBeads, 0.0)};
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j < numOfBeads; j++){
      if (abs(i-j) <= localDist){
	(*localScore)[i] += map->get(i, j);
	(*localCount)[i]++;
      } else {
	(*distalScore)[i] += map->get(i, j);
      }
    }
  }
  
  // Normalise
  double count;
  for (int i {}; i < numOfBeads; i++){
    count = (*localCount)[i];
    (*localScore)[i] /= count;
    (*distalScore)[i] /= (numOfBeads-count);
  }

  // Output
  ofstream writer;
  writer.open(outFile);
  
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  writer << std::setprecision(5) << std::fixed;
  for (int i {}; i < numOfBeads; i++){
    writer << i << " " << (*localScore)[i] << " " 
	   << (*distalScore)[i] << " "
	   << log2((*distalScore)[i]/(*localScore)[i]) << endl;
  }

  // Delete resources
  delete localScore;
  delete distalScore;
  delete localCount;
}
