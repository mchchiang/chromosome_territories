/* Locality.cpp
 *
 * A program that reads the position file and computes
 * the locality of the chromosome based on the average
 * number of contacts between beads
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;

int main(int argc, char* argv[]){

  if (argc < 6){
    cout << "Not enough arguments! Process aborted." << endl;
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

  CMap map = ContactMap::createFromMatrixFile(numOfBeads, numOfBeads, 
					      full, matrixFile);
  
  // Compute locality
  double local {}, distal {}, mean {};
  int localCount {}, distalCount {};
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j < numOfBeads; j++){
      if (abs(i-j) <= localDist){ // Local contact
	local += map->get(i, j);
	mean += map->get(i, j);
	localCount++;
      } else { // Distal contact
	distal += map->get(i, j);
	mean += map->get(i, j);
	distalCount++;
      }
    }
  }
  local /= static_cast<double>(localCount);
  distal /= static_cast<double>(distalCount);
  mean /= static_cast<double>(numOfBeads*numOfBeads);
  cout << "local: " << local << endl;
  cout << "distal: " << distal << endl;
  cout << "mean: " << mean << endl;

  double locality {log2(distal/local)};
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(10) << std::fixed;
  writer << locality << endl;
  writer.close();
}

