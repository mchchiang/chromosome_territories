/* LocalityDistrb.cpp
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

  if (argc < 10){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int localDist {stoi(string(argv[++argi]), nullptr, 10)};
  double min {stod(string(argv[++argi]), nullptr)};
  double max {stod(string(argv[++argi]), nullptr)};
  double binSize {stod(string(argv[++argi]), nullptr)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string localDistrbFile (argv[++argi]);
  string distalDistrbFile (argv[++argi]);

  bool full {true};
  if (mode != "full") full = false;

  CMap map = ContactMap::createFromMatrixFile(numOfBeads, numOfBeads, 
					      full, matrixFile);
  int numOfBins {static_cast<int>(ceil((max-min)/binSize))};
  double range {max-min};
  vector<double> localDistrb (numOfBins, 0.0);
  vector<double> distalDistrb (numOfBins, 0.0);
  int localBeadCount {}, distalBeadCount {};
  
  // Compute locality
  for (int i {}; i < numOfBeads; i++){
    double local {}, distal {};
    int localCount {}, distalCount {};
    for (int j {}; j < numOfBeads; j++){
      if (abs(i-j) <= localDist){ // Local contact
	local += map->get(i, j);
	localCount++;
      } else { // Distal contact
	distal += map->get(i, j);
	distalCount++;
      }
    }
    local /= static_cast<double>(localCount);
    distal /= static_cast<double>(distalCount);
    int localBinIndex {static_cast<int>(floor((local-min)/range))};
    int distalBinIndex {static_cast<int>(floor((distal-min)/range))};
    if (localBinIndex >= 0 && localBinIndex < numOfBins){
      localDistrb[localBinIndex] += 1.0;
      localBeadCount++;
    }
    if (distalBinIndex >= 0 || distalBinIndex < numOfBins){
      distalDistrb[distalBinIndex] += 1.0;
      distalBeadCount++;
    }
  }
  
  for (int i {}; i < numOfBins; i++){
    localDistrb[i] /= static_cast<double>(localBeadCount);
    distalDistrb[i] /= static_cast<double>(distalBeadCount);
  }
  
  ofstream localWriter, distalWriter;
  localWriter.open(localDistrbFile);
  if (!localWriter){
    cout << "Problem with opening the local distrb file!" << endl;
    return 1;
  }
  distalWriter.open(distalDistrbFile);
  if (!distalWriter){
    cout << "Problem with opening the distal distrb file!" << endl;
    return 1;
  }

  localWriter << std::setprecision(10) << std::fixed;
  distalWriter << std::setprecision(10) << std::fixed;
  
  for (int i {}; i < numOfBins; i++){
    localWriter << i*binSize << " " << (i+1)*binSize 
		<< " " << localDistrb[i] << endl;
    distalWriter << i*binSize << " " << (i+1)*binSize 
		 << " " << distalDistrb[i] << endl;
  }
  localWriter.close();
  distalWriter.close();
}

