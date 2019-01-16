/* TypeSpecificContactMap.cpp
 *
 * A program that computes the contact probability as a 
 * function of linear genome distance P(l)
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;
using std::ifstream;
using std::ofstream;

int main(int argc, char* argv[]){

  if (argc < 7){
    cout << "Usage: [numOfBeads] [beadType] " 
	 << "[mode=full/upper] [matrixFile] [beadTypeFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int beadType {stoi(string(argv[++argi]), nullptr, 10)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string beadTypeFile (argv[++argi]);
  string outFile (argv[++argi]);

  bool full {true};
  if (mode != "full") full = false;

  // Read bead type
  ifstream reader;
  int index, type;
  vector<int> beadIndex;
  reader.open(beadTypeFile);

  while (!reader.eof()) {
    reader >> index >> type;
    if (type == beadType){
      beadIndex.push_back(index);
    }
  }

  reader.close();
  
  int numOfTypeSpecificBeads {static_cast<int>(beadIndex.size())};

  CMap map {ContactMap::createFromMatrixFile(numOfBeads, full, matrixFile)};
  
  CMap typeSpecificMap {ContactMap::createZeroMap(numOfTypeSpecificBeads)};

  for (int i {}; i < numOfTypeSpecificBeads; i++) {
    for (int j {}; j < numOfTypeSpecificBeads; j++)  {
      typeSpecificMap->set(i,j,map->get(beadIndex[i],beadIndex[j]));
    }
  }
  
  typeSpecificMap->exportToFile(true, true, true, outFile);
}
