// AverageContactMap.cpp

#include <iostream>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;

int main(int argc, char* argv[]){

  if (argc < 13){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  int lz {stoi(string(argv[++argi]), nullptr, 10)};
  double cutoff {stod(string(argv[++argi]), nullptr)};
  int block {stoi(string(argv[++argi]), nullptr, 10)};
  string contactType (argv[++argi]);
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string contactFile (argv[++argi]);

  // For storing the averaged contact map
  CMap avgMap = ContactMap::createZeroMap(numOfBeads);
  CMap map = ContactMap::createZeroMap(numOfBeads);
  int numOfPosFiles = argc - (++argi);

  for (int i {argi}; i < argc; i++){
    string posFile (argv[i]);
    cout << "Computing contact for " << argv[i] << endl;
    map->importFromPosFile(numOfBeads, lx, ly, lz, cutoff, contactType, 
			   startTime, endTime, timeInc, posFile);
    for (int j {}; j < numOfBeads; j++){
      for (int k {}; k < numOfBeads; k++){
	avgMap->set(j, k, avgMap->get(j, k) + map->get(j, k));
      }
    }
  }
  map = nullptr;

  // Average the contact map
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j < numOfBeads; j++){
      avgMap->set(i, j, avgMap->get(i, j) / numOfPosFiles);
    }
  }

  // Reduce the resolution of the map
  avgMap->reduceByBin(block);
  
  // Output map
  avgMap->exportToFile(true, true, true, contactFile);
}
