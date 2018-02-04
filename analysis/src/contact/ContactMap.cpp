// ContactMap.cpp

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
  string posFile (argv[++argi]);
  string contactFile (argv[++argi]);

  CMap map = ContactMap::createFromPosFile(numOfBeads, lx, ly, lz, cutoff,
					   contactType,
					   startTime, endTime, timeInc, 
					   posFile);
  map->reduceByBin(block);
  map->exportToFile(true, true, true, contactFile);
}
