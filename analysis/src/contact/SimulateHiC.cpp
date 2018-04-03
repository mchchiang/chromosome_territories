/* SimulateHiC.cpp
 *
 * A program that simulates the HiC procedure to generate a contact map
 * from a LAMMPS dump/position file.
 */

#include <iostream>
#include <string>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char* argv[]){
  
  if (argc < 13){
    cout << "Usage: [numOfBeads] [lx] [ly] [lz] [cutoff] [block] "
	 << "[numOfCounts] [startTime] [endTime] [timeInc] "
	 << "[posFile] [contactFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int lx {stoi(string(argv[++argi]), nullptr, 10)};
  int ly {stoi(string(argv[++argi]), nullptr, 10)};
  int lz {stoi(string(argv[++argi]), nullptr, 10)};
  double cutoff {stod(string(argv[++argi]), nullptr)};
  int block {stoi(string(argv[++argi]), nullptr, 10)};
  long numOfCounts {stol(string(argv[++argi]), nullptr, 10)};
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile (argv[++argi]);
  string contactFile (argv[++argi]);

  CMap map = ContactMap::createFromSimulatedHiC(numOfBeads, lx, ly, lz, cutoff,
						numOfCounts,
						startTime, endTime, timeInc, 
						posFile);
  map->reduceByBin(block);
  map->exportToFile(true, true, true, contactFile);
}
