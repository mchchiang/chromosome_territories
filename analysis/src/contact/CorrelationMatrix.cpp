/* CorrelationMatrix.cpp
 *
 * A program that creates a Pearson's correlation matrix 
 * from a contact map that is normalised by the linear
 * contact probability (the contact probability as a 
 * function of linear genome distance)
 */

#include <iostream>
#include <string>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;
using std::ofstream;

int main(int argc, char* argv[]){

  if (argc < 5){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string outFile (argv[++argi]);

  bool full {true};
  if (mode != "full") full = false;

  CMap map = ContactMap::createFromMatrixFile(numOfBeads, full, matrixFile);
  map->linearProbNorm();
  map->convertToCorrelation();
  map->exportToFile(true, true, true, outFile);
}
