/* LinearNorm.cpp
 *
 * A program that normalises the contact matrix by the
 * contact probability as a function of the linear 
 * genome distance P(l).
 */

#include <iostream>
#include <string>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char* argv[]){

  if (argc < 5){
    cout << "Usage: [numOfBeads] [mode=full/upper] "
	 << "[matrixFile] [outFile]" << endl;
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
  map->exportToFile(true, true, true, outFile);
}
