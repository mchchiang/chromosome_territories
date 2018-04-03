/* ConvertMatrix.cpp
 *
 * A program which converts matrix file from one form to another.
 */

#include <iostream>
#include <string>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char* argv[]){

  if (argc < 6){
    cout << "[numOfBeads] [oldMode] [newMode] [matrixFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  string oldMode (argv[++argi]);
  string newMode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string outFile (argv[++argi]);

  bool full1 {true}, full2 {true};
  if (oldMode != "full") full1 = false;
  if (newMode != "full") full2 = false;

  CMap map = ContactMap::createFromMatrixFile(numOfBeads, full1, matrixFile);
  map->exportToFile(full2, true, true, outFile);
}
