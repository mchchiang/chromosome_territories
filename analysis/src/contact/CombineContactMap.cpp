/* CombineContactMap.cpp
 *
 * A program which combines two contact matrices and output them to file
 */

#include <iostream>
#include <string>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char* argv[]){

  if (argc < 7){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  string mode1 (argv[++argi]);
  string matrixFile1 (argv[++argi]);
  string mode2 (argv[++argi]);
  string matrixFile2 (argv[++argi]);
  string outFile (argv[++argi]);

  bool full1 {true}, full2 {true};
  if (mode1 != "full") full1 = false;
  if (mode2 != "full") full2 = false;

  CMap map1 = ContactMap::createFromMatrixFile(numOfBeads, full1, matrixFile1);
  CMap map2 = ContactMap::createFromMatrixFile(numOfBeads, full2, matrixFile2);
  ContactMap::exportCombineMapsToFile(map1, map2, true, true, outFile);
}
