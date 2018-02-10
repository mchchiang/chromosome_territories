/* LinearContactProb.cpp
 *
 * A program that computes the contact probability as a 
 * function of linear genome distance P(l)
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
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
  shared_ptr<vector<double> > prob {map->getLinearProb()};
  
  ofstream writer;
  writer.open(outFile);
  
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  
  writer << std::setprecision(5) << std::fixed;
  for (int i {}; i < numOfBeads; i++){
    writer << i << " " << (*prob)[i] << endl;
  }
  writer.close();  
}
