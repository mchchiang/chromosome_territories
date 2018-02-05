/* EigenAnalysis.cpp
 *
 * Find the largest eigenvalue and the corresponding eigenvector
 * of a contact matrix (that is normalised by the linear genome
 * contact probability) - i.e. the correlation matrix
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
  map->linearProbNorm();
  map->convertToCorrelation();
  
  shared_ptr<vector<double> > eigenvec = make_shared<vector<double> >();
  double convergence {1e-10};
  double eigenval = map->maxEigen(convergence, eigenvec);
  
  ofstream writer;
  writer.open(outFile);
  
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  writer << std::setprecision(5) << std::fixed;
  writer << "# Eigenvalue = " << eigenval << endl;
  for (size_t i {}; i < eigenvec->size(); i++){
    writer << i << " " << (*eigenvec)[i] << endl;
  }
  writer.close();
}
