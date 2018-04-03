/* ICENorm.cpp
 *
 * A program that reads in a contact matrix in the format
 * (bin_i bin_j count_ij) (full or upper triangle) and
 * performs the ICE normalisation procedure as detailed in
 * the paper Imakaev et al. Nature Methods (9) 999-1003, 2012.
 *
 */

#include <iostream>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;


int main(int argc, char* argv[]){

  // Read in arguments
  if (argc < 7){
    cout << "Usage ./ICENorm [size] [numOfIter] [threshold] "
	 << "[mode=full/upper] [matrixFile] [normFile]" << endl;
    return 1;
  }
  int argi {};
  int size {stoi(string(argv[++argi]), nullptr, 10)};
  int iter {stoi(string(argv[++argi]), nullptr, 10)};
  double threshold {stod(string(argv[++argi]), nullptr)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string normFile (argv[++argi]);

  bool full = true;
  if (mode != "full") full = false;

  CMap map = ContactMap::createFromMatrixFile(size, full, matrixFile);
  map->iceNorm(iter, threshold);
  map->exportToFile(true, true, true, normFile);
}
