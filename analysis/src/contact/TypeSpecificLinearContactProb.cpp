/* TypeSpecificLinearContactProb.cpp
 *
 * A program that computes the contact probability as a 
 * function of linear genome distance P(l)
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
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
    cout << "Usage: [numOfBeads] [beadType] [mode=full/upper] "
	 << "[matrixFile] [beadTypeFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int type {stoi(string(argv[++argi]), nullptr, 10)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string beadTypeFile (argv[++argi]);
  string outFile (argv[++argi]);

  bool full {true};
  if (mode != "full") full = false;

  // Read bead type
  ifstream reader;
  int index, t;
  vector<int> beadType (numOfBeads, 0);
  reader.open(beadTypeFile);

  while (!reader.eof()) {
    reader >> index >> t;
    beadType[index] = t;
  }
  reader.close();
  
  CMap map = ContactMap::createFromMatrixFile(numOfBeads, full, matrixFile);

  vector<double> prob (numOfBeads, 0.0);
  vector<int> count (numOfBeads, 0);
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j <= i; j++){
      if (beadType[i] == type && beadType[j] == type) {
	prob[abs(i-j)] += map->get(i,j);
	count[abs(i-j)]++;
      }
    }
  }
  // Normalise
  for (int i {}; i < numOfBeads; i++){
    if (count[i] > 0) {
      prob[i] /= static_cast<double>(count[i]);
    } else {
      prob[i] = 0.0;
    }
  }
  
  ofstream writer;
  writer.open(outFile);
  
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  
  writer << std::setprecision(5) << std::fixed;
  for (int i {}; i < numOfBeads; i++){
    if (prob[i] > 0.0) {
      writer << i << " " << prob[i] << endl;
    }
  }
  writer.close();  
}
