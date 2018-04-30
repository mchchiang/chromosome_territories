/* Get_ChIPSeq.cpp
 * This is a code that outputs the ChIP-seq signal for H3K9me3 and H3K27me3
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <cstdlib>
#include <cmath>
#include "Bead.hpp"
#include "Polymer.hpp"
#include "LAMMPS.hpp"
#include "DataManager.hpp"

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using std::shared_ptr;
using std::make_shared;

int getBeadType(double fracOfLam, double fracOfHet);

int main(int argc, char * argv[]){
  if (argc < 7){
    cout << "Not enough arguments! Generation aborted." << endl;
  }

  int argi {};
  string chromoFile (argv[++argi]);
  string lamFile (argv[++argi]);
  string h3k9File (argv[++argi]);
  string h3k27File (argv[++argi]);
  int chrNum {stoi(string(argv[++argi]), nullptr, 10)}; // Between 1 and 46
  string outFile (argv[++argi]);
  
  const int haploidNum {23}; // For human
  const int numOfFibres {haploidNum*2};

  // Determine the length of the fibre representing each chromosome
  int fibreLength [numOfFibres];
  const int bpPerBead {10000};

  // Read the number of bp in each chromosome
  string header {};
  ifstream chromoReader;
  chromoReader.open(chromoFile);
  
  if (!chromoReader){
    cout << "Unable to read chromo file \""
         << chromoFile << "\"\n"
         << "Aborting reading process" << endl;
    return 1;
  }

  // Skip the header of the file
  getline(chromoReader, header);
  getline(chromoReader, header);

  // Read the number of bp in each chromosome
  int index {};
  double length {};
  for (int i {}; i < numOfFibres; i++){
    chromoReader >> index >> length;
    fibreLength[i] = static_cast<int>(ceil(length / bpPerBead));
  }
  chromoReader.close();

  // Create vectors for storing laminar contact and heterochromatic region
  vector< vector<double> > fracOfLam {};
  vector< vector<double> > fracOfH3K9 {};
  vector< vector<double> > fracOfH3K27 {};

  for (int i {}; i < numOfFibres; i++){
    vector<double> vec (fibreLength[i], 0.0);
    fracOfLam.push_back(vec);
    fracOfH3K9.push_back(vec);
    fracOfH3K27.push_back(vec);
  }

  // Read the bioinformatics files (laminar and HET)
  DataManager dataMan {haploidNum, false, bpPerBead};
  bool readOK {false};

  // Read laminar file
  readOK = dataMan.getFracContent(lamFile, fracOfLam, 1, -1, 1, 2, 3, -1, 4);
  if (!readOK) return 1;

  // Read H3K9me3 (constitutive heterochromatin) file
  //  readOK = dataMan.getFracContent(hetFile, fracOfHet, 1, 0, 1, 2, 3, 5, 10);
  readOK = dataMan.getFracContent(h3k9File, fracOfH3K9, 0, 0, 0, 1, 2, 4, 10);
  if (!readOK) return 1;

  // Read H3K27me3 (facultative heterochromatin) file
  readOK = dataMan.getFracContent(h3k27File, fracOfH3K27, 0, 0, 0, 1, 2, 4, 10);
  if (!readOK) return 1;
  
  // Write the input file
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
  }
  
  for (int i {}; i < fibreLength[chrNum-1]; i++){
    writer << i << " " << fracOfLam[chrNum-1][i] << " "
	   << fracOfH3K9[chrNum-1][i] << " "
	   << fracOfH3K27[chrNum-1][i] << endl;
  }
  writer.close();
}
