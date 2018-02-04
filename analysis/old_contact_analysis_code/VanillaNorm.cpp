/* VanillaNorm.cpp
 *
 * A program that reads in a contact matrix in the format
 * (bin_i bin_j count_ij) (full or upper triangle) and
 * performs the square root vanilla normalisation on the matrix.
 *
 * The square root vanilla normalisation procedure is as follows:
 * M_ij* = M_ij/sqrt(sum_k M_kj * sum_l M_il)
 * For references, see Rao et al. Cell (1509) 1665-1680, 2014
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;


int main(int argc, char* argv[]){

  // Read in arguments
  int argi {};
  if (argc < 6){
    cout << "Not enough arguments!" << endl; 
    cout << "Usage ./VanillaNorm [size] [nx] [ny] [mode=full,upper] "
	 << " [input] [output]" << endl;
    return 1;
  }
  int nx {stoi(string(argv[++argi]), nullptr, 10)};
  int ny {stoi(string(argv[++argi]), nullptr, 10)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string normFile (argv[++argi]);

  // For storing the matrix
  vector<double> zeroVec (ny, 0.0);
  vector< vector<double>>* contact =
    {new vector< vector<double>>(nx, zeroVec)};

  vector<double>* xSum = {new vector<double>(nx, 0.0)};
  vector<double>* ySum = {new vector<double>(ny, 0.0)};

  // Read the contact matrix 
  ifstream reader;
  reader.open(matrixFile);

  if (!reader){
    cout << "Problem in reading position file! Process aborted." << endl;
    return 1;
  }
  
  istringstream iss;
  int i, j;
  double count;
  string line;
  while (!reader.eof()){
    getline(reader, line);

    // Ignore any empty lines or lines begin with space or #
    if (line.size() != 0 && line[0] != ' ' && line[0] != '#'){
      iss.clear();
      iss.str(line);
      iss >> i >> j >> count;
      // Check to make sure bin indices are not out of range
      if (i < 0 || i >= nx){
	cout << "Index i out of range: " << i << endl;
      } else if (j < 0 || j >= ny){
	cout << "Index j out of range: " << j << endl;
      } else if (mode == "full"){
	(*contact)[i][j] = count;
      } else {
	(*contact)[i][j] = count;
	(*contact)[j][i] = count;
      }
    }
  }
  reader.close();
  
  // Sum the rows and cols
  for (i = 0; i < nx; i++){
    for (j = 0; j < ny; j++){
      (*xSum)[i] += (*contact)[i][j];
      (*ySum)[j] += (*contact)[i][j];
    }
  }

  // Do vanilla normalisation
  for (i = 0; i < nx; i++){
    for (j = 0; j < ny; j++){
      if ((*xSum)[i] != 0 && (*ySum)[j] != 0){
	(*contact)[i][j] /= sqrt((*xSum)[i]*(*ySum)[j]);
      } else { // For cases if an entire row/col is zero
	(*contact)[i][j] = 0.0;
      }
    }
  }
  
  // Output normalised matrix to file
  ofstream writer;
  writer.open(normFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  
  writer << std::setprecision(10) << std::fixed;
  // Output full matrix (top)
  for (i = 0; i < nx; i++){
    for (j = 0; j < ny; j++){
      writer << i << " " << j << " " << (*contact)[i][j] << endl;
    }
    writer << endl;
  }

  writer.close();

  // Delete resources
  delete contact;
  delete xSum;
  delete ySum;
}


