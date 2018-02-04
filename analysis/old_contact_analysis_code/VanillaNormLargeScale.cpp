/* VanillaNormLargeScale.cpp
 *
 * A program that reads in a contact matrix in the format
 * (bin_i bin_j count_ij) (full or upper triangle) and
 * performs the square root vanilla normalisation on the matrix.
 * (This program differs from VanillaNorm.cpp as it handles larger
 * sparse contact matrices)
 *
 * The square root vanilla normalisation procedure is as follows:
 * M_ij* = M_ij/sqrt(sum_k M_kj * sum_l M_il)
 * For references, see Rao et al. Cell (1509) 1665-1680, 2014
 *
 * N.B. Code only works with symmetric matrices
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


void setElement(int i, int j, double value, vector<double>* contact,
		vector<int>* rowIndex, vector<int>* colStart);
double getElement(int i, int j, vector<double>* contact,
		vector<int>* rowIndex, vector<int>* colStart);

int main(int argc, char* argv[]){

  // Read in arguments
  int argi {};
  if (argc < 6){
    cout << "Not enough arguments!" << endl; 
    cout << "Usage ./VanillaNorm [size] [nx] [ny]"
	 << " [mode=full,upper] [input] [output]" << endl;
    return 1;
  }
  int nx {stoi(string(argv[++argi]), nullptr, 10)};
  int ny {stoi(string(argv[++argi]), nullptr, 10)};
  string mode (argv[++argi]);
  string matrixFile (argv[++argi]);
  string normFile (argv[++argi]);
  
  if (mode != "full" && nx != ny){
    cout << "Matrix not specified in \"full\" mode must be square!" << endl;
    return 1;
  }

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
	return 1;
      } else if (j < 0 || j >= ny){
	cout << "Index j out of range: " << j << endl;
	return 1;
      }

      (*xSum)[i] += count;
      (*ySum)[j] += count;
 
      if (mode != "full" && i != j){ // If mode is upper/lower
	(*xSum)[j] += count;
	(*ySum)[i] += count;
      }
    }
  }
  reader.close();
  
  // Do vanilla normalisation
  // Re-read the file and write the result (same as original file format)

  ofstream writer;
  writer.open(normFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(10) << std::fixed;

  reader.open(matrixFile);

  while (!reader.eof()){
    getline(reader, line);

    // Ignore any empty lines or lines begin with space or #
    if (line.size() != 0 && line[0] != ' ' && line[0] != '#'){
      iss.clear();
      iss.str(line);
      iss >> i >> j >> count;
      if ((*xSum)[i] != 0 && (*ySum)[j] != 0){
	count /= sqrt((*xSum)[i]*(*ySum)[j]);
      } else { // For cases if an entire row/col is zero
	count = 0.0;
      }
      writer << i << " " << j << " " << count << endl;
    }
  }
  reader.close();
  writer.close();  
  
  // Delete resources
  delete xSum;
  delete ySum;
}

