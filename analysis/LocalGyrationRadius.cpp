/* LocalGyrationRadius.cpp
 *
 * A program that reads the lammpstrj file and produce a contact map
 * which shows the average contact between beads over time
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

  if (argc < 11){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int numOfBeads {stoi(string(argv[1]), nullptr, 10)};
  int lx {stoi(string(argv[2]), nullptr, 10)};
  int ly {stoi(string(argv[3]), nullptr, 10)};
  int lz {stoi(string(argv[4]), nullptr, 10)};
  int startTime {stoi(string(argv[5]), nullptr, 10)};
  int endTime {stoi(string(argv[6]), nullptr, 10)};
  int timeInc {stoi(string(argv[7]), nullptr, 10)};
  int beadLength {stoi(string(argv[8]), nullptr, 10)};
  string posFile (argv[9]);
  string outFile (argv[10]);

  const int dimension {3};
  const int headerLines {2};

  vector<double> zeroVector (dimension,0.0);
  vector< vector<double> >* position = 
    new vector< vector<double> >(numOfBeads, zeroVector);
  vector<int>* type = new vector<int>(numOfBeads, 0);
  
  ifstream reader;
  reader.open(posFile);

  if (!reader){
    cout << "Problem in reading position file! Process aborted." << endl;
    return 1;
  }

  string line {}, sym {};
  istringstream iss {};
  double x {}, y {}, z {};
  int ix {}, iy {}, iz {};
  int t {};
  long time {startTime};

  while (!reader.eof()){
    // Ignore header lines
    for (int i {}; i < headerLines; i++){
      getline(reader, line);
    }

    if (time == endTime){
      // Read bead position data - only store position of polymer beads
      for (int i {}; i < numOfBeads; i++){
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> sym >> x >> y >> z >> ix >> iy >> iz >> t;
	(*position)[i][0] = x + ix*lx;
	(*position)[i][1] = y + iy*ly;
	(*position)[i][2] = z + iz*lz;
	(*type)[i] = t;
      }
    } else if (time < endTime) {
      for (int i {}; i < numOfBeads; i++){
	getline(reader, line);
      }
    } else {
      break;
    }
    time += timeInc;
  }

  reader.close();
  
  // Get the mean position for the first nth beads
  vector<double> meanPos (dimension, 0.0);
  for (int i {}; i < beadLength; i++){
    for (int j {}; j < dimension; j++){
      meanPos[j] += (*position)[i][j];
    }
  }
  for (int i {}; i < dimension; i++){
    meanPos[i] /= static_cast<double>(beadLength);
  }
  
  // Compute local radius of gyration  
  double meanGyrRadius {};
  
  for (int k {}; k <= numOfBeads-beadLength; k++){

    double gyrRadius {}, diff {};
    for (int i {k}; i < k+beadLength; i++){
      for (int j {}; j < dimension; j++){
	diff = (*position)[i][j] - meanPos[j];
	gyrRadius += diff*diff;
      }
    }
    gyrRadius /= static_cast<double>(beadLength);
    meanGyrRadius += sqrt(gyrRadius);  

    // Update mean position
    if (k+beadLength < numOfBeads){
      // Subtract the earliest bead position
      for (int i {}; i < dimension; i++){
	meanPos[i] -= (*position)[k][i] / 
	  static_cast<double>(beadLength);
      }
      
      // Add the latest bead position
      for (int i {}; i < dimension; i++){
	meanPos[i] += (*position)[k+beadLength][i] / 
	  static_cast<double>(beadLength);
      }
    }
  }

  meanGyrRadius /= static_cast<double>(numOfBeads-beadLength+1);
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << meanGyrRadius << endl;
  writer.close();
  
  // Delete resources
  delete position;

}
