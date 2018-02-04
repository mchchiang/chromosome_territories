/* LoopProb.cpp
 *
 * A program that reads the lammpstrj file and determine
 * the looping probability
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

double distance(double x, double y, double z);

int main(int argc, char* argv[]){

  if (argc < 11){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int numOfBeads {stoi(string(argv[1]), nullptr, 10)};
  int lx {stoi(string(argv[2]), nullptr, 10)};
  int ly {stoi(string(argv[3]), nullptr, 10)};
  int lz {stoi(string(argv[4]), nullptr, 10)};
  double cutoff {stod(string(argv[5]), nullptr)};
  int startTime {stoi(string(argv[6]), nullptr, 10)};
  int endTime {stoi(string(argv[7]), nullptr, 10)};
  int timeInc {stoi(string(argv[8]), nullptr, 10)};
  string posFile (argv[9]);
  string outFile (argv[10]);

  const int dimension {3};
  const int headerLines {2};

  vector<double> zeroVector (dimension,0.0);
  vector< vector<double> >* position = 
    new vector< vector<double> >(numOfBeads, zeroVector);
  
  // Contact map matrix
  vector<double> row (numOfBeads, 0.0);
  vector< vector<double> >* contactMap = 
    new vector< vector<double> >(numOfBeads, row);
  
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
  int t {}, count {};
  long time {};

  while (!reader.eof()){
    // Ignore header lines
    for (int i {}; i < headerLines; i++){
      getline(reader, line);
    }

    if (time >= startTime && time <= endTime){
      // Read bead position data - only store position of polymer beads
      for (int i {}; i < numOfBeads; i++){
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> sym >> x >> y >> z >> ix >> iy >> iz >> t;
	position->at(i)[0] = x + ix*lx;
	position->at(i)[1] = y + iy*ly;
	position->at(i)[2] = z + iz*lz;
      }
      
      // Compute contact
      double dx {}, dy {}, dz {};
      for (int i {0}; i < numOfBeads; i++){
	(*contactMap)[i][i] += 1.0;
	for (int j {}; j < i; j++){
	  dx = (*position)[i][0] - (*position)[j][0];
	  dy = (*position)[i][1] - (*position)[j][1];
	  dz = (*position)[i][2] - (*position)[j][2];
	  if (distance(dx, dy, dz) <= cutoff){
	    (*contactMap)[i][j] += 1.0;
	    (*contactMap)[j][i] += 1.0;
	  }
	}
      }
      count++;

    } else if (time < startTime) {
      for (int i {}; i < numOfBeads; i++){
	getline(reader, line);
      }
    } else {
      break;
    }
    time += timeInc;
  }

  reader.close();

  // Normalise contact map
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j < numOfBeads; j++){
      (*contactMap)[i][j] /= count;
    }
  }

  // Compute looping probability
  vector<double>* prob {new vector<double>(numOfBeads, 0.0)};
  for (int i {}; i < numOfBeads; i++){
    (*prob)[0] += (*contactMap)[i][i];
    for (int j {}; j < i; j++){
      (*prob)[i-j] += (*contactMap)[i][j];
    }
  }
  
  for (int i {}; i < numOfBeads; i++){
    (*prob)[i] /= (numOfBeads-i);
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  writer << std::setprecision(10) << std::fixed;
  for (int i {}; i < numOfBeads; i++){
    writer << i+1 << " " << (*prob)[i] << endl;
  }

  writer.close();

  // Delete resources
  delete prob;
  delete position;
  delete contactMap;
}

double distance(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}
