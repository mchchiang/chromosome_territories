/* BeadsInSphere.cpp
 *
 * A program that reads the position file and counts                          
 * the average number of heterochromatin baeds within
 * a cutoff radius from a heterochromatin bead
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

  if (argc < 12){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int numOfBeads {stoi(string(argv[1]), nullptr, 10)};
  int lx {stoi(string(argv[2]), nullptr, 10)};
  int ly {stoi(string(argv[3]), nullptr, 10)};
  int lz {stoi(string(argv[4]), nullptr, 10)};
  double cutoff {stod(string(argv[5]), nullptr)};
  int localDist {stoi(string(argv[6]), nullptr, 10)};
  int startTime {stoi(string(argv[7]), nullptr, 10)};
  int endTime {stoi(string(argv[8]), nullptr, 10)};
  int timeInc {stoi(string(argv[9]), nullptr, 10)};
  string posFile (argv[10]);
  string outFile (argv[11]);

  const int dimension {3};
  const int headerLines {2};

  vector<double> zeroVector (dimension,0.0);
  vector< vector<double> >* position = 
    new vector< vector<double> >(numOfBeads, zeroVector);
  vector<int>* type = new vector<int>(numOfBeads, 0);

  // Beads in sphere
  vector<double>* localBeadsInSphere = new vector<double>(numOfBeads,0.0);
  vector<double>* distalBeadsInSphere = new vector<double>(numOfBeads,0.0);

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
	(*position)[i][0] = x + ix*lx;
	(*position)[i][1] = y + iy*ly;
	(*position)[i][2] = z + iz*lz;
	(*type)[i] = t;
      }
      
      // Compute contact
      double dx {}, dy {}, dz {};
      for (int i {0}; i < numOfBeads; i++){
	// Don't count self-interaction
	for (int j {}; j < i; j++){
	  dx = (*position)[i][0] - (*position)[j][0];
	  dy = (*position)[i][1] - (*position)[j][1];
	  dz = (*position)[i][2] - (*position)[j][2];
	  if (distance(dx, dy, dz) <= cutoff){
	    if (abs(i-j) < localDist){
	      (*localBeadsInSphere)[i] += 1.0;
	      (*localBeadsInSphere)[j] += 1.0;
	    } else {
	      (*distalBeadsInSphere)[i] += 1.0;
	      (*distalBeadsInSphere)[j] += 1.0;
	    }
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
  
  double avgLocalBeadsInSphere {};
  double avgDistalBeadsInSphere {};
  double avgBeadsInSphere {};
  int numOfHetBeads {};
  // Get the average beads in sphere count
  for (int i {}; i < numOfBeads; i++){
    (*localBeadsInSphere)[i] /= static_cast<double>(count);
    (*distalBeadsInSphere)[i] /= static_cast<double>(count);
    if ((*type)[i] == 2){
      avgLocalBeadsInSphere += (*localBeadsInSphere)[i];
      avgDistalBeadsInSphere += (*distalBeadsInSphere)[i];
      numOfHetBeads++;
    }
  }
  
  avgLocalBeadsInSphere /= static_cast<double>(numOfHetBeads);
  avgDistalBeadsInSphere /= static_cast<double>(numOfHetBeads);
  avgBeadsInSphere = avgLocalBeadsInSphere + avgDistalBeadsInSphere;
  
  double fracOfLocalBeads {avgLocalBeadsInSphere / avgBeadsInSphere};
  double fracOfDistalBeads {avgDistalBeadsInSphere / avgBeadsInSphere};

  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  writer << std::setprecision(10) << std::fixed;
  writer << fracOfLocalBeads << endl;

  writer.close();

  // Delete resources
  delete position;
  delete type;
  delete localBeadsInSphere;
  delete distalBeadsInSphere;
}

double distance(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}
