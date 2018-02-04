/* MSD.cpp
 * 
 * Code to compute the mean square displacement
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::vector;
using std::string;

double displacementSq(double x, double y, double z);

int main (int argc, char* argv[]){
  int numOfBeads {stoi(string(argv[1]), nullptr, 10)};
  double lx {stod(string(argv[2]), nullptr)};
  double ly {stod(string(argv[3]), nullptr)};
  double lz {stod(string(argv[4]), nullptr)};
  int startTime {stoi(string(argv[5]), nullptr, 10)};
  int endTime {stoi(string(argv[6]), nullptr, 10)};
  int timeInc {stoi(string(argv[7]), nullptr, 10)};
  string posFile (argv[8]);
  string outFile (argv[9]);

  const int headerLines {9};
  const int dimension {3};
  vector<double> zeroVector (dimension, 0.0);
  vector< vector<double> > displacement (numOfBeads, zeroVector);
  vector< vector<double> > position (numOfBeads, zeroVector);
  vector<double> msd {};
  msd.reserve(static_cast<int>((endTime-startTime)/static_cast<double>(timeInc)));

  ifstream reader;
  reader.open(posFile);

  if (!reader){
    cout << "Problem with reading position file! Process aborted." << endl;
    return 1;
  }

  double x2 {}, y2 {}, z2 {};
  double x1 {}, y1 {}, z1 {};
  double x0 {}, y0 {}, z0 {};
  int ix {}, iy {}, iz {};
  int index {}, type {}; 
  long time {};
  string line {};
  istringstream iss;
  double sum {};

  while (!reader.eof()){
    // Ignore header lines
    for (int i {}; i < headerLines; i++){
      getline(reader, line);
    }
    if (time > startTime && time <= endTime){
      sum = 0;
      for (int i {}; i < numOfBeads; i++){
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> index >> type >> x0 >> y0 >> z0 >> ix >> iy >> iz;
	x2 = x0 + ix*lx;
	y2 = y0 + iy*ly;
	z2 = z0 + iz*lz;
	x1 = position[index-1][0];
	y1 = position[index-1][1];
	z1 = position[index-1][2];
	displacement[index-1][0] += (x2-x1);
	displacement[index-1][1] += (y2-y1);
	displacement[index-1][2] += (z2-z1);
	sum += displacementSq
	  (displacement[index-1][0],
	   displacement[index-1][1],
	   displacement[index-1][2]);
	position[index-1][0] = x2;
	position[index-1][1] = y2;
	position[index-1][2] = z2;
      }
      // Compute MSD for the time step
      sum /= static_cast<double>(numOfBeads);
      msd.push_back(sum);
    } else if (time == startTime){
      for (int i {}; i < numOfBeads; i++){
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> index >> type >> x0 >> y0 >> z0 >> ix >> iy >> iz;
	position[index-1][0] = x0 + ix*lx;
	position[index-1][1] = y0 + iy*ly;
	position[index-1][2] = z0 + iz*lz;
      }
      msd.push_back(0.0);
    } else if (time < startTime){
      for (int i {}; i < numOfBeads; i++){
	getline(reader, line);
      } 
    } else {
      break;
    }
    time += timeInc;
  }
  reader.close();

  // Write MSD output file
  ofstream writer;
  writer.open(outFile);
  
  if (!writer){
    cout << "Problem with opening the output file! Process aborted." << endl;
    return 1;
  }

  for (size_t i {}; i < msd.size(); i++){
    writer << startTime + i * timeInc << " " << msd[i] << endl;
  }

  writer.close();
}

double displacementSq(double x, double y, double z){
  return x*x+y*y+z*z;
}
