/* Distance.cpp
 * A program that reads the lammpstrj file and compute
 * the end to end distance as a function of contour length N
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

double distanceSq(const double r1[3], const double r2[3]);

int main(int argc, char* argv[]){
  if (argc < 10){
    cout << "Usage: Distance [numOfBeads] [lx] [ly] [lz] "
	 << " [startTime] [endTime] [timeInc] [outFile] [posfiles...]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string outFile (argv[++argi]);
  string posFile;

  PositionReader reader;
  int time, sep;
  vector<double> endToEndDistAvg (numOfBeads, 0.0); 
  vector<long> count (numOfBeads, 0);
  double r1[3], r2[3];
  argi++;
  
  for (; argi < argc; argi++) {
    posFile = string(argv[argi]);
    reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
    if (!reader.isOpen()){
      cout << "Problem with reading the position file!" << endl;
      return 1;
    }    
    while (reader.nextFrame()){
      time = reader.getTime();
      if (time >= startTime && time <= endTime){
	for (int i {}; i < numOfBeads; i++){
	  for (int j {}; j <= i; j++){
	    // Calculating (squared) distance between bead i and j
	    for (int k {}; k < 3; k++){
	      r1[k] = reader.getUnwrappedPosition(i,k);
	      r2[k] = reader.getUnwrappedPosition(j,k);
	    }
	    sep = abs(i-j);
	    endToEndDistAvg[sep] += distanceSq(r1,r2);
	    count[sep]++;
	  }
	}
      }
    }  
    reader.close();
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;

  for (int i {}; i < numOfBeads-1; i++){
     endToEndDistAvg[i] /= static_cast<double>(count[i]);
     writer << i << " " << endToEndDistAvg[i] << endl;
  }
  
  writer.close();
}

double distanceSq(const double r1[3], const double r2[3]){
  double distSq {}, diff {};
  for (int i {}; i < 3; i++){
    diff = r1[i] - r2[i];
    distSq += diff*diff;
  }
  return distSq;
}
