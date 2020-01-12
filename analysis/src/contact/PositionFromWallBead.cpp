/* PositionFromWallBead.cpp
 * A program that reads the lammpstrj file and determines
 * the mean position of the beads from the top simulation wall
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::vector;

int main(int argc, char* argv[]){
  if (argc < 10){
    cout << "Usage: PositionFromWall [numOfBeads] [lx] [ly] [lz] "
	 << "[startTime] [endTime] [timeInc] [outFile] [posFiles ...]" << endl;
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
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;  
  
  vector<double> zAvg (numOfBeads, 0.0);
  
  int time;
  double z;
  double wallPos {lz/2.0};
  int count {};

  argi++;
  for (; argi < argc; argi++) {
    PositionReader reader;
    reader.open(string(argv[argi]), numOfBeads, lx, ly, lz, timeInc);
    if (!reader.isOpen()){
      cout << "Problem with reading the position file!" << endl;
      cout << "Skipping to the next file" << endl;
      continue;
    }
    
    while (reader.nextFrame()){
      time = reader.getTime();
      if (time < startTime) {
	continue;
      } else if (time >= startTime && time <= endTime){
	for (int i {}; i < numOfBeads; i++){
	  z = reader.getUnwrappedPosition(i, 2);
	  zAvg[i] += (wallPos-z);
	}
	count++;
      } else {
	break;
      }
    }
    reader.close();   
  }

  // Output data
  for (int i {}; i < numOfBeads; i++) {
    zAvg[i] /= static_cast<double>(count);
    writer << i << " " <<  zAvg[i] << endl;
  }
  writer.close();
}
