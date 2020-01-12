/* PositionFromWallBeadTypeData.cpp
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
  if (argc < 12){
    cout << "Usage: PositionFromWall [numOfBeads] [lx] [ly] [lz] "
	 << "[startTime] [endTime] [timeInc] [euFile] [hcFile] [allFile] "
	 << "[posFiles ...]" << endl;
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
  string euFile (argv[++argi]);
  string hcFile (argv[++argi]);
  string allFile (argv[++argi]);
  
  ofstream euWriter, hcWriter, allWriter;
  euWriter.open(euFile);
  hcWriter.open(hcFile);
  allWriter.open(allFile);
  if (!euWriter || !hcWriter || !allWriter){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  euWriter << std::setprecision(5) << std::fixed; 
  hcWriter << std::setprecision(5) << std::fixed; 
  allWriter << std::setprecision(5) << std::fixed;  
  
  int t, time;
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
      if (time >= startTime && time <= endTime){
	for (int i {}; i < numOfBeads; i++){
	  z = wallPos - reader.getUnwrappedPosition(i, 2);
	  t = reader.getType(i);
	  if (t == 1) {
	    euWriter << z << endl;
	  } else if (t == 2) {
	    hcWriter << z << endl;
	  }
	  allWriter << z << endl;
	}
	count++;
      }
    }
    reader.close();   
  }
  euWriter.close();
  hcWriter.close();
  allWriter.close();
}
