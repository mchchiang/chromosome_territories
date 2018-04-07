/* ContactFraction.cpp
 * A program that reads the lammpstrj file and determine
 * the fraction of beads in contact with the wall
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;

int main(int argc, char* argv[]){
  if (argc != 14){
    cout << "Usage: ContactFraction [numOfBeads] [euBeads] [hetBeads] "
	 << " [lx] [ly] [lz] [wallDist] [beadType] "
	 << " [startTime] [endTime] [timeInc] [posFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int euBeads {stoi(string(argv[++argi]), nullptr, 10)};
  int hetBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
  double wallDist {stod(string(argv[++argi]), nullptr)};
  string beadType (argv[++argi]);
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile (argv[++argi]);
  string outFile (argv[++argi]);
  
  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    cout << "Problem with reading the position file!" << endl;
    return 1;
  }
  
  // Bead selection - selet only hetrochromatin beads or select all beads
  int norm {numOfBeads};
  bool(*isSelectedType)(int){};
  if (beadType == "het"){
    isSelectedType = [](int t){return t==2;};
    norm = hetBeads;
  } else if (beadType == "eu"){
    isSelectedType = [](int t){return t==1;};
    norm = euBeads;
  } else {
    isSelectedType = [](int t){return true;};
    norm = numOfBeads;
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }
  writer << std::setprecision(5) << std::fixed;  
  
  int t, time {}, wallBeadCount {};
  double z, wallPos {lz/2.0};
  
  while (reader.nextFrame()){
    time = reader.getTime();
    if (time >= startTime && time <= endTime){
      wallBeadCount = 0;
      for (int i {}; i < numOfBeads; i++){
	z = reader.getUnwrappedPosition(i, 2);
	t = reader.getType(i);
	if (z > (wallPos-wallDist) && isSelectedType(t)){
	  wallBeadCount++;
	}
      }
      writer << time << " " 
	     << static_cast<double>(wallBeadCount) / 
	static_cast<double>(norm) << endl;
    }
  }

  reader.close();
  writer.close();
}
