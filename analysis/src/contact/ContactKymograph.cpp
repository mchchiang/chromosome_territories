/* ContactKymograph.cpp
 * A program that reads the lammpstrj file and determine
 * when beads are in contact with the wall
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;

int main(int argc, char* argv[]){

  if (argc != 11){
    cout << "Usage: ContactKymograph [numOfBeads] [lx] [ly] [lz] "
	 << "[wallDist] [startTime] [endTime] [timeInc] "
	 << "[outFile] [posFiles]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
  double wallDist {stod(string(argv[++argi]), nullptr)};
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string posFile (argv[++argi]);
  string outFile (argv[++argi]);

  double wallPos {lz/2.0};

  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    cout << "Problem with reading the position file!" << endl;
    return 1;
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  long time {};
  double z, t;
  while (reader.nextFrame()){
    time = reader.getTime();
    if (time >= startTime && time <= endTime){
      for (int i {}; i < numOfBeads; i++){
	z = reader.getUnwrappedPosition(i, 2);
	t = reader.getType(i);
	if (z <= (wallPos-wallDist)){
	  t = 0;
	}
	writer << time << " " << i << " " << t << endl;
      }
      writer << endl;
    } else if (time > endTime){
      break;
    }
  }
  reader.close();
  writer.close();
}
