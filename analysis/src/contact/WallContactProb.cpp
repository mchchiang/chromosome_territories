/* WallContactProb.cpp
 *
 * A program that reads the position file and 
 * computes the probabilty of a bead being at distance z from
 * the wall
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include "PositionReader.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;

int main(int argc, char* argv[]){

  if (argc != 11){
    cout << "Usage: ProbFromWallPerBead [numOfBeads] [lx] [ly] [lz] "
	 << "[wallDist] [startTime] [endTime] [timeInc] "
	 << "[posFiles] [outFile]" << endl;
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
  vector<double>* contactProb {new vector<double>(numOfBeads, 0.0)};

  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    cout << "Problem with reading '"<< posFile << "'!" << endl;
  }

  double z;
  long time {};  
  int count {};
  while (reader.nextFrame()){
    time = reader.getTime();
    if (time >= startTime && time <= endTime){
      for (int i {}; i < numOfBeads; i++){
	z = wallPos-reader.getUnwrappedPosition(i, 2);
	if (z <= wallDist){
	  (*contactProb)[i] += 1.0;
	}
      }
      count++;
    } else if (time > endTime){
      break;
    }
  }
  reader.close();

  for (int i {}; i < numOfBeads; i++){
    (*contactProb)[i] /= static_cast<double>(count);
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }  

  writer << std::setprecision(5);
  for (int i {}; i < numOfBeads; i++){
    writer << i << " " << (*contactProb)[i] << endl;
  }
  writer.close();

  // Delete resources
  delete contactProb;
}
