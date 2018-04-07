/* ProbFromWall.cpp
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

  if (argc < 14){
    cout << "Usage: ProbFromWall [numOfBeads] [lx] [ly] [lz] "
	 << "[min] [max] [binSize] "
	 << "[beadType] [startTime] [endTime] [timeInc] "
	 << "[outFile] [posFiles]" << endl;
    return 1;
  }

  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  double lx {stod(string(argv[++argi]), nullptr)};
  double ly {stod(string(argv[++argi]), nullptr)};
  double lz {stod(string(argv[++argi]), nullptr)};
  double min {stod(string(argv[++argi]), nullptr)};
  double max {stod(string(argv[++argi]), nullptr)};
  double binSize {stod(string(argv[++argi]), nullptr)};
  string beadType (argv[++argi]);
  int startTime {stoi(string(argv[++argi]), nullptr, 10)};
  int endTime {stoi(string(argv[++argi]), nullptr, 10)};
  int timeInc {stoi(string(argv[++argi]), nullptr, 10)};
  string outFile (argv[++argi]);

  // Bead selection - selet only hetrochromatin beads or select all beads
  bool(*isSelectedType)(int){};
  if (beadType == "het"){
    isSelectedType = [](int t){return t==2;};
  } else if (beadType == "eu"){
    isSelectedType = [](int t){return t==1;};
  } else {
    isSelectedType = [](int t){return true;};
  }
  
  double wallPos {lz/2.0};
  int numOfBins {static_cast<int>(ceil((max-min)/binSize))};
  vector<double> distrb (numOfBins, 0.0);

  string posFile;
  double z;
  int col, count {};
  long time {};
  PositionReader reader;
  for (int k {++argi}; k < argc; k++){
    posFile = argv[k];
    reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
    if (!reader.isOpen()){
      cout << "Problem with reading '"<< posFile << "'!" << endl;
    }

    while (reader.nextFrame()){
      time = reader.getTime();
      if (time >= startTime && time <= endTime){
	for (int i {}; i < numOfBeads; i++){
	  // Only consider specific bead types
	  if (isSelectedType(reader.getType(i))){
	    z = wallPos-reader.getUnwrappedPosition(i, 2);
	    col = static_cast<int>(floor(z-min)/binSize);
	    if (col < 0){
	      cout << "The value " << z << " is less than min = "
		   << min << endl;
	    } else if (col >= numOfBins){
	      cout << "The value " << z << " is greater than/equal to max = "
		   << max << endl;
	    } else {
	      distrb[col] += 1.0;
	      count++;
	    }
	  }
	}
      } else if (time > endTime){
	break;
      }
    }
    reader.close();
  }

  for (int i {}; i < numOfBins; i++){
    distrb[i] /= static_cast<double>(count);
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }  

  double left, centre, right;
  writer << std::setprecision(5);
  for (int i {}; i < numOfBins; i++){
    left = i * binSize + min;
    right = (i+1) * binSize + min;
    centre = (right-left)/2.0 + left;
    writer << left << " " << centre << " " << right << " "
	   << distrb[i] << endl;
  }
  writer.close();
}
