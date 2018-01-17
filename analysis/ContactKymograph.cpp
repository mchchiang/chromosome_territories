/* ContactKymograph.cpp
 * A program that reads the lammpstrj file and determine
 * when beads are in contact with the wall
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;

int main(int argc, char* argv[]){
  double cutOff = stod(string(argv[1]), nullptr);
  string posFile (argv[2]);
  string outFile (argv[3]);

  // Simulation parameters
  const double lz {40.0};
  const int numOfBeads {6303};
  const double wallPos {20.0};
  const int timeInc {1000};
  const int endTime {200000};
  const int headerlines {2};

  int time {};

  ifstream reader;
  reader.open(posFile);
  if (!reader){
    cout << "Problem with reading the position file!" << endl;
    return 1;
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return 1;
  }

  string line {};
  double x {}, y {}, z {};
  int ix {}, iy {}, iz {}, t {};
  string sym {};
  istringstream iss;

  writer << std::setprecision(5) << std::fixed;

  while (!reader.eof()){
    if (time > endTime){
      break;
    }
    
    // Skip header lines
    for (int i {}; i < headerlines; i++){
      getline(reader, line);
    }
    
    // Read bead position
    for (int i {1}; i <= numOfBeads; i ++){
      getline(reader, line);
      iss.clear();
      iss.str(line);
      iss >> sym >> x >> y >> z >> ix >> iy >> iz >> t;
      z += iz * lz;
      
      if (z <= (wallPos-cutOff)){
	t = 0;
      }
      writer << time << " " << i << " " << t << endl;
    }
    writer << endl;
    time += timeInc;
  }

  reader.close();
  writer.close();
}
