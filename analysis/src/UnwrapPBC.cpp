/* UnwrapPBC.cpp
 * A code which unwraps the periodic boundary box
 * to get the absolute coordinates for each bead
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
  // Read input arguments
  if (argc != 7){
    cout << "Usage: UnwrapPBC [numOfBeads] [lx] [ly] [lz]"
	 << "[oldPosFile] [newPosFile]" << endl;
    return 1;
  }
  
  int argi {};
  const int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  const double lx {stod(string(argv[++argi]), nullptr)};
  const double ly {stod(string(argv[++argi]), nullptr)};
  const double lz {stod(string(argv[++argi]), nullptr)};
  const string oldPosFile (argv[++argi]);
  const string newPosFile (argv[++argi]);

  string line;
  istringstream iss;
  
  // Read and convert position file
  ifstream posReader;
  ofstream posWriter;
  posReader.open(oldPosFile);
  posWriter.open(newPosFile);

  if (!posReader){
    cout << "Problem with opening position file!" << endl;
    return 1;
  }
  if (!posWriter){
    cout << "Problem wiht opening output file!" << endl;
    return 1;
  }

  string element;
  double x, y, z;
  int t, ix, iy, iz;
  
  while (!posReader.eof()){
    // Skip header lines
    getline(posReader, line);
    posWriter << line << endl;
    getline(posReader, line);
    posWriter << line << endl;
    
    // Read bead position and type
    for (int i {}; i < numOfBeads && !posReader.eof(); i++){
      getline(posReader, line);
      iss.clear();
      iss.str(line);
      iss >> element >> x >> y >> z >> ix >> iy >> iz >> t;
      
      posWriter << element << " "
		<< x+lx*ix << " " 
		<< y+ly*iy << " " 
		<< z+lz*iz << " "
		<< endl;
    }
  }
  
  posReader.close();
  posWriter.close();
}
