/* CombinePosXYZ.cpp
 * A simple code that converts position file to ones compatible
 * with xyz format that is readable by vmd
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::vector;
using std::string;
using std::map;

int main(int argc, char* argv[]){
  // Read input arguments
  if (argc != 5){
    cout << "Usage: TwoToThreeState [numOfBeads] [typeFile] "
	 << "[oldPosFile] [newPosFile]" << endl;
    return 1;
  }
  
  int argi {};
  int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  string typeFile (argv[++argi]);
  string oldPosFile (argv[++argi]);
  string newPosFile (argv[++argi]);

  // xyz element map
  map<int,string> elementMap;
  elementMap[1] = "O";
  elementMap[2] = "N";
  elementMap[3] = "C";
  elementMap[4] = "H";
  elementMap[5] = "F";
  elementMap[6] = "S";

  // Read type file
  string line;
  istringstream iss;
  ifstream typeReader;
  typeReader.open(typeFile);
  if (!typeReader){
    cout << "Problem with opening type file!" << endl;
    return 1;
  }

  int index, t;
  vector<int> type (numOfBeads, 0);
  for (int i {}; i < numOfBeads; i++){
    getline(typeReader, line);
    iss.clear();
    iss.str(line);
    iss >> index >> t;
    type[index] = t;
  }
  typeReader.close();

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
  int ix, iy, iz;
  
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
      
      t = type[i];
      posWriter << elementMap[t] << " "
		<< x << " " << y << " " << z << " "
		<< ix << " " << iy << " " << iz << " " << t << endl;
    }
  }
  
  posReader.close();
  posWriter.close();
}
