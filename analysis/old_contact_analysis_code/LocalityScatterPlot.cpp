/* LocalityScatterPlot.cpp
 *
 * A program that reads the lammpstrj file and produce a contact map
 * which shows the average contact between beads over time
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;

void getAverageContactMap(int numOfBeads, 
			  double lx, double ly, double lz, double cutoff, 
			  int startTime, int endTime, int timeInc,
			  vector< vector<double> >* averageContactMap, 
			  vector<int>* type, vector<string>& files);
double distance(double x, double y, double z);

int main(int argc, char* argv[]){

  if (argc < 14){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int arg {1};
  int numOfBeads {stoi(string(argv[arg++]), nullptr, 10)};
  int lx {stoi(string(argv[arg++]), nullptr, 10)};
  int ly {stoi(string(argv[arg++]), nullptr, 10)};
  int lz {stoi(string(argv[arg++]), nullptr, 10)};
  double cutoff {stod(string(argv[arg++]), nullptr)};
  int localDist {stoi(string(argv[arg++]), nullptr, 10)};
  int startTime {stoi(string(argv[arg++]), nullptr, 10)};
  int endTime {stoi(string(argv[arg++]), nullptr, 10)};
  int timeInc {stoi(string(argv[arg++]), nullptr, 10)};
  string outFile (argv[arg++]);

  int numOfGrowFiles {stoi(string(argv[arg++]), nullptr, 10)};
  vector<string> growFiles;
  for (int i {}; i < numOfGrowFiles; i++){
    growFiles.push_back(string(argv[arg++]));
  }

  int numOfSeneFiles {stoi(string(argv[arg++]), nullptr, 10)};
  vector<string> seneFiles;
  for (int i {}; i < numOfSeneFiles; i++){
    seneFiles.push_back(string(argv[arg++]));
  }

  cout << "cutoff: " << cutoff << endl;
  cout << "localDist: " << localDist << endl;


  // Contact map matrix
  vector<double> row (numOfBeads, 0.0);
  vector< vector<double> >* growContactMap  
  {new vector< vector<double> >(numOfBeads, row)};
  vector< vector<double> >* seneContactMap
  {new vector< vector<double> >(numOfBeads, row)};
  
  vector<int>* type {new vector<int>(numOfBeads, 0)};

  // Compute contacts for growing phase
  getAverageContactMap(numOfBeads, lx, ly, lz, cutoff, 
		       startTime, endTime, timeInc, 
		       growContactMap, type, growFiles);

  // Compute contacts for senescence phase
  getAverageContactMap(numOfBeads, lx, ly, lz, cutoff, 
		       startTime, endTime, timeInc, 
		       seneContactMap, type, seneFiles);
  
  // Output
  vector<double>* growNonLocality 
  {new vector<double>(numOfBeads, 0.0)};
  vector<double>* seneNonLocality 
  {new vector<double>(numOfBeads, 0.0)};

  double growLocal {}, growDistal {};
  double seneLocal {}, seneDistal {};
  //double globalLocal {}, globalDistal {};
  int localCount {}, distalCount {};
  //int globalLocalCount {}, globalDistalCount {};

  for (int i = 0; i < numOfBeads; i++){
    for (int j = 0; j < numOfBeads; j++){
      if (abs(i-j) <= localDist){
	growLocal += (*growContactMap)[i][j];
	seneLocal += (*seneContactMap)[i][j]; 
	//	globalLocal += (*growContactMap)[i][j] + (*seneContactMap)[i][j]; 
	localCount++;
	//	globalLocalCount += 2;
      } else {
	growDistal += (*growContactMap)[i][j];
	seneDistal += (*seneContactMap)[i][j];
	//globalDistal += (*growContactMap)[i][j] + (*seneContactMap)[i][j]; 
	distalCount++;
	//globalDistalCount += 2;
      }
    }
    growLocal /= static_cast<double>(numOfBeads);
    growDistal /= static_cast<double>(numOfBeads);
    seneLocal /= static_cast<double>(numOfBeads);
    seneDistal /= static_cast<double>(numOfBeads);
    (*growNonLocality)[i] = (growDistal - growLocal)/(growDistal + growLocal);
    (*seneNonLocality)[i] = (seneDistal - seneLocal)/(seneDistal + seneLocal);
    growLocal = 0;
    seneLocal = 0;
    localCount = 0;
    distalCount = 0;
  }

  /*globalLocal /= static_cast<double>(globalLocalCount);
  globalDistal /= static_cast<double>(globalDistalCount);
  
  double overallNL = (globalDistal - globalLocal) / 
    (globalDistal + globalLocal);
  
  for (int i {}; i < numOfBeads; i++){
    (*growNonLocality)[i] = ((*growNonLocality)[i] - overallNL) / overallNL;
    (*seneNonLocality)[i] = ((*seneNonLocality)[i] - overallNL) / overallNL;
    }*/

  // Output
  ofstream writer;
  writer.open(outFile);

  if (!writer){
    cout << "Problem with opening output file! Process aborted." << endl;
    return 1;
  }
  
  writer << std::setprecision(5) << std::fixed;
  for (int i {}; i < numOfBeads; i++){
    writer << (*growNonLocality)[i] << " "
	   << (*seneNonLocality)[i] << " "
	   << (*type)[i] << endl;
  }

  writer.close();
  
  // Delete resources
  delete growContactMap;
  delete seneContactMap;
  delete growNonLocality;
  delete seneNonLocality;
}

void getAverageContactMap(int numOfBeads, 
			  double lx, double ly, double lz, double cutoff,
			  int startTime, int endTime, int timeInc,
			  vector< vector<double> >* averageContactMap, 
			  vector<int>* type, vector<string>& files){
  
  const int dimension {3};
  const int headerLines {2};

  vector<double> zeroVector (dimension,0.0);
  vector< vector<double> >* position = 
    new vector< vector<double> >(numOfBeads, zeroVector);
  
  vector<double> row (numOfBeads, 0.0);
  vector< vector<double> >* contactMap = 
    new vector< vector<double> >(numOfBeads, row);
  
  for (int k {}; k < files.size(); k++){

    cout << "Computing contact for " << files[k] << endl;
  
    ifstream reader;
    reader.open(files[k]);
    
    string line {}, sym {};
    istringstream iss {};
    double x {}, y {}, z {};
    int ix {}, iy {}, iz {};
    int t {}, count {};
    long time {};

    while (!reader.eof()){
      // Ignore header lines
      for (int i {}; i < headerLines; i++){
	getline(reader, line);
      }

      if (time >= startTime && time <= endTime){
	// Read bead position data - only store position of polymer beads
	for (int i {}; i < numOfBeads; i++){
	  getline(reader, line);
	  iss.clear();
	  iss.str(line);
	  iss >> sym >> x >> y >> z >> ix >> iy >> iz >> t;
	  (*position)[i][0] = x + ix*lx;
	  (*position)[i][1] = y + iy*ly;
	  (*position)[i][2] = z + iz*lz;
	  (*type)[i] = t;
	}
      
	// Compute contact
	double dx {}, dy {}, dz {};
	for (int i {0}; i < numOfBeads; i++){
	  (*contactMap)[i][i] += 1.0;

	  for (int j {}; j < i; j++){
	    dx = (*position)[i][0] - (*position)[j][0];
	    dy = (*position)[i][1] - (*position)[j][1];
	    dz = (*position)[i][2] - (*position)[j][2];
	    if (distance(dx, dy, dz) <= cutoff){
	      (*contactMap)[i][j] += 1.0;
	      (*contactMap)[j][i] += 1.0;
	    }
	  }
	}
	count++;

      } else if (time < startTime) {
	for (int i {}; i < numOfBeads; i++){
	  getline(reader, line);
	}
      } else {
	break;
      }
      time += timeInc;
    }

    reader.close();

    // Normalise contact map
    for (int i {}; i < numOfBeads; i++){
      for (int j {}; j < numOfBeads; j++){
	(*contactMap)[i][j] /= count;
	(*averageContactMap)[i][j] += (*contactMap)[i][j];
	(*contactMap)[i][j] = 0.0;
      }
    }
  }

  // Perform average
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j < numOfBeads; j++){
      (*averageContactMap)[i][j] /= files.size();
    }
  }

  // Delete resources
  delete position;
  delete contactMap;
}

double distance(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}
