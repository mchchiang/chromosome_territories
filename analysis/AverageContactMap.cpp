/* ContactMap.cpp
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

double distance(double x, double y, double z);

int main(int argc, char* argv[]){

  if (argc < 14){
    cout << "Not enough arguments! Process aborted." << endl;
    return 1;
  }

  int totalBeads {stoi(string(argv[1]), nullptr, 10)};
  int polymerBeads {stoi(string(argv[2]), nullptr, 10)};
  int lx {stoi(string(argv[3]), nullptr, 10)};
  int ly {stoi(string(argv[4]), nullptr, 10)};
  int lz {stoi(string(argv[5]), nullptr, 10)};
  double cutoff {stod(string(argv[6]), nullptr)};
  int block {stoi(string(argv[7]), nullptr, 10)};
  int colour {stoi(string(argv[8]), nullptr, 10)};
  int startTime {stoi(string(argv[9]), nullptr, 10)};
  int endTime {stoi(string(argv[10]), nullptr, 10)};
  int timeInc {stoi(string(argv[11]), nullptr, 10)};
  string contactFile (argv[12]);

  const int dimension {3};
  const int headerLines {2};

  vector<double> zeroVector (dimension,0.0);
  vector< vector<double> >* position = 
    new vector< vector<double> >(polymerBeads, zeroVector);
  vector<int>* type = new vector<int>(polymerBeads, 0);

  // Contact map matrix
  vector<double> row (polymerBeads, 0.0);
  vector< vector<double> >* contactMap = 
    new vector< vector<double> >(polymerBeads, row);
  vector< vector<double> >* averageContactMap = 
    new vector< vector<double> >(polymerBeads, row);
  
  for (int k {13}; k < argc; k++){

    string posFile (argv[k]);

    cout << "Computing contact for " << argv[k] << endl;

    ifstream reader;
    reader.open(posFile);

    if (!reader){
      cout << "Problem in reading position file! Process aborted." << endl;
      return 1;
    }

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
	for (int i {}; i < totalBeads; i++){
	  getline(reader, line);
	  iss.clear();
	  iss.str(line);
	  iss >> sym >> x >> y >> z >> ix >> iy >> iz >> t;
	  if (i < polymerBeads){
	    position->at(i)[0] = x + ix*lx;
	    position->at(i)[1] = y + iy*ly;
	    position->at(i)[2] = z + iz*lz;
	    type->at(i) = t;
	  }
	}
      
	// Compute contact
	double dx {}, dy {}, dz {};
	int t1 {}, t2 {};
	for (int i {0}; i < polymerBeads; i++){
	  if (colour){
	    t1 = type->at(i);
	    if (t1 == 2 || t1 == 3){
	      contactMap->at(i)[i] += 1.0;
	    } else if (t1 == 1){
	      contactMap->at(i)[i] -= 1.0;
	    }
	  } else {
	    contactMap->at(i)[i] += 1.0;
	  }

	  for (int j {}; j < i; j++){
	    dx = position->at(i)[0] - position->at(j)[0];
	    dy = position->at(i)[1] - position->at(j)[1];
	    dz = position->at(i)[2] - position->at(j)[2];
	    t1 = type->at(i);
	    t2 = type->at(j);
	    if (distance(dx, dy, dz) <= cutoff){
	      if (colour){
		if ((t1 == 2 || t1 == 3) && (t2 == 2 || t2 == 3)){
		  contactMap->at(i)[j] += 1.0;
		  contactMap->at(j)[i] += 1.0;
		} else if (t1 == 1 && t2 == 1){
		  contactMap->at(i)[j] -= 1.0;
		  contactMap->at(j)[i] -= 1.0;
		}
	      } else {
		contactMap->at(i)[j] += 1.0;
		contactMap->at(j)[i] += 1.0;
	      }
	    }
	  }
	}
	count++;

      } else if (time < startTime) {
	for (int i {}; i < totalBeads; i++){
	  getline(reader, line);
	}
      } else {
	break;
      }
      time += timeInc;
    }

    reader.close();

    // Normalise contact map
    for (int i {}; i < polymerBeads; i++){
      for (int j {}; j < polymerBeads; j++){
	contactMap->at(i)[j] /= count;
	averageContactMap->at(i)[j] += contactMap->at(i)[j];
	contactMap->at(i)[j] = 0.0;
      }
    }
    
  }

  // Perform average
  for (int i {}; i < polymerBeads; i++){
    for (int j {}; j < polymerBeads; j++){
      averageContactMap->at(i)[j] /= (argc-13);
    }
  }

  // Output contact map
  ofstream writer;
  writer.open(contactFile);

  if (!writer){
    cout << "Problem with opening output file! Process aborted." << endl;
    return 1;
  }

  // Output reduced resolution contact map
  int nRow = static_cast<int>(ceil(polymerBeads/static_cast<double>(block)));
  int nCol = nRow;
  
  for (int i {}; i < nRow; i++){
    for (int j {}; j < nCol; j++){
      // Calculate average of the contacts within the block
      int numOfSites {};
      double sum {};
      for (int k {i*block}; k < (i+1)*block; k++){
	for (int l {j*block}; l < (j+1)*block; l++){
	  if (k < polymerBeads && l < polymerBeads){
	    sum += averageContactMap->at(k)[l];
	    numOfSites++;
	  }
	}
      }
      sum /= numOfSites;
      writer << i*block << " " << j*block << " " << sum << endl;
    }
    writer << endl;
  }
  
  /*  for (int i {}; i < polymerBeads; i++){
    for (int j {}; j < polymerBeads; j++){
      writer << i << " " << j << " " << contactMap->at(i)[j] << endl;
    }
    writer << endl;
    }*/

  writer.close();
  
  

  // Delete resources
  delete position;
  delete contactMap;
  delete averageContactMap;
}

double distance(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}
