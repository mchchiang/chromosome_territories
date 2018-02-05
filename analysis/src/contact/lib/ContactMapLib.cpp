// ContactMapLib.cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <armadillo>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;
using std::unique_ptr;
using std::shared_ptr;
using std::move;
using std::make_shared;
using namespace arma;

// Private constructor
ContactMap::ContactMap(int n, int value){
  size = n;
  unique_ptr<mat> m {new mat(size, size, fill::zeros)}; 
  contact = std::move(m);
  //vector<double> zeroVec (size, value);
  //contact = new vector<vector<double> >(size, zeroVec);
}

// Create contact map
CMap ContactMap::createZeroMap(int n){
  CMap map {new ContactMap(n)};
  return map;
}

CMap ContactMap::createFromArray(vector<vector<double> >* matrix){
  int n = (*matrix).size();
  CMap map {new ContactMap(n)};
  map->importFromArray(matrix);
  return map;
}

CMap
ContactMap::createFromPosFile(int numOfBeads, double lx, double ly, double lz,
			      double cutoff, string contactType,
			      int startTime, int endTime, int timeInc,
			      string file){
  CMap map {new ContactMap(numOfBeads)};
  map->importFromPosFile(numOfBeads, lx, ly, lz, cutoff, contactType,
			 startTime, endTime, timeInc, file);
  return map;
}

CMap ContactMap::createFromMatrixFile(int n, bool full,
				      string file){
  CMap map {new ContactMap(n)};
  map->importFromMatrixFile(n, full, file);
  return map;
}


void ContactMap::importFromPosFile(int numOfBeads, 
				   double lx, double ly, double lz,
				   double cutoff, string contactType,
				   int startTime, int endTime, int timeInc,
				   string file){
  
  // Useful constants for reading the position file
  const int dimensions {3};
  const int headerLines {2};
  
  vector<double> zeroVec (dimensions,0.0);
  vector<vector<double> >* position = 
    new vector<vector<double> >(numOfBeads, zeroVec);
  vector<int>* type = new vector<int>(numOfBeads, 0);

  ifstream reader;
  reader.open(file);
  
  // Check that the file exists and can be read
  if (!reader){
    cout << "Problem in reading position file!" << endl;
    return;
  }

  // Reset the contact map
  reset(numOfBeads);
  
  string line, sym;
  double x, y, z;
  int ix, iy, iz, t;
  istringstream iss;
  int count {};
  long time {};

  while (!reader.eof()){
    // Ignore header lines
    for (int i {}; i < headerLines; i++){
      getline(reader, line);
    }

    if (time >= startTime && time <= endTime){
      // Read bead position data
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
      if (contactType == "colour"){
	computeColourContact(cutoff, position, type);
      } else {
	computeContact(cutoff, position);
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
      set(i, j, get(i, j) / count);
    }
  }
  
  // Delete resources
  delete position;
  delete type;
}

void ContactMap::importFromMatrixFile(int n, bool full, string file){
  ifstream reader;
  reader.open(file);
  
  // Check that the file exists and can be read
  if (!reader){
    cout << "Problem in reading position file!" << endl;
    return;
  }

  // Reset the contact map
  reset(n);

  // Read the contact map values  
  istringstream iss;
  int i, j;
  double count;
  string line;
  while (!reader.eof()){
    getline(reader, line);

    // Ignore any empty lines or lines begin with space or #
    if (line.size() != 0 && line[0] != ' ' && line[0] != '#'){
      iss.clear();
      iss.str(line);
      iss >> i >> j >> count;
      // Check to make sure bin indices are not out of range
      if (i < 0 || i >= size){
	cout << "Index i out of range: " << i << endl;
      } else if (j < 0 || j >= size){
	cout << "Index j out of range: " << j << endl;
      } 

      set(i, j, count);
      if (!full){
	set(j, i, count);
      }
    }
  }
  reader.close();
 
}

void ContactMap::importFromArray(vector<vector<double> >* matrix){
  int n = (*matrix).size();
  
  // Reset the contact map
  reset(n);

  for (int i {}; i < n; i++){
    for (int j {}; j < n; j++){
      set(i, j, (*matrix)[i][j]);
    }
  }
}

// Compute normal contact
void ContactMap::computeContact(double cutoff,
				vector<vector<double> >* position){
  double dx, dy, dz;
  for (int i {}; i < size; i++){
    set(i, i, get(i, i) + 1.0);
    for (int j {}; j < i; j++){
      dx = (*position)[i][0] - (*position)[j][0];
      dy = (*position)[i][1] - (*position)[j][1];
      dz = (*position)[i][2] - (*position)[j][2];
      if (inContact(distance(dx,dy,dz),cutoff)){
	set(i, j, get(i, j) + 1.0);
	set(j, i, get(j, i) + 1.0);
      }
    }
  }
}


// Compute colour contact
void ContactMap::computeColourContact(double cutoff,
				      vector<vector<double> >* position,
				      vector<int>* type){
  double dx, dy, dz;
  for (int i {}; i < size; i++){
    set(i, i, get(i, i) + 1.0);
    for (int j {}; j < i; j++){
      dx = (*position)[i][0] - (*position)[j][0];
      dy = (*position)[i][1] - (*position)[j][1];
      dz = (*position)[i][2] - (*position)[j][2];
      if (inContact(distance(dx,dy,dz),cutoff)){
	set(i, j, get(i, j) + 1.0);
	set(j, i, get(j, i) + 1.0);
      }
    }
  }
}


// Determine if in contact
bool inContact(double distance, double cutoff){
  if (distance < cutoff) return true;
  return false;
}

// Compute the length of a 3D vector
double distance(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}

// Accessor methods
void ContactMap::set(int i, int j, double value){
  contact->at(i,j) = value;
  //(*contact)[i][j] = value;
}

double ContactMap::get(int i, int j){
  return contact->at(i,j);
  //return (*contact)[i][j];
}

// Set all contacts to zero
void ContactMap::setZero(){
  contact->zeros();
  /*for (int i {}; i < size; i++){
    std::fill((*contact)[i].begin(), (*contact)[i].end(), 0.0);
    }*/
}

void ContactMap::reset(int n){
  // Only reset the contact map if the dimensions are valid
  if (n < 0) return;
  
  contact->zeros(n, n);
  /*
  // Erase the current data
  setZero();

  // Resize the contact map to fit new data
  if (size != n){
    (*contact).resize(n);
    for (int i {}; i < n; i++){
      (*contact)[i].resize(n, 0.0);
    }
  }
  */
  size = n;
}

// Normalisation methods
void ContactMap::vanillaNorm(){
  // Sum rows and columns
  vector<double>* xSum {new vector<double>(size, 0.0)};
  vector<double>* ySum {new vector<double>(size, 0.0)};
  for (int i {}; i < size; i++){
    for (int j {}; j < size; j++){
      (*xSum)[i] += get(i, j);
      (*ySum)[j] += get(i, j);
    }
  }
  // Normalise the contact map
  for (int i {}; i < size; i++){
    for (int j {}; j < size; j++){
      if ((*xSum)[i] != 0 && (*ySum)[j] != 0){
	set(i, j, get(i, j) / sqrt((*xSum)[i]*(*ySum)[j]));
      } else {
	set(i, j, 0.0);
      }
    }
  }
  delete xSum;
  delete ySum;
}

// Normalise the contact map by the contact probability function
void ContactMap::correlationNorm(){
  shared_ptr<vector<double> > prob {getGenomeDistContactProb()};
  for (int i {}; i < size; i++){
    for (int j {}; j < size; j++){
      set(i, j, get(i, j) / (*prob)[abs(i-j)]);
    }
  }
}

// Return the contact probability as a function of genome distance
shared_ptr<vector<double> > ContactMap::getGenomeDistContactProb(){
  shared_ptr<vector<double> > prob = make_shared<vector<double> >(size, 0.0);
  for (int i {}; i < size; i++){
    for (int j {}; j <= i; j++){
      (*prob)[abs(i-j)] += get(i,j);
    }
  }
  // Normalise
  for (int i {}; i < size; i++){
    (*prob)[i] /= static_cast<double>(size-i);
  }
  return prob;
}


// Reduce the resolution of the contact map
void ContactMap::reduceByBin(int bin){
  // No need to resize if there is no change
  if (bin <= 1) return;
  
  // Compute the new size of the contact map
  int n = static_cast<int>(ceil(static_cast<double>(size) / bin));

  // Average the original map values
  double sum;
  int count;
  for (int i {}; i < n; i++){
    for (int j {}; j < n; j++){
      sum = 0.0;
      count = 0;
      for (int k {i*bin}; k < (i+1)*bin && k < size; k++){
	for (int l {j*bin}; l < (j+1)*bin && l < size; l++){
	  sum += get(k, l);
	  count++;
	}
      }
      set(i, j, sum / static_cast<double>(count));
    }
  }
  
  // Resize the contact map to the new size
  contact->resize(n,n);
  /*  for (int i {}; i < n; i++){
    (*contact)[i].resize(n);
  }
  (*contact).resize(n);
  */
  size = n;
}

// Output method
void ContactMap::exportToFile(bool full, bool dense, bool space, string file){
  ofstream writer;
  writer.open(file);
  
  // Check that the file exists and can be written
  if (!writer){
    cout << "Problem with opening the output file!" << endl;
    return; 
  }
  
  writer << std::setprecision(10) << std::fixed;

  double value;
  const double tol {1e-15};
  for (int i {}; i < size; i++){
    for (int j {}; j < size; j++){
      if (!full && j > i) break;
      value = get(i, j);
      if (!dense && fabs(value) < tol)	continue;
      writer << i << " " << j << " " << value << endl;
    }
    if (!space) continue;
    writer << endl;
  }
}
