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
#include "PositionReader.hpp"

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
  PositionReader reader;
  reader.open(file, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    return;
  }

  // Reset the contact map
  reset(numOfBeads);

  int count {};
  long time {};
  while(reader.nextFrame()){
    time = reader.getTime();
    if (time >= startTime && time <= endTime){
      if (contactType == "colour"){
	computeColourContact(cutoff, reader);
      } else {
	computeContact(cutoff, reader);
      }
      count++;
    } else if (time > endTime) {
      break;
    }
  }
  reader.close();
  // Normalise contact map
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j < numOfBeads; j++){
      set(i, j, get(i, j) / count);
    }
  }
}



/*void ContactMap::importFromPosFile(int numOfBeads, 
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
  }*/

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
void ContactMap::computeContact(double cutoff, const PositionReader& reader){
  cout << "Computing contacts" << endl;
  double dx, dy, dz;
  for (int i {}; i < size; i++){
    set(i, i, get(i, i) + 1.0);
    for (int j {}; j < i; j++){
      dx = reader.getUnwrappedPosition(i,0) - reader.getUnwrappedPosition(j,0);
      dy = reader.getUnwrappedPosition(i,1) - reader.getUnwrappedPosition(j,1);
      dz = reader.getUnwrappedPosition(i,2) - reader.getUnwrappedPosition(j,2);
      if (inContact(distance(dx,dy,dz),cutoff)){
	set(i, j, get(i, j) + 1.0);
	set(j, i, get(j, i) + 1.0);
      }
    }
  }
}


// Compute colour contact
void ContactMap::computeColourContact(double cutoff, const PositionReader& reader){
  double dx, dy, dz;
  for (int i {}; i < size; i++){
    set(i, i, get(i, i) + 1.0);
    for (int j {}; j < i; j++){
      dx = reader.getUnwrappedPosition(i,0) - reader.getUnwrappedPosition(j,0);
      dy = reader.getUnwrappedPosition(i,1) - reader.getUnwrappedPosition(j,1);
      dz = reader.getUnwrappedPosition(i,2) - reader.getUnwrappedPosition(j,2);
      if (inContact(distance(dx,dy,dz),cutoff)){
	set(i, j, get(i, j) + 1.0);
	set(j, i, get(j, i) + 1.0);
      }
    }
  }
}

/*
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
*/

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
}

double ContactMap::get(int i, int j){
  return contact->at(i,j);
}

// Return the linear size of the contact map
int ContactMap::getSize(){
  return size;
}

// Set all contacts to zero
void ContactMap::setZero(){
  contact->zeros();
}

void ContactMap::reset(int n){
  // Only reset the contact map if the dimensions are valid
  if (n < 0) return;
  
  contact->zeros(n, n);
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
void ContactMap::linearProbNorm(){
  shared_ptr<vector<double> > prob {getLinearProb()};
  double value;
  for (int i {}; i < size; i++){
    for (int j {}; j < size; j++){
      value = (*prob)[abs(i-j)];
      if (fabs(value) < 1e-10){
	set(i, j, 0.0);
      } else {
	value = get(i, j) / value;
	set(i, j, value);
      }
    }
  }
}

// Convert the contact map to a Pearon's correlation matrix
void ContactMap::convertToCorrelation(){
  // Correlation matrix is symmetric
  (*contact) = cor((*contact));
  // Check for division by zero in computing correlation -
  // set those entries to zero
  for (int i {}; i < size; i++){
    for (int j {}; j < size; j++){
      if (isnan(contact->at(i,j))){
	contact->at(i,j) = 0.0;
      }
    }
  }
}

// Return the contact probability as a function of genome distance
shared_ptr<vector<double> > ContactMap::getLinearProb(){
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

double ContactMap::maxEigen(double conv, shared_ptr<vector<double> > evec){
  // Do power iteration to find the largest eigenvalue and eigenvector
  // Init random vector
  vec* v {new vec(size, fill::randu)};
  vec* w {new vec(size)};

  (*v) = (*v) / norm(*v); // Make sure the vector is normalised

  cout << "Finding maximum eigenvalue ..." << endl;
  int iter {};
  double delta;
  do {
    (*w) = (*contact) * (*v);
    (*w) = (*w) / norm(*w);
    delta = norm((*w)-(*v));
    (*v) = (*w);
    iter++;
    cout << "After iteration " << iter << ": delta = " << delta << endl;
  } while (delta > conv);

  // Store the eigenvector
  (*evec) = conv_to<vector<double> >::from(*v);

  // Get the eigenvalue 
  double vnorm {norm(*v)};
  double eigenval {as_scalar((*v).t()*(*contact)*(*v)) / (vnorm*vnorm)};
  delete v;
  delete w;
  return eigenval;
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

void ContactMap::exportCombineMapsToFile(CMap map1, CMap map2, 
					 bool dense, bool space, string file){
  
  // Check that both contact matrices are the same size
  const int size {map1->getSize()};
  if (size != map2->getSize()){
    cout << "Contact maps must be the same size!" << endl;
    return;
  }

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
      if (i <= j){
	value = map2->get(i,j);
      } else {
	value = map1->get(i,j);
      }
      if (!dense && fabs(value) < tol)	continue;
      writer << i << " " << j << " " << value << endl;
    }
    if (!space) continue;
    writer << endl;
  }
}
