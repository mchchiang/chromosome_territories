// ContactMapLib.h
#ifndef CONTACTMAP_HPP
#define CONTACTMAP_HPP

#include <memory>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::shared_ptr;

class ContactMap;
using CMap = shared_ptr<ContactMap>;

// Determine if the distance is less than the cutoff
bool inContact(double distance, double cutoff);

// Compute the length of a 3D vector
double distance(double x, double y, double z);

class ContactMap {

private:
  vector< vector<double>>* contact;
  int nx;
  int ny;
  
  // Private constructor
  ContactMap(int nx, int ny, int value = 0.0);
  
  // Helper methods for computing contacts
  void computeContact(double cutoff, vector< vector<double> >* position);
  void computeColourContact(double cutoff, vector< vector<double> >* position, 
			   vector<int>* type);

public:
  // Static factory methods
  // Create a contact map with zero interaction
  static CMap createZeroMap(int nx, int ny);

  // Create a contact map based on a 2D matrix
  static CMap createFromArray(vector< vector<double>>* matrix);

  // Create a contact map from a position file
  static CMap 
  createFromPosFile(int numOfBeads, double lx, double ly, double lz,
		    double cutoff, string contactType,
		    int startTime, int endTime, int timeInc, string file);


  // Create a contact map from a matrix file
  static CMap createFromMatrixFile(int nx, int ny, bool full, string file);
  
  // Accessor methods
  // Set the contact probability of the interaction pair (i,j)
  void set(int i, int j, double value);

  // Get the contact probability of the interaction pair (i,j)
  double get(int i, int j);

  // Map manipulation methods
  // Normalise the contact map using the vanilla normalisation method
  void vanillaNorm();

  // Normalise the contact map by the contact probability function
  void correlationNorm();

  // Return the contact probability as a function of genome distance
  shared_ptr<vector<double> > getGenomeDistContactProb();
  
  // Reduce the resolution of the contact map
  void reduceByBin(int binx, int biny);

  // Set all contacts to zero
  void setZero();

  // Set the contact map to the defined size with all contacts being zero
  void reset(int nx, int ny);

  // Import contact map from position file
  void importFromPosFile(int numOfBeads, double lx, double ly, double lz,
			 double cutoff, string contactType,
			 int startTime, int endTime, int timeInc, string file);

  // Import contact map from a matrix file
  void importFromMatrixFile(int nx, int ny, bool full, string file);

  void importFromArray(vector< vector<double>>* matrix);

  // Output contact map to a file
  void exportToFile(bool full, bool dense, bool space, string file);
};

#endif
