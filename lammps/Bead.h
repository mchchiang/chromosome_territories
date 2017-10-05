// Bead.h

#ifndef BEAD_H
#define BEAD_H

#include <map>
#include <memory>

using std::map;
using std::shared_ptr;

class Bead : std::enable_shared_from_this<Bead> {
  
private:
  double position [3];
  double velocity [3];
  int boundaryCount [3];
  int type;
  int label;

public:
  // Inner classes - For storing bonds and angles
  class Bond : std::enable_shared_from_this<Bond> {
	
  private:
	const int numOfBeads {2};
	int type;
	shared_ptr<Bead> beads[2];

  public:
	// Constructor
	Bond(int type, shared_ptr<Bead> bead1, shared_ptr<Bead> bead2);

	// Accessor methods
	shared_ptr<Bead> getBead(int id);
	void setType(int type);
	int getType();
	
	void unbond();
  };

  class Angle : std::enable_shared_from_this<Angle> {
	
  private:
	const int numOfBeads {3};
	int type;
	shared_ptr<Bead> beads[3];
	
  public:
	// Constructor
	Angle(int type, shared_ptr<Bead> bead1, 
		  shared_ptr<Bead> bead2, shared_ptr<Bead> bead3);
	
	// Accessor methods
	shared_ptr<Bead> getBead(int id);
	void setType(int type);
	int getType();
	
	void unbond();
  };

  // Constructors
  Bead(double x, double y, double z,
       double vx, double vy, double vz,
       int nx, int ny, int nz, int type, int label);
  Bead(double x, double y, double z);
  Bead();

  // Accesor methods
  void setPosition(int dim, double value);
  double getPosition(int dim);
  void setVelocity(int dim, double value);
  double getVelocity(int dim);
  void setBoundaryCount(int dim, int count);
  int getBoundaryCount(int dim);
  void setType(int type);
  int getType();
  void setLabel(int label);
  int getLabel();

  // For modifying bonds and angles
  void addBond(shared_ptr<Bond> bond);
  void addBondWith(int type, shared_ptr<Bead> bead);
  void removeBond(shared_ptr<Bond> bond);
  void removeBondWith(shared_ptr<Bead> bead);
  void addAngle(shared_ptr<Angle> angle);
  void addAngleWith(int type, shared_ptr<Bead> bead1, shared_ptr<Bead> bead2);
  void removeAngle(shared_ptr<Angle> angle);
  void removeAngleWith(shared_ptr<Bead> bead1, shared_ptr<Bead> bead2);
  void removeAllBonds();
  void removeAllAngles();

private:
  map< shared_ptr<Bond>, shared_ptr<Bond> > bondList {};
  map< shared_ptr<Angle>, shared_ptr<Angle> > angleList {};

};

#endif
