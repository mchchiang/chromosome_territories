// Bead.h

#ifndef BEAD_H
#define BEAD_H

#include <map>
#include <memory>
#include "Bond.h"
#include "Angle.h"

using std::map;
using std::shared_ptr;

class Bead {
  
private:
  const int dimension {3};
  double position [dimension] {};
  double velocity [dimension] {};
  int boundaryCount [dimension] {};
  int type {};
  int label {};
  map< shared_ptr<Bond>, shared_ptr<Bond> > bondList {};
  map< shared_ptr<Angle>, shared_ptr<Angle> > angleList {};

public:
  // Constructors
  Bead(double x, double y, double z,
       double vx, double vy, double vz,
       int nx, int ny, int nz, int type, int label);
  Bead(double x, double y, double z);
  Bead(const Bead& b);
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
  void removeBond(shared_ptr<Bond> bond);
  void addAngle(shared_ptr<Angle> angle);
  void removeAngle(shared_ptr<Angle> angle);
  void removeAllBonds();
  void removeAllAngles();

};

#endif
