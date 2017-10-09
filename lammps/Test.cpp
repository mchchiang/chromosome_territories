// Test.cpp

#include <iostream>
#include <memory>
#include "LAMMPS.h"
#include "Polymer.h"
#include "Bead.h"

using std::cout;
using std::endl;
using std::shared_ptr;
using std::make_shared;

int main(){
  double lx {100.0};
  double ly {100.0};
  double lz {100.0};
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>(lx, ly, lz);

  lammps->setTypesOfBeads(1);
  lammps->setTypesOfBonds(1);
  lammps->setTypesOfAngles(1);
  
  shared_ptr<Polymer> p1 = lammps->createRandomWalkPolymer(0,1000);
  shared_ptr<Polymer> p2 = lammps->createRandomWalkPolymer(1,1000);
  shared_ptr<Bead> b1 = lammps->createBead(0);
  shared_ptr<Bead> b2 = lammps->getPolymer(1)->getBead(324);
  shared_ptr<Bead>  b3 = lammps->getPolymer(1)->getBead(483);
  b1->addBondWith(0, b2);
  b1->addAngleWith(0, b2, b3);
  b2->removeBondWith(b1);
  
  lammps->getPolymer(0)->getBead(0)->setPosition(0,0.23498);
  lammps->getPolymer(0)->getBead(0)->setPosition(1,2.6574);
  lammps->getPolymer(0)->getBead(0)->setPosition(2,7.13);
  lammps->exportData("mydata.dat", "map.dat");
  
  lammps->importData("mydata.dat", "map.dat");
  cout << lammps->getPolymer(0)->getBead(0)->getPosition(2) << endl;
  cout << lammps->getPolymer(0)->getBead(0)->getPosition(1) << endl;
  cout << lammps->getPolymer(0)->getBead(0)->getPosition(0) << endl;
}
