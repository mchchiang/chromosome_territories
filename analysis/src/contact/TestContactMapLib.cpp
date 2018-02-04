/* TestContactMapLib.cpp
 *
 * A set of unit tests for testing the methods in the
 * ContactMapLib class
 */

#define BOOST_TEST_MODULE TestContactMapLib
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "ContactMapLib.hpp"

using std::cout;
using std::endl;

const double tol {1e-5};

BOOST_AUTO_TEST_CASE(TestDistance){
  double x, y, z, expect;

  // Simple case
  x = 0.0;
  y = 0.0;
  z = 1.0;
  expect = 1.0;
  BOOST_CHECK_CLOSE(distance(x,y,z), expect, tol);
  
  x = 3.0;
  y = 4.0;
  z = 0.0;
  expect = 5.0;
  BOOST_CHECK_CLOSE(distance(x,y,z), expect, tol);

  // Non trivial case
  x = 4.5;
  y = 2.6;
  z = 8.7;
  expect = 10.13410;
  BOOST_CHECK_CLOSE(distance(x,y,z), expect, tol);

}

BOOST_AUTO_TEST_CASE(TestInContact){
  double distance, cutoff;
  distance = 5.1;
  cutoff = 8.9;
  BOOST_CHECK_EQUAL(inContact(distance, cutoff), true);
  
  distance = 4.8;
  cutoff = 2.3;
  BOOST_CHECK_EQUAL(inContact(distance, cutoff), false);
}

BOOST_AUTO_TEST_CASE(TestCreateZeroMap){
  int nx {3}, ny {4};
  CMap map = ContactMap::createZeroMap(nx, ny);
  for (int i {}; i < nx; i++){
    for (int j {}; j < ny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), 0.0, tol);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestCreateFromArray){
  int nx {3}, ny {4};
  vector< vector<double>>* matrix  = new vector< vector<double>>
    {
      {0.2, 0.3, 0.1, 0.4},
      {0.1, 0.2, 0.5, 0.3},
      {0.4, 0.1, 0.7, 0.6}
    };
  CMap map = ContactMap::createFromArray(matrix);
  for (int i {}; i < nx; i++){
    for (int j {}; j < ny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), (*matrix)[i][j], tol);
    }
  }
  delete matrix;
}

BOOST_AUTO_TEST_CASE(TestCreateFromPosFile){
  string file {"ContactMapPosFileTest1.dat"};
  int l {10};
  int numOfBeads {5};
  double cutoff {2.0};
  int startTime {0}, endTime {0}, timeInc {1};
  string contactType {"simple"};
  
  CMap map = ContactMap::createFromPosFile(numOfBeads, l, l, l, cutoff, 
					   contactType, startTime, endTime, 
					   timeInc, file);
  vector< vector<double>> expect
    {
      {1, 1, 0, 0, 0},
      {1, 1, 1, 1, 0},
      {0, 1, 1, 0, 1},
      {0, 1, 0, 1, 0},
      {0, 0, 1, 0, 1}
    };
  for (int i {}; i < numOfBeads; i++){
    for (int j {}; j < numOfBeads; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), expect[i][j], tol);
    }
  }
}


BOOST_AUTO_TEST_CASE(TestReset1){
  // Test the case for expanding the contact map
 int nx {3}, ny {4};
  vector< vector<double>>* matrix  = new vector< vector<double>>
    {
      {0.2, 0.3, 0.1, 0.4},
      {0.1, 0.2, 0.5, 0.3},
      {0.4, 0.1, 0.7, 0.6}
    };
  
  CMap map = ContactMap::createFromArray(matrix);
  int newnx {nx+10}, newny {ny+9};
  map->reset(newnx,newny);
  for (int i {}; i < newnx; i++){
    for (int j {}; j < newny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), 0.0, tol);
    }
  }
  delete matrix; 
}

BOOST_AUTO_TEST_CASE(TestReset2){
  // Test the case for reducing the contact map
 int nx {3}, ny {4};
  vector< vector<double>>* matrix  = new vector< vector<double>>
    {
      {0.2, 0.3, 0.1, 0.4},
      {0.1, 0.2, 0.5, 0.3},
      {0.4, 0.1, 0.7, 0.6}
    };
  
  CMap map = ContactMap::createFromArray(matrix);
  int newnx {nx-1}, newny {ny-2};
  map->reset(newnx,newny);
  for (int i {}; i < newnx; i++){
    for (int j {}; j < newny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), 0.0, tol);
    }
  }
  delete matrix; 
}


BOOST_AUTO_TEST_CASE(TestReduceByBin1){
  int nx {4}, ny {4};
  vector< vector<double>>* matrix  = new vector< vector<double>>
    {
      {0.2, 0.3, 0.1, 0.4},
      {0.1, 0.2, 0.5, 0.3},
      {0.4, 0.1, 0.7, 0.6},
      {0.5, 0.1, 0.9, 0.8}
    };

  vector< vector<double>> expect 
    {
      {0.2, 0.325},
      {0.275, 0.75}
    };
  
  CMap map = ContactMap::createFromArray(matrix);
  int binx {2}, biny {2};
  int newnx {nx/binx}, newny {ny/biny};
  map->reduceByBin(binx,biny);
  for (int i {}; i < newnx; i++){
    for (int j {}; j < newny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), expect[i][j], tol);
    }
  }
  delete matrix; 
}

BOOST_AUTO_TEST_CASE(TestReduceByBin2){
  vector< vector<double>>* matrix  = new vector< vector<double>>
    {
      {0.2, 0.3, 0.1},
      {0.1, 0.2, 0.5},
      {0.4, 0.1, 0.7},
      {0.5, 0.1, 0.9},
      {0.6, 0.8, 0.1}
    };

  vector< vector<double>> expect 
    {
      {0.2, 0.3},
      {0.275, 0.8},
      {0.7, 0.1}
    };
  
  CMap map = ContactMap::createFromArray(matrix);
  int binx {2}, biny {2};
  int newnx {3}, newny {2};
  map->reduceByBin(binx,biny);
  for (int i {}; i < newnx; i++){
    for (int j {}; j < newny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), expect[i][j], tol);
    }
  }
  delete matrix; 
}

BOOST_AUTO_TEST_CASE(TestVanillaNorm1){
  int nx {4}, ny {4};
  vector< vector<double>>* matrix  = new vector< vector<double>>
    {
      {0.2, 0.3, 0.1, 0.4},
      {0.1, 0.2, 0.5, 0.3},
      {0.4, 0.1, 0.7, 0.6},
      {0.5, 0.1, 0.9, 0.8}
    };

  vector< vector<double>> expect 
    {
      {0.18257419, 0.35856858, 0.06741999, 0.27602622},
      {0.08703883, 0.22792115, 0.32141217, 0.19738551},
      {0.27216553, 0.08908708, 0.35176324, 0.30860670},
      {0.30096463, 0.07881104, 0.40009880, 0.36401260}
    };
  
  CMap map = ContactMap::createFromArray(matrix);
  map->vanillaNorm();
  for (int i {}; i < nx; i++){
    for (int j {}; j < ny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), expect[i][j], tol);
    }
  }
  delete matrix;  
}

BOOST_AUTO_TEST_CASE(TestVanillaNorm2){
  int nx {5}, ny {3};
  vector< vector<double>>* matrix  = new vector< vector<double>>
    {
      {0.2, 0.3, 0.1},
      {0.1, 0.2, 0.5},
      {0.4, 0.1, 0.7},
      {0.5, 0.1, 0.9},
      {0.6, 0.8, 0.1}
    };

  vector< vector<double>> expect 
    {
      {0.19245009, 0.31622777, 0.08512565},
      {0.08333333, 0.18257419, 0.36860489},
      {0.27216553, 0.07453560, 0.42135049},
      {0.30429031, 0.06666667, 0.48454371},
      {0.36514837, 0.53333333, 0.05383819}
    };
  
  CMap map = ContactMap::createFromArray(matrix);
  map->vanillaNorm();
  for (int i {}; i < nx; i++){
    for (int j {}; j < ny; j++){
      BOOST_CHECK_CLOSE(map->get(i,j), expect[i][j], tol);
    }
  }
  delete matrix; 
}


