// GenreateSphericalPoints

#include<iomanip>
#include<fstream>
#include<cmath>

using namespace std;

int main(int argc, char * argv[])
{

  int numOfPoints = stoi(argv[1], NULL);
  double maxRadius = stod(argv[2], NULL);
  
  srand(time(NULL)); // Initialise the random generator
  
  // Generate each bead's position randomly within the laminar
  double pi {M_PI};
  
  double u {}, v {}, w {};
  double r {}, theta {}, phi {};
  double x {}, y {}, z {};

  ofstream unbiased, biased;
  unbiased.open("spherical_points_unbiased.dat");
  biased.open("spherical_points_biased.dat");

  for (int i {}; i < numOfPoints; i++)
	{
	  u = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	  v = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	  w = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
	  r = maxRadius * pow(u, 1.0/3.0);
	  theta = acos(1 - 2 * v);
	  phi = 2 * pi * w;
	  x = r * sin(theta) * cos(phi);
	  y = r * sin(theta) * sin(phi);
	  z = r * cos(theta);
	  unbiased << fixed << setprecision(16);
	  unbiased << x << " " << y << " " << z << endl;
	  r = u * maxRadius;
	  theta = v * pi;
	  phi = w * 2 * pi;
	  x = r * sin(theta) * cos(phi);
	  y = r * sin(theta) * sin(phi);
	  z = r * cos(theta);
	  biased << fixed << setprecision(16);
	  biased << x << " " << y << " " << z << endl;
	}

  unbiased.close();
  biased.close();
}
