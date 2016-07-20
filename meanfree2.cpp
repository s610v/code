#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <random>
using namespace std;

#define Nx 10E22 // number of soft photons
#define Ny 10000 // number of gamma ray photons
#define E 1./511. // energy of soft photons
#define Ey 2045 // energy of gamma ray
#define mu 0.5
#define PI 3.14159265

int main(){
	srand(time(NULL));
	double l[Ny];
	double lmean;
	double w = sqrt(E*Ey*(1-mu)/2.);
    double sigma = PI*(7.9408E-26)/pow(w,6)*((2*pow(w,4)+2*w*w-1)*log(w+sqrt(w*w-1))-w*(w*w-1)*sqrt(w*w-1)); // Cross section
    double Sigma = Nx*sigma; // Total cross section
    ofstream f;
    f.open("meanfreepathlength.txt");
    for (int i=0; i<Ny; i++) { // loop over number of particles
        double tau = (double) rand()/RAND_MAX;
        l[i] = -1./Sigma*log(tau); // transformation
        f << l[i] << endl;
    }
    cout << "Success" << endl;
    f.close();
}
