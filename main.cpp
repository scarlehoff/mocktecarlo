#include <iostream>
#include <fstream>
#include <cuba.h>
#include <LHAPDF/LHAPDF.h>

using namespace std;
using namespace LHAPDF;

#include "crossSection.h"


int testRecipe(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata);

int main() {
   // ----- Vegas variables
   int userdata = 0;
   int ncomp = 1;
   int nvec = 1;
   int spin = 0;
   int neval, fail;
   int verbose = 2;
   char statefile[1] = "";
   // Dimensionality
   int ndim = 9; // base, n = 6
   ndim += (NPARTICLES-6)*3; // 7: 12, 8: 15
   // Error Tolerance
   double epsrel = 1e-7;
   double epsabs = 1e-7;
   // Integration parameters
   int seed = 0;
   int mineval = 0;
   // cluster
   int maxeval   =100000000;
   int nstart    = 5000000;
   int nincrease = 2000000;
   int nbatch    = 600000;
   // desktop
//   int maxeval   = 5000000;
//   int nstart    = 200000;
//   int nincrease = 800000;
//   int nbatch    = 200000;
   // --- 
   int gridno = 0;
   cubareal integral[ncomp], error[ncomp], prob[ncomp];

   // PDF for alpha_s
   const string setname = "MSTW2008nnlo90cl";
   const int imem = 0;
   const PDF* pdf = mkPDF(setname, imem);
   double alpha_s = pdf->alphasQ2(pow(125.0,2)); // necessary for LHAPDF to work multicore?

   Vegas(ndim, ncomp,
         crossSection, &pdf, nvec,
//         testRecipe, &userdata, nvec,
         epsrel, epsabs,
         verbose, seed,
         mineval, maxeval,
         nstart, nincrease, nbatch,
         gridno, statefile, &spin,
         &neval, &fail,
         integral, error, prob);

   // fbtoGeV

   return 0;
}

#include "FourVector.h"
#include "phaseSpace.h"

int testRecipe(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
   cout << "Debug system, press enter to continue" << endl;
   cin.ignore();
   cout << "Testing four vectors" << endl;
//   FourVector test = FourVector(x[0], x[1], x[2], x[3]);
//   FourMomentum test = FourMomentum(x[0], x[1], x[2], x[3]);
//   cout << test << endl;
//   cout << test.sq() << endl;
//   cout << 4.0*test << endl;
//   cout << test*4.0 << endl;
//
//   vector <FourMomentum> test2 ;
//   test2.reserve(2);
//   test2.emplace_back(x[0],x[1],x[2],x[3]);
//   test2.emplace_back(x[3],x[1],x[2],x[3]);
//   cout << test2[0] << endl;
//   cout << test2[0].sq() << endl;
//   cout << test2[1] << endl;
//   cout << test2[1].sq() << endl;

   cout << "Phase Space test" << endl;
   MomentumSet pset = phaseSpace(6, 1000.0, x);
}
