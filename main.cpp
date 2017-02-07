#include <iostream>
#include <fstream>
#include <cuba.h>
#include <LHAPDF/LHAPDF.h>


using namespace std;
using namespace LHAPDF;

#include "crossSection.h"

int main() {
   // "User defined input"
   double mu = 125.0;
   // ----- Vegas variables
   int ncomp = 1;
   int userdata = 0;
   int nvec = 1;
   int spin = 0;
   int neval, fail;
   int verbose = 2;
   char statefile[1] = "";
   // Dimensionality
   int ndim = 2;
   // Error Tolerance
   double epsrel = 1e-3;
   double epsabs = 1e-4;
   // Integration parameters
   int seed = 0;
   int mineval = 0;
   int maxeval = 500000;
   int nstart = 10000;
   int nincrease = 5000;
   int nbatch = 10000;
   int gridno = 0;
   cubareal integral[ncomp], error[ncomp], prob[ncomp];

   // PDF for alpha_s
//   const string setname = "MSTW2008nnlo90cl";
//   const int imem = 0;
//   const PDF* pdf = mkPDF(setname, imem);
//   double alpha_s = pdf->alphasQ2(pow(mu,2));

   Vegas(ndim, ncomp,
         crossSection, &userdata, nvec,
         epsrel, epsabs,
         verbose, seed,
         mineval, maxeval,
         nstart, nincrease, nbatch,
         gridno, statefile, &spin,
         &neval, &fail,
         integral, error, prob);

   return 0;
}

