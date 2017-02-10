#include <cuba.h>
#include <math.h>
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include "phaseSpace.h" 
#include "MomentumSet.h" 
#include "matrixElement.h"

#define FBGEV2 389379365600
#define NC 3.0
#define AMZ 7.5563839072873536E-003

using namespace std;
using namespace LHAPDF;

int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *pdf) {
   int n = 6;
   double roots = 8000.0;
   double muR   = 125.0 ; // muR = muF = mH
   double s     = pow(roots, 2);

   // Generate phase space point from cuba x[]
   MomentumSet pset = phaseSpace(n, s, x);
   if (pset.ifail) {
      f[0] = 0.0;
      return 0;
   }
   
   // Check whether the point goes through all cuts (Check the momentum is not null, if it is then return)
   
   // Compute Matrix Element from it
   double mesq = matrixElement(&pset);
   
   // Compute the PDFs
   double f1      = (*((PDF **) pdf))->xfxQ2(2, pset.x1, pow(muR,2));
               f1+= (*((PDF **) pdf))->xfxQ2(4, pset.x1, pow(muR,2));
   double f2      = (*((PDF **) pdf))->xfxQ2(1, pset.x2, pow(muR,2));
               f2+= (*((PDF **) pdf))->xfxQ2(3, pset.x2, pow(muR,2));
   double pdfval  = f1*f2/pset.x1/pset.x2;
   double alpha_s = (*((PDF **) pdf))->alphasQ2(pow(muR,2));
   
   // Compute flux factor for this ps point and QCD factor (cte)
   pset.weight = pset.weight/pset.s(1,2);
   double average = (1.0/NC)*(1.0/NC)/4.0;
   double qcdborn = pow(4.0*M_PI*AMZ, 3)*NC*NC/2.0;
   double qcdfactor = 1.0;
   switch(n) {
      case 6:
         qcdfactor = qcdborn;
         break;
      case 7:
         qcdfactor = qcdborn*NC*(4.0*M_PI*alpha_s);
         break;
   }
   double flux = average*qcdfactor*FBGEV2/2.0;

   // Put everything together
   if (UNIT_PHASE) {
      f[0] = pset.weight;
   } else {
      f[0] = mesq*flux*pset.weight*pdfval;
   }
   
   return 0;
}
