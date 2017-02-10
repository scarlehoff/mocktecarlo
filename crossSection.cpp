#include <cuba.h>
#include <math.h>
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include "phaseSpace.h" 
#include "MomentumSet.h" 

#define FBGEV2 389379365600

using namespace std;
using namespace LHAPDF;

int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *pdf) {
   int n = 6;
   double roots = 8000.0;
   double s     = pow(roots, 2);

   // Generate phase space point from cuba x[]
   MomentumSet pset = phaseSpace(n, s, x);
   if (pset.ifail) {
      f[0] = 0.0;
      return 0;
   }
   double f1 = (*((PDF **) pdf))->xfxQ2(2, pset.x1, pow(125.0,2));
   f1 += (*((PDF **) pdf))->xfxQ2(4, pset.x1, pow(125.0,2));
   double f2 = (*((PDF **) pdf))->xfxQ2(1, pset.x2, pow(125.0,2));
   f2 += (*((PDF **) pdf))->xfxQ2(3, pset.x2, pow(125.0,2));
   double pdfval = f1*f2/pset.x1/pset.x2;
   
   // Check whether the point goes through all cuts (Check the momentum is not null, if it is then return)
   
   // Compute Matrix Element from it
   double mesq = 1.0;
   
   // Compute the PDFs
   
   // Compute flux factor for this ps point
   double flux = FBGEV2/pset.s(1,2)/2.0;

   // Put everything together
   if (UNIT_PHASE) {
      f[0] = pset.weight;
   } else {
      f[0] = mesq*flux*pset.weight*pdfval;
   }
   
   return 0;
}
