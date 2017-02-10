#include <cuba.h>
#include <math.h>
#include <iostream>
#include "phaseSpace.h" 
#include "MomentumSet.h" 

using namespace std;

int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
   int n = 6;
   double roots = 8000.0;
   double s     = pow(roots, 2);

   // Generate phase space point from cuba x[]
   MomentumSet pset = phaseSpace(n, s, x);
   if (pset.ifail) {
      f[0] = 0.0;
      return 0;
   }
   
   // Check whether the point goes through all cuts (Check the momentum is not null, if it is then return)
   
   // Compute Matrix Element from it
   
   // Compute the PDFs
   
   // Compute flux factor for this ps point
   double flux = 1.0/pset.s(1,2);
   
   // Put everything together
   if (UNIT_PHASE) {
      f[0] = pset.weight;
   } else {
      f[0] = pset.x1 * pset.x2;
   }
   
   return 0;
}
