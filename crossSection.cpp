#include <cuba.h>
#include <math.h>
#include "phaseSpace.h" 
//#include "MomentumSet.h" included in phaseSpace.h

int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
   double fun;
   fun = x[0] * x[1]*4.0;

   int n = 6;
   double roots = 8000.0;
   double s     = pow(roots, 2);

   // Generate phase space point from cuba x[]
   MomentumSet pset = phaseSpace(n, s, x);
   
   // Check whether the point goes through all cuts
   
   // Compute Matrix Element from it
   
   // Compute the PDFs
   
   // Compute flux factor for this ps point
   double flux = 1.0/pset.shat;
   
   // Put everything together
   
   f[0] = fun;
   return 0;
}
