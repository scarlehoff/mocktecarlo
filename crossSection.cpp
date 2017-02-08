#include "Cuba-4.2/include/cuba.h"

int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
   double fun;
   fun = x[0] * x[1]*4.0;

   // Generate phase space point from cuba x[]
   
   // Check whether the point goes through all cuts
   
   // Compute Matrix Element from it
   
   // Compute the PDFs
   
   // Put everything together
   
   f[0] = fun;
   return 0;
}
