#include "Cuba-4.2/include/cuba.h"

int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata) {
   double fun;
   fun = x[0] * x[1]*4.0;
   f[0] = fun;
   return 0;
}
