#include "MomentumSet.h"
#include <math.h>
#include <cuba.h>
#include <stdlib.h>

#define UNIT_PHASE 1.0
#define COSTHMIN -1.0
#define COSTHMAX 1.0
#define PHIMAX 2.0*M_PI
#define PHIMIN 0.0

#define y0 1e-7


typedef boost::multi_array<double, 2> matrix;

MomentumSet phaseSpace(const int npar, const double s_input, const cubareal x[]);

// Sampling functions
double generateInwa(const int itype, const double r, const double shat, double *wtps);
double pickRand(const int itype, const double r, const double smax, const double smin, double *wtps);

// Generate 4-Momenta
momentum_t *p3generic_nj(const cubareal x[], const double s_input[], double *wtps);
void makePs2cm_nj(const double s, const double cos12, const double phi, const double s1, const double s2, momentum_t *p1, momentum_t *p2);

// Generic
double dlambda(const double s1, const double s2, const double s3);
double dot(momentum_t *p1, momentum_t *p2);
matrix unboostrest(momentum_t *pin);
