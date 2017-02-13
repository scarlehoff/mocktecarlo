#pragma once
#include <math.h>
#include <cuba.h>
#include <stdlib.h>

#include "MomentumSet.h"
#include "FourVector.h"



#define UNIT_PHASE 1
#define DEBUG 0
#define COSTHMIN -1.0
#define COSTHMAX 1.0
#define PHIMAX 2.0*M_PI
#define PHIMIN 0.0

#define y0 1e-7


typedef boost::multi_array<double, 2> matrix;

MomentumSet phaseSpace(const int npar, const double s_input, const cubareal x[]);

// Sampling functions
double generateInwa(const int itype, const double r, const double shat, const double mh, double *wtps);
double pickRand(const int itype, const double r, const double smax, const double smin, double *wtps);

// Generate 4-Momenta
int p3generic_nj(const cubareal x[], const int n_i, const double shat, const double s1, const double s2, const double s3, std::vector<FourMomentum> *pset, double *wtps);
void makePs2cm_nj(const double s, const double cos12, const double phi, const double s1, const double s2, FourMomentum *p1, FourMomentum *p2);

// G-Determinant
void glimits(const double x, const double y, const double z, const double u, const double v, const double w, double *ymin, double *ymax);
double gres(const double x, const double y, const double z, const double u, const double v, const double w);


// Generic
double dlambda(const double s1, const double s2, const double s3);
