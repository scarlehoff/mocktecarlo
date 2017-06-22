#pragma once
#include <cuba.h>

double matrixElement(MomentumSet *pset);

// 2 -> 4
double C0g0WFH(MomentumSet *pset);

// 2 -> 5
double C1g0WFH(MomentumSet *pset);
double C1g0WFHup(MomentumSet *pset);
double C1g0WFHdown(MomentumSet *pset);
double C1g0WFHs0(int i1, int i5, int i3, int i2, int i4, MomentumSet *pset);

// 2 -> 6
double C2g0WFH(MomentumSet *pset);
double C2g0WFHnadj(int i1, int i5, int i4, int i2, int i6, int i3, MomentumSet *pset);
double C2g0WFHadj(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset);
cplx C2g0WFHadj1p2p(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset);
cplx C2g0WFHadj1p2m(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset);
cplx C2g0WFHadj1m2p(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset);
cplx C2g0WFHadj1m2m(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset);

// Virtual functions
double vertexCorrectionC0g1(MomentumSet *, double);
double integratedDipolesC0g1(MomentumSet *, double, const int, const double, const double);

// Common
double propagatorVBF(const double s1, const double s2, const int iboson);
double propagatorWFH(const double s1, const double s2);
double propagator(const double s, const double mass, const double width);
double helperC1(double x);
