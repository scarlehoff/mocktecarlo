#pragma once

double matrixElement(MomentumSet *pset);

// 2 -> 4
double C0g0WFH(MomentumSet *pset);

// 2 -> 5
double C1g0WFH(MomentumSet *pset);
double C1g0WFHs0(int i1, int i5, int i3, int i2, int i4, MomentumSet *pset);

// 2 -> 6
double C2g0WFH(MomentumSet *pset);
double C2g0WFHnadj(int i1, int i5, int i4, int i3, int i6, int i2, MomentumSet *pset);
double C2g0WFHadj(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset);

// Common
double propagatorVBF(const double s1, const double s2, const int iboson);
double propagatorWFH(const double s1, const double s2);
double propagator(const double s, const double mass, const double width);
