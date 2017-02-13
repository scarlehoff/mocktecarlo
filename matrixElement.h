#pragma once

double matrixElement(MomentumSet *pset);
double C0g0WFH(MomentumSet *pset);
double C1g0WFH(MomentumSet *pset);
double C1g0WFHs0(int i1, int i5, int i3, int i2, int i4, MomentumSet *pset);
double propagatorVBF(const double s1, const double s2, const int iboson);
double propagatorWFH(const double s1, const double s2);
double propagator(const double s, const double mass, const double width);
