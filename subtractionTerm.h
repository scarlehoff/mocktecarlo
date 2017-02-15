#pragma once
#include "MomentumSet.h"

double subtractionTerm(MomentumSet *pset, const double ptcut, const double rkt, const int minjets);

// Sub terms
double C1g0WFHS(MomentumSet *pset, const double ptcut, const double rkt, const int minjets);
double C2g0WFHS(MomentumSet *pset, const double ptcut, const double rkt, const int minjets);

// Antennae
double A30(MomentumSet *pset, const int i1, const int i2, const int i3);
