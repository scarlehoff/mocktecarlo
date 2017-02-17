#pragma once
#include "MomentumSet.h"

double subtractionTerm(MomentumSet *pset, const double ptcut, const double rkt, const int minjets);

// Sub terms
double C1g0WFHS(MomentumSet *pset, const double ptcut, const double rkt, const int minjets);
double C2g0WFHS(MomentumSet *pset, const double ptcut, const double rkt, const int minjets);
// wrappers
double adjLimitIF(MomentumSet *pset, int i1, int i2, int i3, double (*MatrixE)(MomentumSet*), const double ptcut, const double rkt, const int minjets);
double adjLimitFF(MomentumSet *pset, int i1, int i2, int i3, double (*MatrixE)(MomentumSet*), const double ptcut, const double rkt, const int minjets);
double nadjLimitIF(MomentumSet *pset, int i1, int i2, int i3, double (*MatrixE)(MomentumSet*), const double ptcut, const double rkt, const int minjets);
// Antennae
double A30(MomentumSet *pset, const int i1, const int i2, const int i3);
double d30(MomentumSet *pset, const int i1, const int i2, const int i3);
