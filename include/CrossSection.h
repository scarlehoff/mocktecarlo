#pragma once

#include "Runcard.h"

// Prototypes
int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *pdf); 
double pdfValue(double x1, double x2, double muR2, void *pdf);
