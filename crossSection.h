#pragma once
// Eventually part of runcard:
#define NPARTICLES 7

// Prototypes
int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *pdf); 
