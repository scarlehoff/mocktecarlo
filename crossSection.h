#pragma once
// Eventually part of runcard:
#define NPARTICLES 7
#define MINJETS_   2

// Prototypes
int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *pdf); 
