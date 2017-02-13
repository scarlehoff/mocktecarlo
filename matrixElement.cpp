#include "MomentumSet.h"
#include "matrixElement.h"

double matrixElement(MomentumSet *pset) {
   switch (pset->npar) {
      case 6:
         return C0g0WFH(pset);
      case 7:
         return C1g0WFH(pset);
   }
}

// Amplitudes
// 2 -> 4
double C0g0WFH(MomentumSet *pset) {
   int i1 = 2; int i2 = 1; int i3 = 3; int i4 = 4;

   double s1i = pset->s(i1,i4);
   double s2j = pset->s(i2,i3);
   double prop = propagatorVBF(s1i, s2j, 1);

   // Amplitude
   cplx zamp = pset->zA(i1, i2)*pset->zA(i3,i4);

   double amp = pow(abs(zamp), 2);
   return 2.0*amp*prop;
}
// 2 -> 5
double C1g0WFH(MomentumSet *pset) {
   int i1 = 2; int i2 = 1; int i3 = 4; int i4 = 5; int i5 = 3;
   double a1 = C1g0WFHs0(i1, i5, i3, i2, i4, pset);
   double a2 = C1g0WFHs0(i2, i5, i4, i1, i3, pset);
   return a1 + a2;
}

double C1g0WFHs0(int i1, int i5, int i3, int i2, int i4, MomentumSet *pset) {
   double s1i = pset->s(i1,i4);
   double s1k = pset->s(i1,i5);
   double sik = pset->s(i4,i5);
   double s2j = pset->s(i2,i3);
   double q1  = s1i + s1k + sik;

   double prop = propagatorVBF(q1, s2j, 1);

   // Amplitude
   // g+
   cplx z1p   = pset->zA(i1,i4) * pset->zB(i3,i4);
   cplx z2p   = pset->zA(i1,i5) * pset->zB(i3,i5);
   cplx ztp   = pset->zA(i1,i2) / pset->zA(i4,i5) / pset->zA(i1,i5);
   cplx zampp = ztp *(z1p + z2p);
   // g-
   cplx z1m   = pset->zB(i1,i4) * pset->zA(i1,i2);
   cplx z2m   = pset->zB(i4,i5) * pset->zA(i2,i5);
   cplx ztm   = pset->zB(i4,i3) / pset->zB(i4,i5) / pset->zB(i1,i5);
   cplx zampm = ztm *(z1m + z2m);

   double amp = pow(abs(zampp),2) + pow(abs(zampm),2);
   
   return 2.0*amp*prop;
}

// 2 -> 6


// Propagators

double propagatorVBF(const double s1, const double s2, const int iboson) {
   // Only W boson to Higgs right now so iboson is irrelevant
   return propagatorWFH(s1, s2);
}

double propagatorWFH(const double s1, const double s2) { 
   double emw     = 80.398;
   double ewwidth = 2.1054;
   double stw     = 0.22264585341299603;
   double p1      = propagator(s1, emw, ewwidth);
   double p2      = propagator(s2, emw, ewwidth);
   return emw*emw/p1/p2/pow(stw,3);
}

double propagator(const double s, const double mass, const double width) {
   double t1 = pow(s - mass*mass, 2);
   double t2 = pow(mass*width, 2);
   return t1 + t2;
}
