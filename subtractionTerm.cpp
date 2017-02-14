#include "subtractionTerm.h"
#include "matrixElement.h"

double subtractionTerm(MomentumSet *pset, const double ptcut, const double rkt, const int minjets) {
   switch (pset->npar) {
      case(7) :
         return C1g0WFHS(pset, ptcut, rkt, minjets);
      case(8): 
         return C2g0WFHS(pset, ptcut, rkt, minjets);
   }
   return 0.0;
}

// npar = 7
double C1g0WFHS(MomentumSet *pset, const double ptcut, const double rkt, const int minjets) {
   int i1, i2, i3, i4, i5;
   pset->getID(&i1, &i2, &i3, &i4, &i5);
   // We have 2 possible limits:
   double lim1 = 0.0;
   double lim2 = 0.0;
   // i1 - i5 - i4
   MomentumSet reducedpmap1 = pset->mapIF(i1, i5, i4);
   int ifail = reducedpmap1.apply_cuts(ptcut, rkt, minjets);
   if (!ifail) {
      // set up parton identities manually, safer than guessing
      reducedpmap1.setID(1,3,4,2);
//       double antenae   = A30(pset, i1, i5, i4);
       double reducedME = C0g0WFH(&reducedpmap1); 
//       lim1 = antenae*reducedME;
   }
   // 
   
   // i2 - i5 - i3
   // MomentumSet reducedpmap2 = pset->mapIF(i2, i5, i3);
   // int ifail = reducedpmap2.apply_cuts(ptcut, rkt, minjets);
   // if (!ifail) = {
   //    double antenae   = A30(pset, i2, i5, i3);
   //    double reducedME = C0g0WFH(reducedpmap2); 
   //    lim2 = antenae*reducedME;
   // }
   return lim1 + lim2;
}

// npar = 8
double C2g0WFHS(MomentumSet *pset, const double ptcut, const double rkt, const int minjets) {
   // We have many possible limits:
   double lim1 = 0.0;
   double lim2 = 0.0;
   return lim1 + lim2;
}
