#include "SubtractionTerm.h"
#include "MatrixElement.h"

using namespace std;

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

   // i1 - i5 - i4
   double lim1 = nadjLimitIF(pset, i1, i5, i4, C0g0WFH, ptcut, rkt, minjets);
   
   // i2 - i5 - i3
   double lim2 = nadjLimitIF(pset, i2, i5, i3, C0g0WFH, ptcut, rkt, minjets);

   return lim1 + lim2;
}

// npar = 8
double C2g0WFHS(MomentumSet *pset, const double ptcut, const double rkt, const int minjets) {
   int i1, i2, i3, i4, i5, i6;
   pset->getID(&i1, &i2, &i3, &i4, &i5, &i6);
   int nadj = 0, adj = 1;
   double limadj = 0.0, limnadj = 0.0;
   if (nadj) {
      // i1 - i5 - i4
      double lim1 = nadjLimitIF(pset, i1, i5, i4, C1g0WFHdown, ptcut, rkt, minjets);
      // i2 - i5 - i3
      double lim2 = nadjLimitIF(pset, i2, i5, i3, C1g0WFHup, ptcut, rkt, minjets);
      // i1 - i6 - i4
      double lim3 = nadjLimitIF(pset, i1, i6, i4, C1g0WFHdown, ptcut, rkt, minjets);
      // i2 - i6 - i3
      double lim4 = nadjLimitIF(pset, i2, i6, i3, C1g0WFHup, ptcut, rkt, minjets);
      limnadj = lim1 + lim2 + lim3 + lim4;
   }

   if (adj) {
      // i1 - i5 - i6
      double lim1 = adjLimitIF(pset, i1, i5, i6, C1g0WFHup, ptcut, rkt, minjets);
      // i1 - i6 - i5
      double lim2 = adjLimitIF(pset, i1, i6, i5, C1g0WFHup, ptcut, rkt, minjets);
      // i2 - i5 - i6
      double lim3 = adjLimitIF(pset, i2, i5, i6, C1g0WFHdown, ptcut, rkt, minjets);
      // i2 - i6 - i5
      double lim4 = adjLimitIF(pset, i2, i6, i5, C1g0WFHdown, ptcut, rkt, minjets);
      // i4 - i5 - i6
      double lim5 = adjLimitFF(pset, i4, i5, i6, C1g0WFHup, ptcut, rkt, minjets);
      // i4 - i6 - i5
      double lim6 = adjLimitFF(pset, i4, i6, i5, C1g0WFHup, ptcut, rkt, minjets);
      // i3 - i5 - i6
      double lim7 = adjLimitFF(pset, i3, i5, i6, C1g0WFHdown, ptcut, rkt, minjets);
      // i3 - i6 - i5
      double lim8 = adjLimitFF(pset, i3, i6, i5, C1g0WFHdown, ptcut, rkt, minjets);
      limadj = lim1 + lim2 + lim3 + lim4 + lim5 + lim6 + lim7 + lim8;
//      limadj = lim1 + lim6;
   }

   return limadj + limnadj;


}

// adj limit IF (i1 initial, i2, i3 soft)
double adjLimitIF(MomentumSet *pset, int i1, int i2, int i3, double (*MatrixE)(MomentumSet*), const double ptcut, const double rkt, const int minjets) {
   MomentumSet reducedpmap1 = pset->mapIF(i1, i2, i3);
   int ifail = reducedpmap1.apply_cuts(ptcut, rkt, minjets);
   if (!ifail) {
      double antenae   = d30(pset, i1, i2, i3);
      double reducedme = MatrixE(&reducedpmap1);
      return antenae*reducedme;
   }
   return 0.0;
}
// adj limit FF (i1 hard, i2, i3 soft)
double adjLimitFF(MomentumSet *pset, int i1, int i2, int i3, double (*MatrixE)(MomentumSet*), const double ptcut, const double rkt, const int minjets) {
   MomentumSet reducedpmap1 = pset->mapFF(i1, i2, i3);
   int ifail = reducedpmap1.apply_cuts(ptcut, rkt, minjets);
   if (!ifail) {
      double antenae   = d30(pset, i1, i2, i3);
      double reducedme = MatrixE(&reducedpmap1);
      return antenae*reducedme;
   }
   return 0.0;
}

// nadj limit (i1 initial, i2 soft)
double nadjLimitIF(MomentumSet *pset, int i1, int i2, int i3, double (*MatrixE)(MomentumSet*), const double ptcut, const double rkt, const int minjets) {
   MomentumSet reducedpmap1 = pset->mapIF(i1, i2, i3);
   int ifail = reducedpmap1.apply_cuts(ptcut, rkt, minjets);
   if (!ifail) {
      double antenae   = A30(pset, i1, i2, i3);
      double reducedme = MatrixE(&reducedpmap1);
      return antenae*reducedme;
   }
   return 0.0;
}

// Antennae

double A30(MomentumSet *pset, const int i1, const int i2, const int i3) {
   double s12  = pset->s(i1, i2);
   double s13  = pset->s(i1, i3);
   double s23  = pset->s(i2, i3);
   double s123 = s12 + s13 + s23;

   double FullAnt = s12/s23 + s23/s12 + 2.0*s13*s123/s12/s23;
   return FullAnt/s123;
}

double d30(MomentumSet *pset, const int i1, const int i2, const int i3) {
   double s12  = pset->s(i1, i2);
   double s13  = pset->s(i1, i3);
   double s23  = pset->s(i2, i3);
   double s123 = s12 + s13 + s23;


   double FullAnt = 2.0*s123*s123*s13/s12/s23 + s12*s13/s23+s13*s23/s12
                    + s23*s23/s12 + 5.0*s123/2.0+s23/2.0;
   return FullAnt/pow(s123,2);
}

