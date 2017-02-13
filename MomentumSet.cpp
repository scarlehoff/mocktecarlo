#include "MomentumSet.h"
#include <string.h>
#include <iostream>
#include "fastjet/ClusterSequence.hh"

using namespace std;
// Constructor
MomentumSet::MomentumSet(const int error) {
   ifail = 1;
}

MomentumSet::MomentumSet(const int n_in, const vector<FourMomentum> p_in, const double wt, const double x_1, const double x_2) {
   npar    = n_in;
   weight  = wt;
   x1      = x_1;
   x2      = x_2;

   pset.reserve(p_in.size());
   copy(p_in.begin(), p_in.end(), back_inserter(pset));

//   compute_spinors();
   ifail = 0;
}

// Getters
const cplx MomentumSet::zA(const int i, const int j) {
   int ii = i - 1;
   int ij = j - 1;
//   cplx_array spinorA(boost::extents[npar][npar]);
//   return spinorA[i][j];
   return eval_zA(ii,ij);
}
const cplx MomentumSet::zB(const int i, const int j) {
   double ss = 1.0;
   if (i < 3) ss = - ss;
   if (j < 3) ss = - ss;
// This is broken at the moment I don't know why: Todo
//   cplx_array spinorA(boost::extents[npar][npar]);
//   return - ss*conj(spinorA[i][j]);
   int ii = i - 1;
   int ij = j - 1;
   return - ss*conj(eval_zA(ii,ij));
}
const double MomentumSet::s(const int i, const int j) {
   cplx res = (zA(i,j)*zB(j,i));
   return res.real();
}

// Debug
void MomentumSet::printAll() {
   cout << "----DEBUG-----" << endl;
   cout << "Momentum Set:" << endl;
   for (int i = 0; i < npar; i++) {
      cout << "p[" << (i+1) << "] : " << pset[i] << endl;
   }

   cout << "Invariants s_ij:" << endl;
   for (int i = 0; i < npar - 1; i++) {
      for (int j = i+j; j < npar; j++) {
         cout << "s_" << i << j << " = " << s(i,j) << endl;
      }
   }
   cout << "Spinors:" << endl;
   for (int i = 0; i < npar - 1; i++) {
      for (int j = i+j; j < npar; j++) {
         cout << "zA_" << i << j << " = " << zA(i,j) << endl;
         cout << "zB_" << i << j << " = " << zB(i,j) << endl;
      }
   }
}

// Private functions: Spinors
void MomentumSet::compute_spinors() {
   cplx_array spinorA(boost::extents[npar][npar]);
   for (int i = 0; i < (npar-1) ; i++) {
      spinorA[0][0] = 0.0;
      for (int j = (i + 1); j < npar; j++) {
         spinorA[i][j] = eval_zA(i,j);
//         cout << i << j << " " << spinorA[i][j] << endl;
         spinorA[j][i] = - spinorA[i][j];
      }
   }
}
const cplx MomentumSet::eval_zA(int i, int j) {
   cplx zi = cplx(0.0, 1.0);
   FourMomentum *a, *b;
   a = &(pset[i]);
   b = &(pset[j]);

//   cout << i << ": " << (*a) << endl;
//   cout << j << ": " << (*b) << endl;

   double at2, at, ap, am;
   cplx zea;
   at2 = pow(a->px, 2) + pow(a->py, 2);
   at  = sqrt(at2);
   ap  = a->E + a->pz;
   am  = a->E - a->pz;

   if ( ap < (a->E/2.0) ) ap = at2 / am ;
   if ( am < (a->E/2.0) ) am = at2 / ap ;

   if ( at == 0.0 ) {
      zea = {1.0, 0.0} ;
   } else {
      zea = {a->px, a->py};
      zea = zea/at;
   }

   double bt2, bt, bp, bm;
   cplx zeb;
   bt2 = pow(b->px, 2) + pow(b->py, 2);
   bt  = sqrt(bt2);
   bp  = b->E + b->pz;
   bm  = b->E - b->pz;

   if ( bp < (b->E/2.0) ) bp = bt2 / bm ;
   if ( bm < (b->E/2.0) ) bm = bt2 / bp ;

   if ( bt == 0.0 ) {
      zeb = {1.0, 0.0} ;
   } else {
      zeb = {b->px, b->py};
      zeb = zeb/bt;
   }

   cplx res = sqrt(am*bp)*zea - sqrt(ap*bm)*zeb;
   if ( i == 1  || j == 1) res = res*zi ;
   if ( i == 2  || j == 2) res = -res*zi ;

   return res;
}

// Public functions: cuts
int MomentumSet::apply_cuts(const double ptcut, const  double rkt, const int minjet) {
   using namespace fastjet;
   // Copy the jets we care about to fastjet
   vector <PseudoJet> particles;
   for (int i = 2; i < (npar - 2); i++) {
      particles.push_back(PseudoJet(pset[i].px, pset[i].py, pset[i].pz, pset[i].E));
   }
   // Define jets
   JetDefinition jet_def(antikt_algorithm, rkt);

   // Run the clustering and get the jets back sorted by pt
   ClusterSequence cs(particles, jet_def);
   vector <PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

   // Do we have enough jets?
   njets = jets.size();
   if (njets < minjet) return 1;

   // And do they have enough pt?
   for (int i = 0; i < jets.size(); i++) {
      if (jets[i].pt() < ptcut) njets -= 1;
   }
   if (njets < minjet) return 1;

   return 0;
   // anti-kt
   // Create the diB array:
//   vector <double> diB[2];
//   njets = 0;
//   for (int i = 2; i < (npar-2); i++) {
//      pset[i].computeKin();
//      diB = 1.0/pset[i].pt2;
//      if (pset[i].pt >= ptcut) njets += 1;
//   } // Check whether there is enough pt to go through the cuts while you are at it
//   if (minjet > njets) return 1;
//
//   // Now let's go through the actual antikt
//   double dij;
//   double dmin = *min_element(diB.begin(), diB.end());
//   for (int i = 2; i < (npar-3); i++) {
//      for (int j = 3; i < (npar-2); i++) {
//         deltaij = pow(pset[i].yrap-pset[j].yrap,2) + pow(pset[i].phi-pset[j].phi,2);
//         dij     = min(1.0/pset[i].pt2, 1.0/pset[j].pt2) * deltaij/pow(rmin,2);
//         // Is the smallest dij smaller than diB?
//         if (dij < dmin) {
//            njets -= 1;
//         }
//      }
//   }
//   if (minjet > njets) return 1;
//   return 0;
}
