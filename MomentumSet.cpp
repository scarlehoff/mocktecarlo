#include "MomentumSet.h"
#include <string.h>
#include <iostream>
#include "fastjet/ClusterSequence.hh"

using namespace std;
// Constructor
MomentumSet::MomentumSet(const int error) {
   ifail = 1;
}

MomentumSet::MomentumSet(const int n_in, const vector<FourMomentum> p_in) :
   MomentumSet(n_in, p_in, 0.0, 0.0, 0.0 ){
      ifail = 0;
}

MomentumSet::MomentumSet(const int n_in, const vector<FourMomentum> p_in, const double wt, const double x_1, const double x_2)
  : spinorA(boost::extents[n_in][n_in]) {
   npar    = n_in;
   weight  = wt;
   x1      = x_1;
   x2      = x_2;

   pset.reserve(p_in.size());
   copy(p_in.begin(), p_in.end(), back_inserter(pset));

   compute_spinors();
   ifail = 0;
}

// Setter
//void MomentumSet::setID(const int ii1, const int ii2, const int ii3, const int ii4, const int ii5 = 0, const int ii6 = 0) {
void MomentumSet::setID(int ii1, int ii2, int ii3, int ii4, int ii5, int ii6) {
   i1 = ii1; i2 = ii2;
   i3 = ii3; i4 = ii4;
   i5 = ii5; i6 = ii6;
}

// Getter
const void MomentumSet::getID(int *ii1, int *ii2, int *ii3, int *ii4, int *ii5, int *ii6) {
  *ii1 = i1; *ii2 = i2;
  *ii3 = i3; *ii4 = i4;
  if (ii5 != 0) *ii5 = i5; 
  if (ii6 != 0) *ii6 = i6; 
}

const cplx MomentumSet::zA(const int i, const int j) {
   int ii = i - 1;
   int ij = j - 1;
   return spinorA[ii][ij];
}
const cplx MomentumSet::zB(const int i, const int j) {
   double ss = 1.0;
   if (i < 3) ss = - ss;
   if (j < 3) ss = - ss;
   return - ss * conj(zA(i,j));
}
const double MomentumSet::s(const int i, const int j) {
   cplx res = (zA(i,j)*zB(j,i));
   return res.real();
}

// Mapping
MomentumSet MomentumSet::mapIF(const int ii1, const int ii3, const int ii4) {
   // In order to use fortran notation in matrix elements we need to do this . . . 
   int vi1 = ii1 - 1;
   int vi3 = ii3 - 1;
   int vi4 = ii4 - 1;
   // Map particle i3 unto i1 and i4, i1 always initial hard radiator
   int new_n = npar - 1;
   vector <FourMomentum> p_out;
   p_out.reserve(new_n);

   double omx2 = -s(ii3, ii4)/( s(ii1,ii3) + s(ii1,ii4) );
   double   x2 = 1.0 - omx2;

   // set up initial particles
   if (ii1 == 1) {
      p_out.emplace_back(x2*pset[0]);
      p_out.emplace_back(pset[1]);
   } else if (ii1 == 2) {
      p_out.emplace_back(pset[0]);
      p_out.emplace_back(x2*pset[1]);
   } else {
      cout << "Particle ii1 not initial in mapIF" << endl;
   }

   // Map the rest, ii3 is pinched out and ii4 gets a contribution from p1
   for (int i = 2; i < npar; i++) {
      if (i == vi3) continue;
      if (i == vi4) {
         FourMomentum pini   = omx2*pset[vi1];
         FourMomentum mapped = pset[i] + pset[vi3] - pini;
         p_out.emplace_back(mapped);
      } else {
         p_out.emplace_back(pset[i]);
      }
   }
   return MomentumSet(new_n, p_out);
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
   for (int i = 0; i < (npar-1) ; i++) {
      spinorA[0][0] = 0.0;
      for (int j = (i + 1); j < npar; j++) {
         spinorA[i][j] = eval_zA(i,j);
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
