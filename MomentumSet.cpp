#include "MomentumSet.h"
#include <string.h>
#include <iostream>


using namespace std;
// Constructor
MomentumSet::MomentumSet( int n_in, momentum_t p_in[] ) {
   npar    = n_in;
   // not convinced
   pset    = new momentum_t[npar];
   memcpy(pset, p_in, sizeof(pset));
   // -- 
   cplx_array spinorA(boost::extents[npar][npar]);
   compute_spinors();
}

// Getters
cplx MomentumSet::zA(int i, int j) {
    return spinorA[i][j];
}
cplx MomentumSet::zB(int i, int j) {
    double ss = 1.0;
    if (i < 2) ss = - ss;
    if (j < 2) ss = - ss;
    return - ss*conj(spinorA[i][j]);
}
double MomentumSet::s(int i, int j) {
    cplx res = (zA(i,j)*zB(i,j));
    return res.real();
}

// Debug
void MomentumSet::printAll() {
    cout << "----DEBUG-----" << endl;
    cout << "Momentum Set:" << endl;
    momentum_t a = { 1.0,2.0,3.0,5.0 }; 
    for (int i = 0; i < npar; i++) {
        cout << "p[" << (i+1) << "] = ";
        cout << pset[i].px << " " << pset[i].py << " " << pset[i].pz << " " << pset[i].E << endl;
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

// Destructor
void MomentumSet::cleanmem() { 
    delete(pset);
}

// Private functions: Spinors
void MomentumSet::compute_spinors() {
   for (int i = 0; i < npar-1 ; i++) {
      spinorA[i][i] = 0.0;
      for (int j = i + 1; j < npar; j++) {
         spinorA[i][j] = eval_zA(i,j);
         spinorA[j][i] = - spinorA[i][j];
      }
   }
}
cplx MomentumSet::eval_zA(int i, int j) {
    cplx zi = cplx(0.0, 1.0);
    momentum_t *a, *b;
    a = &(pset[i]);
    b = &(pset[j]);

    double at2, at, ap, am;
    cplx zea;
    at2 = pow((*a).px,2) + pow((*a).py,2);
    at  = sqrt(at2);
    ap  = (*a).E + (*a).pz;
    am  = (*a).E - (*a).pz;

    if ( ap < ((*a).E/2.0) ) ap = at2 / am ;
    if ( am < ((*a).E/2.0) ) am = at2 / ap ;

    if ( at == 0.0 ) {
        zea = 0.0;
    } else {
        zea = ( (*a).px + zi*(*a).py ) /at;
    }

    double bt2, bt, bp, bm;
    cplx zeb;
    bt2 = pow((*b).px,2) + pow((*b).py,2);
    bt  = sqrt(bt2);
    bp  = (*b).E + (*b).pz;
    bm  = (*b).E - (*b).pz;

    if ( bp < ((*b).E/2.0) ) bp = bt2 / bm ;
    if ( bm < ((*b).E/2.0) ) bm = bt2 / bp ;

    if ( bt == 0.0 ) {
        zeb = 0.0 ;
    } else {
        zeb = ( (*b).px + zi*(*b).py ) /bt ;
    }

    cplx res = sqrt(am*bp)*zea - sqrt(ap*bm)*zeb;
    if ( i == 1  || j == 1) res = res*zi ;
    if ( i == 2  || j == 2) res = -res*zi ;

	return res;
}

// Private functions: cuts
