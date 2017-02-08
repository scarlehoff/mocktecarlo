#include <algorithm>

#include "MomentumSet.h"


// Constructor
MomentumSet::MomentumSet( int n_in, momentum_t p_in[] ) {
   npar    = n_in;
   pset    = new momentum_t[npar];
   cplx_array spinorA(boost::extents[npar][npar]);
   compute_spinors();
}

void MomentumSet::compute_spinors() {
   int ss[npar];
   std::fill_n(ss, npar, 1);
   ss[0] = -1;
   ss[1] = -1;
   for (int i = 0; i < npar-1 ; i++) {
      spinorA[i][i] = 0.0;
      for (int j = i + 1; j < npar; j++) {
         spinorA[i][j] = eval_zA(i,j);
         spinorA[j][i] = - spinorA[i][j];
      }
   }
}

cplx MomentumSet::eval_zA(int i, int j) {	
	return 0.0;
}
	
