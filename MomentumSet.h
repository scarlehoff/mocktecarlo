#include <complex>
#include "boost/multi_array.hpp"

typedef double momentum_t[4];
typedef std::complex <double> cplx;
typedef boost::multi_array<cplx, 2> cplx_array;
  
cplx_array spinorA(boost::extents[3][4]);

class MomentumSet {
   public:
      // Constructor
      MomentumSet(int n_in, momentum_t p_in[]);

      int npar, njets;
      double x1, x2, shat;
      momentum_t *pset;

      cplx zA(int i, int j);
      cplx zB(int i, int j);
      double   s(int i, int j);

      // Destructor
      void cleanmem();

   private:
      cplx_array spinorA;
      void compute_spinors();
      void apply_cuts();
		cplx eval_zA(int i, int j);
};


   
