#include <complex>
#include "boost/multi_array.hpp"
#include <algorithm>

//typedef double momentum_t[4];
struct momentum_t  {
    // todo: Make it a child of "4-vector"
    double E;
    double px;
    double py;
    double pz;

};

typedef std::complex <double> cplx;
typedef boost::multi_array<cplx, 2> cplx_array;
  
//cplx_array spinorA(boost::extents[3][4]);

class MomentumSet {
   public:
      // Constructor
      MomentumSet(int n_in, momentum_t p_in[]);

      int npar, njets;
      double x1, x2, shat;
      momentum_t *pset;

      // Getter
      cplx zA(int i, int j);
      cplx zB(int i, int j);
      double s(int i, int j);

      // Debug
      void printAll();

      // Destructor
      void cleanmem();

   private:
      cplx_array spinorA;
      void compute_spinors();
      void apply_cuts();
	  cplx eval_zA(int i, int j);
};


   