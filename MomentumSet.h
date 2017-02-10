#pragma once
#include <complex>
#include <algorithm>
#include "boost/multi_array.hpp"
#include "FourVector.h"

typedef std::complex <double> cplx;
typedef boost::multi_array<cplx, 2> cplx_array;
  
//cplx_array spinorA(boost::extents[3][4]);

class MomentumSet {
   public:
      // Constructor
      MomentumSet(const int error);
      MomentumSet(const int n_in, const std::vector <FourMomentum> p_in, const double wt, const double x_1, const double x_2);

      int npar, njets;
      int ifail;
      double x1, x2, weight;
      std::vector <FourMomentum> pset;

      // Getter
      const cplx zA(const int i, const int j);
      const cplx zB(const int i, const int j);
      const double s(const int i, const int j);

      // Debug
      void printAll();

      // Destructor
      void cleanmem();

   private:
      cplx_array spinorA;
      void compute_spinors();
      void apply_cuts();
	   const cplx eval_zA(int i, int j);
};


   
