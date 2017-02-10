#pragma once
#include <complex>
#include <algorithm>
#include "boost/multi_array.hpp"
#include "FourVector.h"

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
      MomentumSet(const int error);
      MomentumSet(const int n_in, const std::vector <FourMomentum> p_in, const double wt, const double x_1, const double x_2);

      int npar, njets;
      int ifail;
      double x1, x2, weight;
      std::vector <FourMomentum> pset;

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


   
