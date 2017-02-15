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
      MomentumSet(const int n_in, const std::vector <FourMomentum> p_in, int i1, int i2, int i3, int i4, int i5 = 0, int i6 = 0);
      MomentumSet(const int n_in, const std::vector <FourMomentum> p_in, const double wt, const double x_1, const double x_2);

      int npar, njets;
      int ifail;
      double x1, x2, weight;
      std::vector <FourMomentum> pset;

      int apply_cuts(const double ptcut, const double rkt, const int minjets);

      // Setter
      void setID(int i1, int i2, int i3, int i4, int i5 = 0, int i6 = 0);

      // Getter
      const cplx zA(const int i, const int j);
      const cplx zB(const int i, const int j);
      const double s(const int i, const int j);
      const void getID(int *i1, int *i2, int *i3, int *i4, int *i5 = 0, int *i6 = 0);

      // Mapping
      MomentumSet mapIF(const int i1, const int i3, const int i4);

      // Debug
      void printAll();
      void fortranOutput();

      // Destructor
      void cleanmem();

   private:
      cplx_array spinorA;
      void compute_spinors();
      void boostToLab();
	   const cplx eval_zA(int i, int j);
      int i1, i2, i3, i4, i5, i6;
};


   
