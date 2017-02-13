#pragma once
#include <string>
#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 2> matrix;

class FourVector {
   public:
      double x;
      double y;
      double z;
      double t;

      // Constructor
      FourVector();
      FourVector(double xin, double yin, double zin, double ti);
      
      // Utilities
      double sq();

      // Operator overloads:
      double operator * (FourVector &v) const;
      FourVector operator * (const double k) const;
      FourVector operator / (const double k) const;

      friend std::ostream& operator << (std::ostream& os, const FourVector &v);

      // Names of the variables for debugging
      std::string xname;
      std::string yname;
      std::string zname;
      std::string tname;
};  
FourVector operator * (const double k, FourVector &v);

class FourMomentum : public FourVector {
   public:
      // Reference to members of FourVector (C++ 11)
      double &px = FourVector::x;
      double &py = FourVector::y;
      double &pz = FourVector::z;
      double &E  = FourVector::t;
      double *p[4]; // Array of pointers towards [E, px, py, pz]
      double pt, pt2, yrap, phi;

      // Constructor
      FourMomentum();
      FourMomentum(double xin, double yin, double zin, double ti);
      FourMomentum(const FourMomentum &p);

      // Some more overloads
      FourMomentum operator * (const double k) const;
      FourMomentum operator / (const double k) const;

      // Utilities
      matrix boost_matrix();
      matrix unboost_matrix();

      // Kinematical variables
      double computeKin();

      // Todo, overload matrix * fourvector instead
      void transformation(matrix mat);
};
FourMomentum operator * (const double k, FourMomentum &v);
