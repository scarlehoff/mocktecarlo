#include "FourVector.h"
#include <iostream>

using namespace std;

// Four Vector Super class

// Constructor
FourVector::FourVector() {
   xname = " x";
   yname = " y";
   zname = " z";
   tname = " t";
} 
FourVector::FourVector(double xin, double yin, double zin, double tin) {
   x = xin;
   y = yin;
   z = zin;
   t = tin;
}

// Utilities
double FourVector::sq() {
   return (*this)*(*this);
}

// Operator overloading
FourVector FourVector::operator + (FourVector &v) const {
   return FourVector(x + v.x, y + v.y, z + v.z, t + v.t);
}
FourVector FourVector::operator - (FourVector &v) const {
   return FourVector(x - v.x, y - v.y, z - v.z, t - v.t);
}

double FourVector::operator * (FourVector &v) const {
   return v.t*t - (v.x*x + v.y*y + v.z*z);
}

FourVector FourVector::operator * (const double k) const {
   return FourVector(k*x,k*y,k*z,k*t);
}
FourVector FourVector::operator / (const double k) const {
   return FourVector(x/k,y/k,z/k,t/k);
}

// Global
FourVector operator * (const double k, FourVector &v) { return v*k; }
FourMomentum operator * (const double k, FourMomentum &v) { return v*k; }

ostream& operator << (ostream& os, const FourVector &v) {
   cout << v.xname << " = " << v.x;
   cout << v.yname << " = " << v.y;
   cout << v.zname << " = " << v.z;
   cout << v.tname << " = " << v.t;
   cout << endl;
}

//////////// Four Momentum
FourMomentum::FourMomentum() {
   xname = " px";
   yname = " py";
   zname = " pz";
   tname = "  E";
   p[0]  = &E;
   p[1]  = &px;
   p[2]  = &py;
   p[3]  = &pz;
}
FourMomentum::FourMomentum(double xin, double yin, double zin, double tin) :
   FourVector(xin, yin, zin, tin) {
   xname = " px";
   yname = " py";
   zname = " pz";
   tname = "  E";
}
FourMomentum::FourMomentum(const FourMomentum &p) :
   FourVector(p.x, p.y, p.z, p.t) {
   xname = " px";
   yname = " py";
   zname = " pz";
   tname = "  E";
}
matrix FourMomentum::unboost_matrix() {
   double v[4];
   matrix res(boost::extents[4][4]);

   double s = sq();
   if (s < 0.0) {
       std::cout << "Undue boost" << std::endl;
       std::cin.ignore();
   }
   double roots = sqrt(s);
   double gamma = E / roots;
   v[0] = -1.0;
   v[1] = -px / E;
   v[2] = -py / E;
   v[3] = -pz / E;
   double v2 = 0.0;
   for (int i = 1; i < 4; i++) v2 += v[i]*v[i];

   res[0][0] = gamma;
   for (int i = 1; i < 4; i++) {
       res[i][0] = -gamma*v[i];
       res[0][i] = -gamma*v[i];
       for (int j = 1; j < 4; j++)  res[j][i] = (gamma - 1.0)*v[i]*v[j]/v2;
       res[i][i] = res[i][i] + 1.0;
   }
   return res;
}

void FourMomentum::transformation(matrix mat) {
   double tmp, pt[4];
   for (int i = 0; i < 4; i++) {
      tmp = 0.0;
      for (int j = 0; j < 4; j++) {
         tmp += mat[i][j]*(*p[j]);
      }
      pt[i] = tmp;
   }

   for (int i = 0; i < 4; i++) {*p[i] = pt[i]; }
   return;
}

double FourMomentum::computeKin() {
   // pt
   double ptx = px;
   double pty = py;
   pt2 = px*px + py*py;
   pt  = sqrt(pt2);
   // rapidity
   double yp = E + pz;
   double ym = E - pz;
   yrap = log(yp/ym)/2.0;
   // azimth
   phi  = atan(py/px);

}

// Overrides
FourMomentum FourMomentum::operator * (const double k) const { return FourMomentum(k*x,k*y,k*z,k*t); }
FourMomentum FourMomentum::operator / (const double k) const { return FourMomentum(x/k,y/k,z/k,t/k); }
FourMomentum FourMomentum::operator + (FourMomentum &v) const { return FourMomentum(x + v.x, y + v.y, z + v.z, t + v.t); }
FourMomentum FourMomentum::operator - (FourMomentum &v) const { return FourMomentum(x - v.x, y - v.y, z - v.z, t - v.t);}


