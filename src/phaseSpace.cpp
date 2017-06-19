#include "Runcard.h"
#include "PhaseSpace.h"
#include <iostream>
#include <string.h>

using namespace std;

MomentumSet phaseSpace(int npar, double s_input, const cubareal x[]) {
   // Parameter
   int nrand_p3 = 5;
   int nrand_out = nrand_p3 + 4;
   // Weight
   double wtps = 1.0;
   double mh   = 125.0;
   double mh2  = pow(mh, 2);

   // Generate shat from x[0], x[1] 
   double taumin = mh2 / s_input;
   double taumax = 1.0;
   double tau    = pickRand(2, taumin, taumax, x[0], &wtps);
   double x1     = pow(tau, x[1]);
   double x2     = tau/x1;
   wtps = -wtps*log(tau);
   if (UNIT_PHASE) {
      x1   = 1.0;
      x2   = 1.0;
      wtps = 1.0;
   }
   double shat = x1*x2*s_input;

   // Generate Higgs mass and photons 1 and 2 from x[2], x[3], x[4]
   if ( mh2 >= shat ) {
      return MomentumSet(0);
   }
   double shiggs = generateInwa(4, x[2], shat, mh, &wtps);
   double gamcos = pickRand(1, COSTHMIN, COSTHMAX, x[3], &wtps);
   double gamphi = pickRand(1, PHIMIN, PHIMAX, x[4], &wtps);
   wtps = wtps/2.0/M_PI;
   wtps = wtps/32.0/M_PI/M_PI;

   // generate the 2 -> 3 system from x[5], x[6], x[7], x[8]
   vector <FourMomentum> pset;
   double s1, s2, smin, smax;
   pset.reserve(5);
   switch (npar) {
      case 6:
         s1 = 0.0;
         s2 = 0.0;
         break;
      case 7:
         smin = Y0*shat;
         smax = pow(sqrt(shat) - sqrt(shiggs), 2) - smin;
         if (smax <= smin) return MomentumSet(0);
         s1 = pickRand(2, smin, smax, x[nrand_out], &wtps);
         s2 = 0.0;
         break;
      case 8: // Triplecolinear case, this is the same for 7 and 8
         smin = Y0*shat;
         smax = pow(sqrt(shat) - sqrt(shiggs), 2) - smin;
         if (smax <= smin) return MomentumSet(0);
         s1   = pickRand(2, smin, smax, x[nrand_out], &wtps);
         s2   = 0.0;
         break;
   }
   // input x[5:8], s[1:4], 
   int ipass = p3generic_nj(x, nrand_p3, shat, s1, shiggs, s2, &pset, &wtps);
   if ( ipass != 1 ) {
      // Bad phasespace point
      return MomentumSet(ipass);
   }
   
   if (DEBUG) {
      cout << "Back in phaseSpace" << endl;
      cout << "Pset[0]: " <<  pset[0] << endl;
      cout << "Pset[1]: " <<  pset[1] << endl;
      cout << "Pset[2]: " <<  pset[2] << endl;
      cout << "Pset[3]: " <<  pset[3] << endl;
      cout << "Pset[4]: " <<  pset[4] << endl;
   }

   // Decay the Higgs (pset[3]) into two pgamma
   FourMomentum pg1, pg2;
   makePs2cm_nj(shiggs, gamcos, gamphi, 0.0, 0.0, &pg1, &pg2);

   if (DEBUG) {
      cout << "Decay of the Higgs" << endl;
      cout << pg1 << endl;
      cout << pg2 << endl;
      cout << "cos: " << gamcos << " phi: " << gamphi << endl;
   }
   pg1.transformation(pset[3].unboost_matrix());
   pg2.transformation(pset[3].unboost_matrix());
   if (DEBUG) {
      cout << "after boost: " << endl;
      cout << pg1 << endl;
      cout << pg2 << endl;
   }

   // All done, generate the momentum set
   vector <FourMomentum> p_out;
   p_out.reserve(npar);
   p_out.emplace_back(pset[0]);
   p_out.emplace_back(pset[1]);

   // Todo: refactoring
   double cos14, phi14;
   FourMomentum pout1, pout2;
   switch (npar) {
      case 6:
         p_out.emplace_back(pset[2]);
         p_out.emplace_back(pset[4]);
         break;
     // For n>6 decay of the blobs coming out from p3generic (if necessary)
      case 7:
         cos14 = pickRand(1, COSTHMIN, COSTHMAX, x[nrand_out + 1], &wtps);
         phi14 = pickRand(1, PHIMIN, PHIMAX, x[nrand_out + 2], &wtps);

         wtps = wtps/pow(4.0*M_PI, 3);
         makePs2cm_nj(s1, cos14, phi14, 0.0, 0.0, &pout1, &pout2);
         pout1.transformation(pset[2].unboost_matrix());
         pout2.transformation(pset[2].unboost_matrix());

         // Store 
         p_out.emplace_back(pout1);
         p_out.emplace_back(pout2);
         p_out.emplace_back(pset[4]);
         break;
      case 8:
         FourMomentum pout3, ptemp;
         cos14 = pickRand(1, COSTHMIN, COSTHMAX, x[nrand_out + 1], &wtps);
         phi14 = pickRand(1, PHIMIN, PHIMAX, x[nrand_out + 2], &wtps);
         // And now generate the second blob
         smax = s1 - smin;
         if (smax <= smin) return MomentumSet(0);
         s2   = pickRand(2, smin, smax, x[nrand_out+3], &wtps);

         // pinch out one of the particles (pout1) and generate the second blob (ptemp)
         wtps = wtps/pow(4.0*M_PI, 3);
         makePs2cm_nj(s1, cos14, phi14, s2, 0.0, &ptemp, &pout1);
         ptemp.transformation(pset[2].unboost_matrix());
         pout1.transformation(pset[2].unboost_matrix());

         // separate the second blob into massless partons pout2, pout3
         double cos56 = pickRand(1, COSTHMIN, COSTHMAX, x[nrand_out + 4], &wtps);
         double phi56 = pickRand(1, PHIMIN, PHIMAX, x[nrand_out + 5], &wtps);

         wtps = wtps*sqrt(dlambda(s1, 0.0, s2))/s1/pow(4.0*M_PI, 3);
         makePs2cm_nj(s2, cos56, phi56, 0.0, 0.0, &pout2, &pout3);
         pout2.transformation(ptemp.unboost_matrix());
         pout3.transformation(ptemp.unboost_matrix());

         // Store 
         p_out.emplace_back(pout1);
         p_out.emplace_back(pout2);
         p_out.emplace_back(pout3);
         p_out.emplace_back(pset[4]);
         break;
   }

   p_out.emplace_back(pg1);
   p_out.emplace_back(pg2);

   return MomentumSet(npar, p_out, wtps, x1, x2);

}

// Sampling functions
double generateInwa(const int itype, const double r, const double shat, double mh, double *wtps) {
   double s12 = 0.0;
   if (UNIT_PHASE) {
      double smin = Y0*shat;
      double smax = shat;
      s12 = pickRand(1, smin, smax, r, wtps);
      //        cout << "GenerateInwa" << endl;
      //        cout << smin << endl;
      //        cout << smax << endl;
      //        cout << r << endl;
      //        cout << s12 << endl;
   } else {
      s12 = pow(mh, 2);
      *wtps = (*wtps) * pow(4.0*M_PI, 2);
   }
   return s12;
}

double pickRand(const int itype,  const double smin, const double smax, const double r, double *wtps) {
   double res;
   switch(itype) {
      case 1: // Linear
         res   = smin + r*(smax - smin);
         if(DEBUG) {
            cout << "Pickrand" << endl;
            cout << "  smin: " << smin << endl;
            cout << "  smax: " << smax << endl;
            cout << "  res: " << res << endl;
         }
         *wtps = (*wtps)*(smax - smin);
         break;
      case 2:
         res   = smin*pow((smax/smin), r);
         *wtps = (*wtps)*log(smax/smin)*res;
         break;
   }
   return res;
} 

// Generate 4-momenta
int p3generic_nj(const cubareal x[], const int n_i, const double shat, const double s1, const double s2, const double s3, vector<FourMomentum> *pset, double *wtps) {
   using namespace std;
   // Parameters
   double rma2 = 0.0;
   double rmb2 = 0.0;
   // Parse input
   double roots= sqrt(shat);
   double Eab  = roots/2.0;
   double pin  = roots/2.0;
   double x1   = x[n_i];
   double x2   = x[n_i+1];
   double x3   = x[n_i+2];
   double x4   = x[n_i+3];

   //    cout << "Momentum in p3generic_nj" << endl;
   //    cout << roots << endl;
   //    cout << s1    << endl;
   //    cout << s2    << endl;
   //    cout << s3    << endl;

   // Incoming momenta
   FourMomentum pa = FourMomentum(0.0, 0.0, pin, Eab);
   FourMomentum pb = FourMomentum(0.0, 0.0, -pin, Eab);

   // Sample t1 (cos1a)
   double s23 = pow(sqrt(s2) + sqrt(s3), 2);
   double t1max, t1min;
   glimits(shat, 0.0, s23, rma2, rmb2, s1, &t1max, &t1min);
   t1max = min(t1max, Y0*t1min); 
   double t1 = - pickRand(2,-t1max,-t1min, x1, wtps);
   if (DEBUG) { 
      cout << "t1 sampling:" << endl;
      cout << t1max << endl;
      cout << t1min << endl;
      cout << t1 << endl;
   }

   // Sample t2 (cos3b)
   s23 = (s1*t1 - shat*t1 - t1*t1) / (s1 - t1);
   double t2min, t2max;
   glimits(s23, 0.0, s2, rmb2, t1, s3, &t2max, &t2min);
   t2max = min(t2max, Y0*t2min);
   double t2 = - pickRand(2,-t2max,-t2min, x2, wtps);
   if (DEBUG) {
      cout << "t2 sampling:" << endl;
      cout << "   " << t2max << endl;
      cout << "   " << t2min << endl;
      cout << "   " << t2 << endl;
   }

   // Sample s23
   double rnum1  =  s3*(s2 - t1) - (s2 + s3 + t1)*t2 + t2*t2;
   double rnum2  =  (-s3 + t2)*sqrt(dlambda(s2,t1,t2));
   double s23min =  - (rnum1 - rnum2)/2.0/t2;
   double s23max =  s23;
   s23min        =  max(s23min, pow(sqrt(s2) + sqrt(s3),2));
   s23           =  pickRand(1, s23min, s23max, x3, wtps);
   if (DEBUG) {
      cout << "s23 sampling:" << endl;
      cout << "   " << s23max << endl;
      cout << "   " << s23min << endl;
      cout << "   " << s23 << endl;
   }

   // Sample s12
   double detcom = 2.0*(s2*shat*(s23 - t1) +
         t1*(2.0*s3*shat - s23*(s3 + shat) + (s3 + shat)*t1) + 
         s1*(s23 - t1)*(s3 - t2) - 
         (-s23*s23 + shat*t1 + s23*(shat + t1))*t2) ;
   double detsqr = 16.0*shat*(s1*(s23 - t1) + t1*(-s23 + shat + t1))*
      (s2*(s23 - t1)*(s3 - t2) + 
       (s23 - s3 - t1 + t2)*(-(s3*t1) + s23*t2)) ;
   double rden1 = 2.0*pow(s23 - t1, 2);
   //    double rden1  = 2.0*(s23**2 - 2d0*s23*t1 + t1**2)
   double s12max = (detcom + sqrt(detsqr)) / rden1;
   double s12min = (detcom - sqrt(detsqr)) / rden1;
   double s12    = pickRand(1, s12min, s12max, x4, wtps);

   // Check all s_max > s_min before filling the rest of the thing
   if ( (s12min >= s12max) || (s23min >= s23max) 
         || (t2min >= t2max)   || (t1min >= t1max) ) {
      cout << "max greater than min!" << endl;
      cout << "s12: " << s12min << " " << s12max << endl;
      cout << "s23: " << s23min << " " << s23max << endl;
      cout << "t2 : " << t2min << " " << t2max << endl;
      cout << "t1 : " << t1min << " " << t1max << endl;
      return 0;
   }

   double dl1    = dlambda(s23,t1,rmb2);
   double rg1    = gres(shat,t1,s23,rma2,rmb2,s1);
   double rg2    = gres(s23,t2,s3,t1,rmb2,s2);
   double r1     = sqrt(rg1*rg2)*2.0;
   double det3   = -(s23*(s2*shat + (s3 + shat)*t1)) -
      s1*(s23 - t1)*(s3 - t2) + 
      s23*s23*(s3 + shat - t2) + s23*(shat + t1)*t2 + 
      shat*t1*(s2 - 2.0*s3 + t2);
   double cp  = (s12 - shat - s3 + 1.0/dl1*det3)*dl1/r1;
   double rph = acos(cp);

   // Generate p1 and ps23
   double E1    = (shat + s1 - s23)/2.0/roots;
   double E23   = (shat + s23 - s1)/2.0/roots;
   double dlin  = sqrt(dlambda(shat, rma2, rmb2));
   double dlout = sqrt(dlambda(shat, s23, s1));
   double p1out = dlout/2.0/roots;
   double cos1a = (t1 - s1 - rma2 + 2.0*Eab*E1)/(2.0*pin*p1out);
   double sin1a = sqrt(1.0 - pow(cos1a,2));
   double rph1  = pickRand(1, PHIMIN, PHIMAX, 0.0, wtps);

   FourMomentum p1 = FourMomentum(
         p1out*sin1a*cos(rph1),
         p1out*sin1a*sin(rph1),
         p1out*cos1a,
         E1);
   FourMomentum p23 = FourMomentum(-p1.px, -p1.py, -p1.pz, E23);

   // Generate p2 and p3
   double dl2   = sqrt(dlambda(s23,s3,s2));
   double dl3   = sqrt(dlambda(s23, rmb2, t1));
   double cos3b = (t2-rmb2-s3+(s23+rmb2-t1)*(s23+s3-s2)/2.0/s23 )
      * 2.0*s23/dl3/dl2;
   double sin3b = sqrt(1.0 - pow(cos3b, 2));

   FourMomentum p2, p3;
   makePs2cm_nj(s23, cos3b, rph, s2, s3, &p2, &p3);

   // Boost p2 and p3 back to c.o.m. frame
   //    matrix unboost_mat = unboostrest(&p23);

   if(DEBUG) {
      cout << "p2 and p3 before unboost" << endl;
      cout << p2 << endl;
      cout << p3 << endl;
      cout << "p23^2: " << p23.sq() << endl;
   }
   p2.transformation(p23.unboost_matrix());
   p3.transformation(p23.unboost_matrix());

   // Compute weight
   double r3 = pow(s12*s23 - s1*s3 - s12*t1 + s3*t1 + s1*t2 - s23*t2,2);
   double r4 = shat*shat*(s2*s2 + pow(t1 - t2,2) - 2.0*s2*(t1 + t2));
   double r5 = shat*(2.0*(s3*(s1 - t1)*(-s2 + t1) + 
            (s2*(s23 - 2.0*t1) + (s23 + s3)*t1 + 
             s1*(s2 - 2.0*s23 + s3 + t1))*t2 - (s1 + s23)*t2*t2 + 
            s12*(s2*(-s23 + t1) + t1*(-2.0*s3 - t1 + t2) +
               s23*(t1 + t2))));
   double r2 = (r3+r4+r5) / 16.0;
   if (r2 >= 0) {
      if (DEBUG) cout << "r2: " << r2 << endl;
      return 0;
   }

   *wtps = (*wtps)/16.0/2.0/sqrt(dlambda(shat, rma2, rmb2))/sqrt(-r2)/pow(2.0*M_PI,5);

   // p23 ---> p2 p3
   (*pset).emplace_back(pa);
   (*pset).emplace_back(pb);
   (*pset).emplace_back(p1);
   (*pset).emplace_back(p2);
   (*pset).emplace_back(p3);

   if (DEBUG) {
      cout << "List of momenta out of p3generic" << endl;
      cout << "pa: " << pa << endl;
      cout << "pb: " <<  pb << endl;
      cout << "p1: " <<  p1 << endl;
      cout << "p2: " <<  p2 << endl;
      cout << "p3: " <<  p3 << endl;
      cout << p23 << endl;
   }

   // Save all pa, pb, p1, p2, p3 to pres and return
   return 1;
}

void makePs2cm_nj(const double s, const double cos12, const double phi, const double s1, const double s2, FourMomentum *p1, FourMomentum *p2) {
   double sin12 =  sqrt(1.0 - pow(cos12,2));
   double sph   =  sin(phi);
   double cph   =  cos(phi);
   double smin  =  s1 + s2;
   double roots =  sqrt(s);

   if (s < smin) std::cout << "Improper call of makePs2cm_nj" << std::endl;

   double E1 = (s + s1 - s2) / 2.0 / roots;
   double E2 = (s + s2 - s1) / 2.0 / roots;
   double pp = sqrt(E1*E1 - s1);
   if (DEBUG) {
      cout << "s1 and s2: " ;
      cout << s1 << " " << s2 << endl;
      cout << E1 << " " << E2 << endl;
      cout << "s: " << s << endl;
   }

   p1->E  = E1 ;
   p1->px = pp*sin12*cph ; 
   p1->py = pp*sin12*sph ; 
   p1->pz = pp*cos12;

   p2->E  = E2 ;
   p2->px = - (p1->px);
   p2->py = - (p1->py);
   p2->pz = - (p1->pz);

   return;
}

// G-Determinant
void glimits(const double x, const double y, const double z, const double u, const double v, const double w, double *ymax, double *ymin) {
   double dl1 = sqrt(dlambda(x, u, v));
   double dl2 = sqrt(dlambda(x, w, z));

   double ycom = u + w - 1.0/2.0*(x+u-v)*(x+w-z)/x;
   double yvar = dl1*dl2 / 2.0 / x;

   *ymax = ycom + yvar;
   *ymin = ycom - yvar;
}

double gres(const double x, const double y, const double z, const double u, const double v, const double w) {
   double yp, ym;

   glimits(x, y, z, u, v, w, &yp, &ym);
   return x*(y-yp)*(y-ym);
}

// General utilities
double dlambda(const double a, const double b, const double c) {
   using namespace std;
   if (abs(c) < min(abs(a),abs(b))) {
      return pow(a-b,2) + c*c - 2.0*(a+b)*c;
   } else if ( abs(c) < min(abs(a),abs(c)) ) {
      return pow(a-c,2) + b*b - 2.0*(a+c)*b;
   } else {
      return pow(b-c,2) + a*a - 2.0*(b+c)*a;
   }
}
