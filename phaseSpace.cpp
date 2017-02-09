#include "phaseSpace.h"
#include <iostream>
#include <string.h>

momentum_t *p3generic_nj() {
}

MomentumSet phaseSpace(int npar, double s_input, const cubareal x[]) {
    // Weight
    double wtps;
    
    // Generate shat from x[1], x[2] "naively"
    double x1 = x[1];
    double x2 = x[2];
    if (UNIT_PHASE) {
        x1 = 1.0;
        x2 = 1.0;
    }
    double shat = x1*x2*s_input;
    
    // Generate Higgs mass and photons 1 and 2 from x[3], x[4], x[5]
    double shiggs = generateInwa(4, x[3], shat, &wtps);
    double gamphi = pickRand(1, x[4], COSTHMIN, COSTHMAX, &wtps);
    double gamcos = pickRand(1, x[5], PHIMIN, PHIMAX, &wtps);
    wtps = wtps/2.0/M_PI;
    wtps = wtps/32.0/M_PI/M_PI;
    
    // generate the 2 -> 5 system from x[6], x[7], x[8], x[9]
    momentum_t pset[5];
    // input x[6:9], s[1:4], 
//    pset = p3generic_nj();
//    check whether the return momentum set is null before continuing

}

// Sampling functions
double generateInwa(const int itype, const double r, const double shat, double *wtps) {
    if (UNIT_PHASE) {
        double smin = y0*shat;
        double smax = shat;
        return pickRand(1, r, smax, smin, wtps);
    } else {
        double s12 = pow(125.0, 2);
        *wtps = (*wtps) * pow(4.0*M_PI, 2);
        return s12;
    }
}

double pickRand(const int itype, const double r, const double smax, const double smin, double *wtps) {
    double res;
    switch(itype) {
        case 1: // Linear
            res   = smin + r*(smax - smin);
            *wtps = (*wtps)*(smax - smin);
            break;
        case 2:
            res   = smin*pow((smax/smin), r);
            *wtps = (*wtps)*log(smax/smin)*res;
            break;
    }
} 

// Generate 4-momenta
momentum_t *p3generic_nj(const cubareal x[], const double s_input[], double *wtps) {
    using namespace std;
    momentum_t pa, pb, p1, p2, p3, p23;
    // Parameters
    double rma2 = 0.0;
    double rmb2 = 0.0;
    // Parse input
    double shat = s_input[0];
    double s1   = s_input[1];
    double s2   = s_input[2];
    double s3   = s_input[3];
    double roots= sqrt(shat);
    double Eab  = roots/2.0;
    double pin  = roots/2.0;

    // Incoming momenta
    pa.E = Eab;
    pb.E = Eab;
    pa.px = 0.0 ; pb.px = 0.0;
    pa.py = 0.0 ; pb.py = 0.0;
    pa.pz = pin ; pb.pz = pin;

    // Sample t1 (cos1a)
    double s23 = pow(sqrt(s2) + sqrt(s3), 2);
    double t1max, t1min;
    // glimits for t1max t1min
    t1max = min(t1max, y0*t1min); 
    double t1 = - pickRand(2, t1max, t1min, x[0], wtps);

    // Sample t2 (cos3b)
    s23 = (s1*t1 - shat*t1 - t1*t1) / (s1 - t1);
    double t2min, t2max;
    // glimits for t2max t2min
    t2max = min(t2max, y0*t2min);
    double t2 = - pickRand(2, t2max, t2min, x[1], wtps);

    // Sample s23
    double rnum1  =  s3*(s2 - t1) - (s2 + s3 + t1)*t2 + t2*t2;
    double rnum2  =  (-s3 + t2)*sqrt(dlambda(s2,t1,t2));
    double s23min =  - (rnum1 - rnum2)/2.0/t2;
    double s23max =  s23;
    s23min        =  max(s23min, pow(sqrt(s2) + sqrt(s3),2));
    s23           =  pickRand(2, s23min, s23max, x[2], wtps);

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
    double s12    = pickRand(2, s12min, s12max, x[3], wtps);

    // Check all s_max > s_min before filling the rest of the thing
    if ( (s12min >= s12max) || (s23min >= s23max) 
      || (t2min >= t2max)   || (t1min >= t1max) ) return 0;

    double dl1    = dlambda(s23,t1,rmb2);
    double rg1    = 1.0; // gres(shat,t1,s23,rma2,rmb2,s1)
    double rg2    = 1.0;// gres(s23,t2,s3,t1,rmb2,s2)
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

    p1.E  = E1;
    p1.px = p1out*sin1a*cos(rph1);
    p1.py = p1out*sin1a*sin(rph1);
    p1.pz = p1out*cos1a;

    p23.E = E23;
    p23.px = -p1.px; 
    p23.py = -p1.py;
    p23.pz = -p1.pz;

    // Generate p2 and p3
    double dl2   = sqrt(dlambda(s23,s3,s2));
    double dl3   = sqrt(dlambda(s23, rmb2, t1));
    double cos3b = (t2-rmb2-s3+(s23+rmb2-t1)*(s23+s3-s2)/2.0/s23 )
                 * 2.0*s23/dl3/dl2;
    double sin3b = sqrt(1.0 - pow(cos3b, 2));

    makePs2cm_nj(s23, cos3b, rph, s2, s3, &p2, &p3);

    // Boost p2 and p3 back to c.o.m. frame
    matrix unboost_mat = unboostrest(&p23);

    // Compute weight
    double r3 = pow(s12*s23 - s1*s3 - s12*t1 + s3*t1 + s1*t2 - s23*t2,2);
    double r4 = shat*shat*(s2*s2 + pow(t1 - t2,2) - 2.0*s2*(t1 + t2));
    double r5 = shat*(2.0*(s3*(s1 - t1)*(-s2 + t1) + 
                (s2*(s23 - 2.0*t1) + (s23 + s3)*t1 + 
                s1*(s2 - 2.0*s23 + s3 + t1))*t2 - (s1 + s23)*t2*t2 + 
                s12*(s2*(-s23 + t1) + t1*(-2.0*s3 - t1 + t2) +
                s23*(t1 + t2))));
    double r2 = (r3+r4+r5) / 16.0;
    if (r2 <= 0) return 0;

    *wtps = (*wtps)/16.0/2.0/sqrt(dlambda(shat, rma2, rmb2))/sqrt(-r2)/pow(2.0*M_PI,5);

    // p23 ---> p2 p3
    momentum_t *pres = new momentum_t[5];
    size_t size_p = sizeof(momentum_t);
    memcpy(&pres[0], &pa, size_p);
    memcpy(&pres[1], &pb, size_p);
    memcpy(&pres[2], &p1, size_p);
    memcpy(&pres[3], &p2, size_p);
    memcpy(&pres[4], &p3, size_p);

    // Save all pa, pb, p1, p2, p3 to pres and return
    return pres;
}

void makePs2cm_nj(const double s, const double cos12, const double phi, const double s1, const double s2, momentum_t *p1, momentum_t *p2) {
    double sin12 =  sqrt(1.0 - pow(cos12,2));
    double sph   =  sin(phi);
    double cph   =  cos(phi);
    double smin  =  s1 + s2;
    double roots =  sqrt(s);

    if (s < smin) std::cout << "Improper call of makePs2cm_nj" << std::endl;

    double E1 = (s + s1 - s2) / 2.0 / roots;
    double E2 = (s + s2 - s1) / 2.0 / roots;
    double pp = sqrt(E1*E1 - s1);

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

double dot(momentum_t *p1, momentum_t *p2) {
    // todo: Overload momentum_t  *
    return (p1->E*p2->E - (p1->px*p2->px + p1->py*p2->py + p1->pz*p2->pz));
}

matrix unboostrest(momentum_t *pin) {
    // todo: make it part of the momentum_t future class
   double v[4];
   matrix res(boost::extents[4][4]);

   double s = dot(pin,pin);
   if (s < 0.0) {
       std::cout << "Undue boost" << std::endl;
       abort();
   }
   double roots = sqrt(s);
   double gamma = pin->E / roots;
   v[0] = -1.0;
   v[1] = -pin->px / pin->E;
   v[2] = -pin->py / pin->E;
   v[3] = -pin->pz / pin->E;
   double v2 = 0.0;
   for (int i = 1; i < 4; i++) v2 += v[i]*v[i];

   res[0][0] = gamma;
   for (int i = 1; i < 4; i++) {
       res[i][0] = -gamma*v[i];
       res[0][i] = -gamma*v[i];
       for (int j = 1; j < 4; j++)  res[i][j] = (gamma - 1.0)*v[i]*v[j]/v2;
       res[i][i] = res[i][i] + 1.0;
   }

   return res;
}

momentum_t boostWrapper(matrix *mat, momentum_t *pin) {
    // Part of momentum class!
}
