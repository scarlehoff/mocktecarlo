#include "MomentumSet.h"
#include "MatrixElement.h"
#include <iostream>
// dilogarithm
#include <gsl/gsl_sf_dilog.h>

double matrixElement(MomentumSet *pset) {
    switch (pset->npar) {
        case 6:
            return C0g0WFH(pset);
        case 7:
            return C1g0WFH(pset);
        case 8:
            return C2g0WFH(pset);
    }
}

// Amplitudes
// 2 -> 4
double C0g0WFH(MomentumSet *pset) {
    int i1, i2, i3, i4;
    pset->getID(&i1, &i2, &i3, &i4);
    double s1i = pset->s(i1,i4);
    double s2j = pset->s(i2,i3);
    double prop = propagatorVBF(s1i, s2j, 1);

    // Amplitude
    cplx zamp = pset->zA(i1, i2)*pset->zA(i3,i4);

    double amp = pow(abs(zamp), 2);
    return 2.0*amp*prop;
}

// 2 -> 5
double C1g0WFH(MomentumSet *pset) {
    int i1, i2, i3, i4, i5;
    pset->getID(&i1, &i2, &i3, &i4, &i5);
    double a1 = C1g0WFHs0(i1, i5, i3, i2, i4, pset);
    double a2 = C1g0WFHs0(i2, i5, i4, i1, i3, pset);
    return a1 + a2;
}

double C1g0WFHup(MomentumSet *pset) {
    int i1, i2, i3, i4, i5;
    pset->getID(&i1, &i2, &i3, &i4, &i5);
    double a1 = C1g0WFHs0(i1, i5, i3, i2, i4, pset);
    return a1;
}

double C1g0WFHdown(MomentumSet *pset) {
    int i1, i2, i3, i4, i5;
    pset->getID(&i1, &i2, &i3, &i4, &i5);
    double a2 = C1g0WFHs0(i2, i5, i4, i1, i3, pset);
    return a2;
}



double C1g0WFHs0(int i1, int i5, int i3, int i2, int i4, MomentumSet *pset) {
    double s1i = pset->s(i1,i4);
    double s1k = pset->s(i1,i5);
    double sik = pset->s(i4,i5);
    double s2j = pset->s(i2,i3);
    double q1  = s1i + s1k + sik;

    double prop = propagatorVBF(q1, s2j, 1);

    // Amplitude
    // g+
    cplx z1p   = pset->zA(i1,i4) * pset->zB(i3,i4);
    cplx z2p   = pset->zA(i1,i5) * pset->zB(i3,i5);
    cplx ztp   = pset->zA(i1,i2) / pset->zA(i4,i5) / pset->zA(i1,i5);
    cplx zampp = ztp *(z1p + z2p);
    // g-
    cplx z1m   = pset->zB(i1,i4) * pset->zA(i1,i2);
    cplx z2m   = pset->zB(i4,i5) * pset->zA(i2,i5);
    cplx ztm   = pset->zB(i4,i3) / pset->zB(i4,i5) / pset->zB(i1,i5);
    cplx zampm = ztm *(z1m + z2m);

    double amp = pow(abs(zampp),2) + pow(abs(zampm),2);

    return 2.0*amp*prop;
}

// 2 -> 6
double C2g0WFH(MomentumSet *pset) {
    int i1, i2, i3, i4, i5, i6;
    pset->getID(&i1, &i2, &i3, &i4, &i5, &i6);

    //   double nadj1 = C2g0WFHnadj(i1, i5, i4,     i2, i6, i3, pset);
    //   double nadj2 = C2g0WFHnadj(i1, i6, i4,     i2, i5, i3, pset);
    //   return nadj1 + nadj2;
    double nadj1 = 0.0;
    double nadj2 = 0.0;
    double adj11 = C2g0WFHadj (i1, i5, i6, i3, i2, i4, pset);
    double adj12 = C2g0WFHadj (i1, i6, i5, i3, i2, i4, pset);
    double adj21 = C2g0WFHadj (i2, i5, i6, i4, i1, i3, pset);
    double adj22 = C2g0WFHadj (i2, i6, i5, i4, i1, i3, pset);
    return nadj1 + nadj2 + adj11 + adj12 + adj21 + adj22;
}

double C2g0WFHnadj(int i1, int i5, int i4,     int i2, int i6, int i3, MomentumSet *pset) {
    double s1i = pset->s(i1,i4);
    double s1k = pset->s(i1,i5);
    double sik = pset->s(i4,i5);
    double s2j = pset->s(i2,i3);
    double s2l = pset->s(i2,i6);
    double sjl = pset->s(i3,i6);
    double q1  = s1i + s1k + sik;
    double q2  = s2j + s2l + sjl;

    double prop = propagatorVBF(q1, q2, 1);

    // Amplitude
    // + +
    cplx z1 = - pset->zA(i1,i5)*pset->zA(i2,i6)*pset->zB(i5,i6)
        - pset->zA(i1,i5)*pset->zA(i2,i3)*pset->zB(i5,i3)
        - pset->zA(i1,i4)*pset->zA(i2,i6)*pset->zB(i4,i6)
        - pset->zA(i1,i4)*pset->zA(i2,i3)*pset->zB(i4,i3);
    cplx z2 = pset->zA(i1,i2);
    cplx zd = pset->zA(i1,i5)*pset->zA(i5,i4)*pset->zA(i2,i6)*pset->zA(i6,i3);

    cplx zamp1p2p = (z1*z2) / zd;

    // + -
    z1 = + pset->zA(i1,i5)*pset->zA(i1,i2)*pset->zB(i3,i5)*pset->zB(i3,i2)
        + pset->zA(i1,i5)*pset->zA(i1,i6)*pset->zB(i3,i5)*pset->zB(i3,i6)
        + pset->zA(i1,i4)*pset->zA(i1,i2)*pset->zB(i3,i4)*pset->zB(i3,i2)
        + pset->zA(i1,i4)*pset->zA(i1,i6)*pset->zB(i3,i4)*pset->zB(i3,i6);
    zd = pset->zA(i1,i5)*pset->zA(i5,i4)*pset->zB(i2,i6)*pset->zB(i6,i3);

    cplx zamp1p2m = z1/zd;

    // - + 
    z1 = - pset->zA(i1,i2)*pset->zA(i2,i6)*pset->zB(i4,i1)*pset->zB(i4,i6)
        - pset->zA(i1,i2)*pset->zA(i2,i3)*pset->zB(i4,i1)*pset->zB(i4,i3)
        - pset->zA(i5,i2)*pset->zA(i2,i6)*pset->zB(i4,i5)*pset->zB(i4,i6)
        - pset->zA(i5,i2)*pset->zA(i2,i3)*pset->zB(i4,i5)*pset->zB(i4,i3);
    zd = pset->zB(i1,i5)*pset->zB(i5,i4)*pset->zA(i2,i6)*pset->zA(i6,i3);

    cplx zamp1m2p = z1/zd;

    // - -
    z1 = + pset->zA(i1,i2)*pset->zB(i4,i1)*pset->zB(i2,i3)
        + pset->zA(i1,i6)*pset->zB(i4,i1)*pset->zB(i6,i3)
        + pset->zA(i5,i2)*pset->zB(i4,i5)*pset->zB(i2,i3)
        + pset->zA(i5,i6)*pset->zB(i4,i5)*pset->zB(i6,i3);
    z2 = pset->zB(i4,i3);
    zd = pset->zB(i1,i5)*pset->zB(i5,i4)*pset->zB(i2,i6)*pset->zB(i6,i3);

    cplx zamp1m2m = (z1*z2)/zd;

    // total:

    double amp = pow(abs(zamp1p2p), 2) + pow(abs(zamp1p2m), 2)
        + pow(abs(zamp1m2p), 2) + pow(abs(zamp1m2m), 2);

    return 2.0 * amp * prop;
}

double C2g0WFHadj(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset) {
    // i5 and i6 radiated between it and i4
    double s1i = pset->s(i1,i4);
    double s1k = pset->s(i1,i5);
    double s1l = pset->s(i1,i6);
    double sik = pset->s(i4,i5);
    double sil = pset->s(i4,i6);
    double skl = pset->s(i5,i6);
    double s2j = pset->s(i2,i3);
    double q1  = s1k + s1l + s1i + skl + sik + sil;

    double prop = propagatorVBF(q1, s2j, 1);

    // Amplitude
    cplx zamp1p2p = C2g0WFHadj1p2p(i1, i5, i6, i3, i2, i4, pset);
    cplx zamp1p2m = C2g0WFHadj1p2m(i1, i5, i6, i3, i2, i4, pset);
    cplx zamp1m2p = C2g0WFHadj1m2p(i1, i5, i6, i3, i2, i4, pset);
    cplx zamp1m2m = C2g0WFHadj1m2m(i1, i5, i6, i3, i2, i4, pset);

    double amp = pow(abs(zamp1p2p), 2) + pow(abs(zamp1p2m), 2)
        + pow(abs(zamp1m2p), 2) + pow(abs(zamp1m2m), 2);

    return 2.0 * amp * prop;
}

cplx C2g0WFHadj1p2p(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset) {
    cplx z1 = pset->zA(i1, i4)*pset->zB(i4, i3);
    cplx z2 = pset->zA(i1, i5)*pset->zB(i3, i5);
    cplx z3 = pset->zA(i1,i6)*pset->zB(i3, i6);
    cplx z4 = pset->zA(i1,i2);
    cplx zd = pset->zA(i1,i5)*pset->zA(i5,i6)*pset->zA(i6,i4);
    return (+z1 -z2 -z3)*z4/zd;

}
cplx C2g0WFHadj1p2m(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset) {
    cplx z1d = + pow(pset->zA(i1,i5),2)*pset->zA(i5,i6)*pset->zB(i5,i4)*pset->sijk(i1,i6,i5)
        + pset->zA(i1,i5)*pset->zA(i1,i6)*pset->zA(i5,i6)*pset->zB(i6,i4)*pset->sijk(i1,i6,i5);
    cplx z2d = + pset->zA(i6,i1)*pow(pset->zB(i6,i4),2)*pset->zB(i6,i5)*pset->sijk(i4,i6,i5)
        - pset->zA(i5,i1)*pset->zB(i6,i4)*pset->zB(i6,i5)*pset->zB(i4,i5)*pset->sijk(i4,i6,i5);
    cplx z1 = 1.0/z1d;
    cplx z2 = 1.0/z2d;
    cplx zT = + pow(pset->zA(i1,i6),3)*pset->zA(i1,i2)*pset->zB(i1,i4)*pset->zB(i4,i3)*z1
        + pow(pset->zA(i1,i6),3)*pset->zA(i6,i2)*pset->zB(i6,i4)*pset->zB(i4,i3)*z1
        - pow(pset->zA(i1,i6),3)*pset->zA(i5,i2)*pset->zB(i4,i3)*pset->zB(i4,i5)*z1
        - pset->zA(i1,i6)*pset->zA(i1,i2)*pset->zB(i6,i3)*pow(pset->zB(i4,i5),3)*z2
        - pset->zA(i1,i4)*pset->zA(i1,i2)*pset->zB(i4,i3)*pow(pset->zB(i4,i5),3)*z2
        + pset->zA(i1,i5)*pset->zA(i1,i2)*pow(pset->zB(i4,i5),3)*pset->zB(i3,i5)*z2;
    return zT;
}
cplx C2g0WFHadj1m2p(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset) {
    cplx z1d = + pset->zA(i5,i4)*pow(pset->zB(i1,i5),2)*pset->zB(i5,i6)*pset->sijk(i1,i6,i5)
        + pset->zA(i6,i4)*pset->zB(i1,i5)*pset->zB(i1,i6)*pset->zB(i5,i6)*pset->sijk(i1,i6,i5);
    cplx z2d = - pset->zA(i5,i6)*pset->zA(i5,i4)*pset->zA(i6,i4)*pset->zB(i5,i1)*pset->sijk(i4,i6,i5)
        - pset->zA(i5,i6)*pow(pset->zA(i6,i4),2)*pset->zB(i6,i1)*pset->sijk(i4,i6,i5);
    cplx z1 = 1.0/z1d;
    cplx z2 = 1.0/z2d;
    cplx zT = - pset->zA(i1,i6)*pset->zA(i1,i2)*pow(pset->zB(i1,i6),3)*pset->zB(i6,i3)*z1
        + pset->zA(i1,i6)*pset->zA(i5,i2)*pow(pset->zB(i1,i6),2)*pset->zB(i6,i3)*pset->zB(i6,i5)*z1
        + pset->zA(i1,i5)*pset->zA(i1,i2)*pow(pset->zB(i1,i6),3)*pset->zB(i3,i5)*z1
        - pset->zA(i1,i5)*pset->zA(i1,i2)*pow(pset->zB(i1,i6),2)*pset->zB(i1,i3)*pset->zB(i6,i5)*z1
        - pset->zA(i1,i5)*pset->zA(i5,i2)*pow(pset->zB(i1,i6),2)*pset->zB(i6,i5)*pset->zB(i3,i5)*z1
        + pset->zA(i1,i5)*pset->zA(i5,i2)*pset->zB(i1,i6)*pset->zB(i1,i3)*pow(pset->zB(i6,i5),2)*z1
        + pset->zA(i1,i4)*pset->zA(i1,i2)*pow(pset->zB(i1,i6),3)*pset->zB(i3,i4)*z1
        - pset->zA(i1,i4)*pset->zA(i5,i2)*pow(pset->zB(i1,i6),2)*pset->zB(i6,i5)*pset->zB(i3,i4)*z1
        - pset->zA(i1,i2)*pow(pset->zA(i6,i5),2)*pset->zA(i5,i4)*pset->zB(i1,i6)*pset->zB(i6,i3)*z2
        - pset->zA(i1,i2)*pset->zA(i6,i5)*pow(pset->zA(i5,i4),2)*pset->zB(i1,i6)*pset->zB(i3,i4)*z2
        + pset->zA(i1,i2)*pset->zA(i6,i5)*pow(pset->zA(i5,i4),2)*pset->zB(i1,i4)*pset->zB(i6,i3)*z2
        - pset->zA(i1,i2)*pset->zA(i6,i5)*pow(pset->zB(i1,i6),2)*pset->zB(i6,i3)*pset->zB(i6,i5)*z1
        + pset->zA(i1,i2)*pow(pset->zA(i5,i4),3)*pset->zB(i1,i4)*pset->zB(i3,i4)*z2
        - pset->zA(i1,i2)*pset->zA(i5,i4)*pow(pset->zB(i1,i6),2)*pset->zB(i6,i5)*pset->zB(i3,i4)*z1
        + pow(pset->zA(i6,i5),2)*pset->zA(i5,i4)*pset->zA(i5,i2)*pset->zB(i6,i3)*pset->zB(i6,i5)*z2
        + pow(pset->zA(i6,i5),2)*pset->zA(i5,i4)*pset->zA(i4,i2)*pset->zB(i6,i3)*pset->zB(i6,i4)*z2
        + pset->zA(i6,i5)*pset->zA(i6,i2)*pow(pset->zA(i5,i4),2)*pset->zB(i6,i3)*pset->zB(i6,i4)*z2
        + pset->zA(i6,i5)*pow(pset->zA(i5,i4),2)*pset->zA(i5,i2)*pset->zB(i6,i3)*pset->zB(i5,i4)*z2
        + pset->zA(i6,i5)*pow(pset->zA(i5,i4),2)*pset->zA(i5,i2)*pset->zB(i6,i5)*pset->zB(i3,i4)*z2
        + pset->zA(i6,i5)*pow(pset->zA(i5,i4),2)*pset->zA(i4,i2)*pset->zB(i6,i4)*pset->zB(i3,i4)*z2
        + pset->zA(i6,i5)*pset->zA(i5,i2)*pset->zB(i1,i6)*pset->zB(i6,i3)*pow(pset->zB(i6,i5),2)*z1
        + pset->zA(i6,i2)*pow(pset->zA(i5,i4),3)*pset->zB(i6,i4)*pset->zB(i3,i4)*z2
        + pow(pset->zA(i5,i4),3)*pset->zA(i5,i2)*pset->zB(i3,i4)*pset->zB(i5,i4)*z2
        + pset->zA(i5,i4)*pset->zA(i5,i2)*pset->zB(i1,i6)*pow(pset->zB(i6,i5),2)*pset->zB(i3,i4)*z1;
    return zT;
}
cplx C2g0WFHadj1m2m(int i1, int i5, int i6, int i3, int i2, int i4, MomentumSet *pset) {
    cplx z1 = - pset->zA(i1, i2)*pset->zB(i1, i4) - pset->zA(i2, i5)*pset->zB(i4, i5)
        - pset->zA(i2,i6)*pset->zB(i4, i6);
    cplx z2 =   pset->zB(i4,i3);
    cplx zd =   pset->zB(i1,i5)*pset->zB(i5,i6)*pset->zB(i6,i4);
    return (z1 * z2) / zd;
}

// Virtual functions
double vertexCorrectionC0g1(MomentumSet *pset, double scale) {
    double rtree = C0g0WFH(pset);
    int i1, i2, i3, i4;
    pset->getID(&i1, &i2, &i3, &i4);
    // logs
    double rs2 = pow(scale, 2);
    double dls14 = -log(fabs(pset->s(i1,i4)/rs2));
    double dls23 = -log(fabs(pset->s(i2,i3)/rs2));
    // vertex corrections
    // cdr: cvirt = pi^2/6 - 8 when summing both legs, taking out a factor of bar(c(e))
    double cvirt = pow(M_PI,2)/12.0 - 4;
    double vc1 = -3.0/2.0*dls14 - pow(dls14,2)/2.0 + cvirt;
    double vc2 = -3.0/2.0*dls23 - pow(dls23,2)/2.0 + cvirt;
    double total = (vc1 + vc2)*rtree;
    return total;
}

double integratedDipolesC0g1(MomentumSet *pset, double scale, const int ix, const double xx1, const double xx2) {
    double x1 = xx1;
    double x2 = xx2;
    double rtree = C0g0WFH(pset);
    int i1, i2, i3, i4;
    pset->getID(&i1, &i2, &i3, &i4);
    // logs
    double rs2 = pow(scale, 2);
    double dls14 = -log(fabs(pset->s(i1,i4)/rs2));
    double dls23 = -log(fabs(pset->s(i2,i3)/rs2));
//    dls14 = 0.0 ; dls23 = 0.0;
//    dls23 = 4.3801956656574976;
//    dls14 = -0.47943487776582028;
    double total = 0.0;
    double z1 = pset->x1;
    double z2 = pset->x2;
    double pi2 = pow(M_PI, 2);
    // subtraction type
    int subtraction = 3;
    // Auxiliary definitions
//        z2 = 0.12190765265928201;
//        z1 = 6.1001586905563530E-003;
//        x2 = 0.60743837005012802;
//        x1 = 0.56845142508152535;
//        z1 = 8.5132577035339566E-003;
//        z2 = 0.51985205836236692;
//        x1 = 0.95580243388760533;
//        x2 = 0.61075662779675721;
    double omz1 = 1.0 - z1;
    double omx1 = 1.0 - x1;
    double d0x1 = 1.0/omx1;
    double d1x1 = log(omx1)/omx1;
    double omz2 = 1.0 - z2;
    double omx2 = 1.0 - x2;
    double d0x2 = 1.0/omx2;
    double d1x2 = log(omx2)/omx2;
    double a1=0., a1u=0., a1l=0.;
    double c1zl=0., c1zu=0.;
    double c1u=0., c1l=0.;
    double bx=0., cx=0.;
    if (subtraction == 0) { 
        /* Antenna subtraction */
        switch (ix) {
            case 1: // 1, 1
                // A(1)/(1-z)
                a1 = (7.0 - pi2)/4.0 ;
                a1u = a1 + 3.0*dls14/4.0 + pow(dls14,2)/2.0;
                a1l = a1 + 3.0*dls23/4.0 + pow(dls23,2)/2.0;
                total = (a1u + a1l)/omz1/omz2;
                // C(1)*ln(1-z)/(1-z)
                c1zu = (-3.0/4.0 + log(omz1)/2.0 - dls14)*log(omz1);
                c1zl = (-3.0/4.0 + log(omz2)/2.0 - dls23)*log(omz2);
                total += (c1zu+c1zl)/omz1/omz2;
                // C(1)*Dn(1-x)
                c1u = d0x1*(-3.0/4.0 - dls14) + d1x1;
                c1l = d0x2*(-3.0/4.0 - dls23) + d1x2;
                total -= c1u/omz2 + c1l/omz1;
                break;
            case 2: // 1, x2
                bx = log(x2)*( (1.0+x2)/2.0 - 1.0/omx2 );
                bx += -log(omx2)*(1.0+x2)/2.0;
                bx += (3.0 - x2 + dls23 +x2*dls23)/2.0;
                cx = d0x2*(-3.0/4.0 - dls23);
                cx += d1x2;
                total = (bx + cx)/omz1;
                break;
            case 3: // x1, 1
                bx = log(x1)*( (1.0+x1)/2.0 - 1.0/omx1 );
                bx += -log(omx1)*(1.0+x1)/2.0;
                bx += (3.0 - x1 + dls14 +x1*dls14)/2.0;
                cx = d0x1*(-3.0/4.0 - dls14);
                cx += d1x1;
                total = (bx + cx)/omz2;
                break;
        }
//        std :: cout << "ix: " << ix << " |total: " << -total << std::endl;
//        std :: cin.get();
    } else if (subtraction == 1) {
        /* CS paper, DIS example 
         * ics  = (10.0 - 7.0*pow(M_PI,2)/6.0) dd + dls
         * kqq  = 2D1 - 2D0*log(x) - (1+x)*(log(1-x)-log(x)) + (1-x) -dd*(5-pi^2)
         * rest = -3/2*D0 + 3/2*dd - Pqq*(log(x/z???) + log(muF/Q2))
        */
//        if(ix==1) std :: cout << "-------------" << std::endl;
        switch (ix) {
            case 1: // 1, 1
                // A(1)/(1-z)
                a1 = 7.0/2.0 - pi2/6.0;
                a1u = a1 + 3.0*dls14/2.0 + pow(dls14,2)/1.0;
                a1l = a1 + 3.0*dls23/2.0 + pow(dls23,2)/1.0;
                total = (a1u + a1l)/omz1/omz2;
                // C(1)*ln(1-z)/(1-z)
                c1zu = helperC1(z1);
                c1zu += -2.0*log(omz1)*dls14;
                c1zu += pow(log(omz1),2);
                c1zl = helperC1(z2);
                c1zl += -2.0*log(omz2)*dls23;
                c1zl += pow(log(omz2),2);
                total += (c1zu+c1zl)/omz1/omz2;
                // C(1)*Dn(1-x)
                c1u = d0x1*(3.0/2.0); 
                c1u += d0x1*dls14*2.0;
                c1u += d0x1*(2.0*log(z1));
                c1u += -2.0*d1x1;    
                c1l = d0x2*(3.0/2.0);
                c1l += d0x2*dls23*2.0;
                c1l += d0x2*(2.0*log(z2));
                c1l += -2.0*d1x2;
                total += c1u/omz2 + c1l/omz1;
                break;
            case 2: // 1, x2
                bx = omx2 - (1.0+x2)*log(omx2/x2);
                bx += (dls23 +x2*dls23);
                cx = d0x2*(-3./2.0);
                cx += d0x2*(-2.*log(x2));
                cx += log(z2)*( (1.0+x2) - 2.0*d0x2 );
                cx += 2.*log(x2)*d0x2*( 1.0 + pow(x2,2) );
                cx += 2.0*d1x2;
                cx += -d0x2*dls23*2.0;
                total = (bx+cx)/omz1;
                break;
            case 3: // x1, 1
                bx = omx1 - (1.0+x1)*log(omx1/x1);
                bx += (dls14 +x1*dls14);
                cx = d0x1*(-3./2.0);
                cx += d0x1*(-2.*log(x1));
                cx += log(z1)*( (1.0+x1) - 2.0*d0x1 );
                cx += 2.*log(x1)*d0x1*( 1.0 + pow(x1,2) );
                cx += 2.0*d1x1;
                cx += -d0x1*dls14*2.0;
                total = (bx+cx)/omz2;
                break;
        }
//        std :: cout << "ix: " << ix << " |total: " << total/2. << std::endl;
//        std :: cin.get();
        total = total / 2.0;
    } else if (subtraction == 2) {
        switch (ix) {
            case 1: // 1, 1
                // A(1)/(1-z)
                a1 = 7.0/2.0 - 3.0*pi2/6.0;
                a1u = a1 + 3.0*dls14/2.0 + pow(dls14,2)/1.0;
                a1l = a1 + 3.0*dls23/2.0 + pow(dls23,2)/1.0;
                total = (a1u + a1l)/omz1/omz2;
                // C(1)*ln(1-z)/(1-z)
                c1zu = -3.0/2.0*log(omz1);
                c1zu += -2.0*log(omz1)*dls14;
                c1zu += pow(log(omz1),2);
                c1zl = -3.0/2.0*log(omz2);
                c1zl += -2.0*log(omz2)*dls23;
                c1zl += pow(log(omz2),2);
                total += (c1zu+c1zl)/omz1/omz2;
                // C(1)*Dn(1-x)
                c1u = d0x1*(3.0/2.0); 
                c1u += d0x1*dls14*2.0;
                c1u += -2.0*d1x1;    
                c1l = d0x2*(3.0/2.0);
                c1l += d0x2*dls23*2.0;
                c1l += -2.0*d1x2;
                total += c1u/omz2 + c1l/omz1;
                break;
            case 2: // 1, x2
                bx = omx2 - (1.0+x2)*log(omx2);
                bx += (dls23 +x2*dls23);
                cx = d0x2*(-3./2.0);
                cx += 2.0*d1x2;
                cx += -d0x2*dls23*2.0;
                total = (bx+cx)/omz1;
                break;
            case 3: // x1, 1
                bx = omx1 - (1.0+x1)*log(omx1);
                bx += (dls14 +x1*dls14);
                cx = d0x1*(-3./2.0);
                cx += 2.0*d1x1;
                cx += -d0x1*dls14*2.0;
                total = (bx+cx)/omz2;
                break;
        }
//        std :: cout << "ix: " << ix << " |total: " << total/2. << std::endl;
//        std :: cin.get();
        total = total / 2.0;
    } else if (subtraction == 3) {
        double i1 = 0.0, i2 = 0.0, i3 = 0.0;
        double scale = 0.0;
        switch (ix) {
            case 1:
                i1  = csI0(x1, z1, ix)/omz2;
                i1 += csI0(x2, z2, ix)/omz1;
                i2  = csI1(x1, z1, ix)/omz2;
                i2 += csI1(x2, z2, ix)/omz1;
                i3  = csI2(x1, z1, ix)/omz2;
                i3 += csI2(x2, z2, ix)/omz1;
                scale = 3.0*dls14/2.0 + pow(dls14,2);
                scale += 3.0*dls23/2.0 + pow(dls23,2);
                scale += -2.0*log(omz1)*dls14;
                scale += -2.0*log(omz2)*dls23;
                scale = scale/omz1/omz2;
                scale += 2.0*(d0x1*dls14/omz2 + d0x2*dls23/omz1);
                break;
            case 2:
                i1 = csI0(x2, z2, ix)/omz1;
                i2 = csI1(x2, z2, ix)/omz1;
                i3 = csI2(x2, z2, ix)/omz1;
                scale = dls23 + x2*dls23;
                scale += -2.0*dls23*d0x2;
                scale = scale/omz1;
                break;
            case 3:
                i1 = csI0(x1, z1, ix)/omz2;
                i2 = csI1(x1, z1, ix)/omz2;
                i3 = csI2(x1, z1, ix)/omz2;
                scale = dls14 + x1*dls14;
                scale += -2.0*dls14*d0x1;
                scale = scale/omz2;
                break;
        }
        total = (i1 + i2 + i3) / 2.0;
        total += scale / 2.0;
//        std :: cout << "ix: " << ix << " |total: " << total << std::endl;
//        std :: cin.get();
    }
    return (total)*rtree;
}

// Propagators

double propagatorVBF(const double s1, const double s2, const int iboson) {
    // Only W boson to Higgs right now so iboson is irrelevant
    return propagatorWFH(s1, s2);
}

double propagatorWFH(const double s1, const double s2) { 
    double emw     = 80.398;
    double ewwidth = 2.1054;
    double stw     = 0.22264585341299603;
    double p1      = propagator(s1, emw, ewwidth);
    double p2      = propagator(s2, emw, ewwidth);
    return emw*emw/p1/p2/pow(stw,3);
}

double propagator(const double s, const double mass, const double width) {
    double t1 = pow(s - mass*mass, 2);
    double t2 = pow(mass*width, 2);
    return t1 + t2;
}

double helperC1(double x) {
    double li2 = gsl_sf_dilog(1.0 - x);
    double total = -3.0*log(1.-x)/2.0;
    total += 2.0*li2 - pow(M_PI,2)/3.0;
//    total += -log(x)*(x + pow(x,2)/2.0 + 2.0*log(1.0-x));
    total += -4.0*li2 + 4.0*pow(M_PI,2)/6.0 - x/4.*(-x+2.*(x+2.)*log(x)-4.);
    return total;
}

double csI0(const double x, const double z, const int ix) {
    double a = 0.0, b = 0.0, c = 0.0;
    double omz = 1.0-z;
    if (ix == 1) {
        a = (10.0 - 7.0*pow(M_PI,2)/6.0)/omz;
    }
    return a + b + c;
}

double csI1(const double x, const double z, const int ix) {
    double a = 0.0, b = 0.0, c = 0.0;
    double omz = 1.0-z;
    double lz  = log(omz);
    double omx = 1.0-x;
    double d0x = 1.0/omx;
    double d1x = log(omx)/omx;
    if (ix == 1) {
        a = (-3.0/2.0)/omz;
        c = -3.0*lz/omz/2.0 - pow(lz,2)/omz;
        c += -(-3.0/2.0*d0x - 2.0*d1x);
    } else {
        b = 2.0*log(2.0-x)/omx;
        c = d0x*(-3.0/2.0) + d1x*(-2.0);
    }
    return a + b + c;;
}

double csI2(const double x, const double z, const int ix) {
    double a = 0.0, b = 0.0, c = 0.0;
    double omz = 1.0-z;
    double lz  = log(omz);
    double omx = 1.0-x;
    double d0x = 1.0/omx;
    double d1x = log(omx)/omx;
    if (ix == 1) {
        a = (2.0/3.0*pow(M_PI,2) - 5)/omz;
        c = 2.0*pow(lz,2)/omz ;
        c += - 4.0*d1x;
    } else {
        b = -2.0*log(2.0-x)/omx + omx - (1.0+x)*log(omx);
        b += log(x) * (1.0 + x - 2.0/omx);
        c = 4.0*d1x;
    }
    return a + b + c;;
}
