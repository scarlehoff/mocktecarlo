#include "MomentumSet.h"
#include "MatrixElement.h"
#include <iostream>

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

double integratedDipolesC0g1(MomentumSet *pset, double scale) {
    double rtree = C0g0WFH(pset);
    // Insertion I
    double I = 10.0 - 7.0*pow(M_PI,2)/6.0;
    return 0.0;
}
double integratedSplitting(MomentumSet *pset, const double scale, const int ix) {
    double rtree = C0g0WFH(pset);
    // Insertion I
    int i1, i2, i3, i4;
    pset->getID(&i1, &i2, &i3, &i4);
    // logs
    double rs2 = pow(scale, 2);
    double dls14 = -log(fabs(pset->s(i1,i4)/rs2));
    double dls23 = -log(fabs(pset->s(i2,i3)/rs2));

    if (ix == 1) {
    }

    double total = 0.0;

    return total;
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

