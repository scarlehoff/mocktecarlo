#include <cuba.h>
#include <math.h>
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include "PhaseSpace.h" 
#include "MomentumSet.h" 
#include "MatrixElement.h"
#include "SubtractionTerm.h"
#include "CrossSection.h"

#define FBGEV2 389379365600
#define NC 3.0
#define AMZ 7.5563839072873536E-003

using namespace std;
using namespace LHAPDF;

int crossSection(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *pdf) {
    int n        = NPARTICLES;
    double roots = ROOTS_;
    double muR   = SCALE_; // muR = muF = mH
    double s     = pow(roots, 2);

    // Generate phase space point from cuba x[]
    MomentumSet pset = phaseSpace(n, s, x);
    if (pset.ifail) {
        f[0] = 0.0;
        return 0;
    }

    // Check whether the point goes through all cuts (Check the momentum is not null, if it is then return)
    if (UNIT_PHASE) {
        int ifail = pset.sijcuts(Y0);
        if (!ifail) {
            f[0] = pset.weight/abs(pset.s(1,3));
        } else {
            f[0] = 0.0f;
        }
    } else {
        int minjets  = MINJETS_;
        double ptcut = PTCUT_;
        double rkt = RKT_;
        if(n > 6) minjets += 1;
        if(n > 7) minjets += 1;
        int ifail = pset.apply_cuts(ptcut, rkt, minjets);

        // Define particle identities
        switch (n) {
            case 6:
                // int i1 = 2; int i2 = 1; int i3 = 3; int i4 = 4;
                pset.setID(2, 1, 3, 4);
                break;
            case 7:
                //	int i1 = 2; int i2 = 1; int i3 = 4; int i4 = 5; int i5 = 3;
                pset.setID(2, 1, 4, 5,  3);
                break;
            case 8:
                //	int i1 = 2; int i2 = 1; int i3 = 5; int i4 = 6; int i5 = 3; int i6 = 4;
                pset.setID(2, 1, 5, 6,  3, 4);
                break;
        }
        // Compute Matrix Element
        double mesq, vfact, intfact;
        if (!ifail) {
            if (VIRTUAL) {
                vfact = vertexCorrectionC0g1(&pset, muR);
                intfact = integratedDipolesC0g1(&pset, muR);
                mesq = vfact - intfact;
            } else {
                mesq = matrixElement(&pset);
            }
        } else {
            mesq = 0.0;
        }

        // Compute subtraction term if necessary
        if ( (n-4) > minjets ) {
            double subt = subtractionTerm(&pset, ptcut, rkt, minjets);
            mesq = mesq - subt;
        }

        // Compute the PDFs for a u(c) d(s) interaction
        double f1      = (*((PDF **) pdf))->xfxQ2(2, pset.x1, pow(muR,2));
        f1+= (*((PDF **) pdf))->xfxQ2(4, pset.x1, pow(muR,2));
        double f2      = (*((PDF **) pdf))->xfxQ2(1, pset.x2, pow(muR,2));
        f2+= (*((PDF **) pdf))->xfxQ2(3, pset.x2, pow(muR,2));
        double pdfval  = f1*f2/pset.x1/pset.x2;
        double alpha_s = (*((PDF **) pdf))->alphasQ2(pow(muR,2));

        // Compute flux factor for this ps point and QCD factor (cte)
        double average = (1.0/NC)*(1.0/NC)/4.0;
        double qcdborn = pow(4.0*M_PI*AMZ, 3)*NC*NC/2.0;
        double qcdfactor = qcdborn;
        if (n = 6 && VIRTUAL) {
            qcdfactor = qcdfactor*(4.0*M_PI*alpha_s)*(NC*NC-1.0)/NC ;
            qcdfactor = qcdfactor / (8.0*pow(M_PI,2));
        }
        if (n > 6) qcdfactor = qcdfactor*(4.0*M_PI*alpha_s)*(NC*NC-1.0)/NC ;
        if (n > 7) qcdfactor = qcdfactor*(4.0*M_PI*alpha_s)*NC/2.0;

        double flux = average*qcdfactor*FBGEV2/2.0/pset.s(1,2);

        // Put everything together
        f[0] = mesq*flux*pset.weight*pdfval;
    }

    return 0;
}
