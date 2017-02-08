#include "phaseSpace.h"

MomentumSet phaseSpace(int npar, double s_input, const cubareal x[]) {
    // Generate shat from x[1], x[2] "naively"
    double x1 = x[1];
    double x2 = x[2];

    if (UNIT_PHASE) {
        x1 = 1.0;
        x2 = 1.0;
    }

    double shat = x1*x2*s_input;
    
    
    // Generate Higgs mass and photons 1 and 2 from x[3], x[4], x[5]
    
    // generate the 2 -> 5 system from x[6], x[7], x[8], x[9]

}
