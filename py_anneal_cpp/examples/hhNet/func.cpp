#define NX 4  // dim of state variable + number of parameters
#define NS 1  // number of stimuli

using namespace alglib;

void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt);
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac);

// hhNet vector field
void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt)
{
    dxdt[0] = 1.0*sti[it][0] - 0.3*x[0] + 120.0*x[2]*pow(x[1], 3)*(-x[0] + 45.0) + 36.0*pow(x[3], 4)*(-x[0] - 82.0) - 17.8161;
    dxdt[1] = -0.0818723028574019*x[1]*exp(-0.0555555555555556*x[0]) - 0.1*(-x[1] + 1)/(-0.05*x[0] - 1.25);
    dxdt[2] = -1.0*x[2]/(0.0183156388887342*exp(-0.1*x[0]) + 1.0) + 0.0021138168395623*(-x[2] + 1)*exp(-0.05*x[0]);
    dxdt[3] = -0.0521077524598135*x[3]*exp(-0.0125*x[0]) - 0.01*(-x[3] + 1)/(-0.05*x[0] - 2.0);
}

//hhNet Jacobian matrix
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac)
{
    Jac[0][0] = -120.0*x[2]*pow(x[1], 3) - 36.0*pow(x[3], 4) - 0.3;
    Jac[0][1] = 360.0*x[2]*pow(x[1], 2)*(-x[0] + 45.0);
    Jac[0][2] = 120.0*pow(x[1], 3)*(-x[0] + 45.0);
    Jac[0][3] = 144.0*pow(x[3], 3)*(-x[0] - 82.0);
    Jac[1][0] = 0.00454846126985566*x[1]*exp(-0.0555555555555556*x[0]) - 0.005*(-x[1] + 1)/pow(-0.05*x[0] - 1.25, 2);
    Jac[1][1] = -0.0818723028574019*exp(-0.0555555555555556*x[0]) + 0.1/(-0.05*x[0] - 1.25);
    Jac[1][2] = 0;
    Jac[1][3] = 0;
    Jac[2][0] = -0.00183156388887342*x[2]*exp(-0.1*x[0])/pow(0.0183156388887342*exp(-0.1*x[0]) + 1.0, 2) - 0.000105690841978115*(-x[2] + 1)*exp(-0.05*x[0]);
    Jac[2][1] = 0;
    Jac[2][2] = -0.0021138168395623*exp(-0.05*x[0]) - 1.0/(0.0183156388887342*exp(-0.1*x[0]) + 1.0);
    Jac[2][3] = 0;
    Jac[3][0] = 0.000651346905747669*x[3]*exp(-0.0125*x[0]) - 0.0005*(-x[3] + 1)/pow(-0.05*x[0] - 2.0, 2);
    Jac[3][1] = 0;
    Jac[3][2] = 0;
    Jac[3][3] = -0.0521077524598135*exp(-0.0125*x[0]) + 0.01/(-0.05*x[0] - 2.0);
}
