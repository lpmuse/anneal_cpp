#define NX 3  // dim of state variable + number of parameters
#define NS 0  // number of stimuli

using namespace alglib;

void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt);
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac);

// lorenz63 vector field
void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt)
{
    dxdt[0] = -10.0*x[0] + 12.5*x[1] - 1.25;
    dxdt[1] = -1.0*x[1] + 0.02*(40.0*x[0] - 20.0)*(-50.0*x[2] + 28.0) + 0.5;
    dxdt[2] = -2.66666666666667*x[2] + 0.02*(40.0*x[0] - 20.0)*(50.0*x[1] - 25.0);
}

//lorenz63 Jacobian matrix
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac)
{
    Jac[0][0] = -10.0000000000000;
    Jac[0][1] = 12.5000000000000;
    Jac[0][2] = 0;
    Jac[1][0] = -40.0*x[2] + 22.4;
    Jac[1][1] = -1.00000000000000;
    Jac[1][2] = -40.0*x[0] + 20.0;
    Jac[2][0] = 40.0*x[1] - 20.0;
    Jac[2][1] = 40.0*x[0] - 20.0;
    Jac[2][2] = -2.66666666666667;
}
