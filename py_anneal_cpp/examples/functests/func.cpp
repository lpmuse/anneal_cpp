#define NX 1  // dim of state variable + number of parameters
#define NS 0  // number of stimuli

using namespace alglib;

void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt);
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac);

// auxtest vector field
void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt)
{
    dxdt[0] = 5*x[0];
}

//auxtest Jacobian matrixvoid func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac)
{
    Jac[0][0] = 5;
}
