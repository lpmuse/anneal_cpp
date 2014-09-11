#define NX 6	// dim of state variable + number of parameters
#define NS 0	// number of stimulus
using namespace alglib;
void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt);
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac);


void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt){
	//lorenz96 vector field
	dxdt[0]=x[5] - x[0] + x[4]*(x[1] - x[3]);
	dxdt[1]=x[5] + x[0]*(x[2] - x[4]) - x[1];
	dxdt[2]=x[5] + x[1]*(-x[0] + x[3]) - x[2];
	dxdt[3]=x[5] + x[2]*(-x[1] + x[4]) - x[3];
	dxdt[4]=x[5] + x[3]*(x[0] - x[2]) - x[4];
	dxdt[5]=0;

}

void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac){
	//lorenz96 Jacobian matrix
	Jac[0][0]=-1;
	Jac[0][1]=x[4];
	Jac[0][2]=0;
	Jac[0][3]=-x[4];
	Jac[0][4]=x[1] - x[3];
	Jac[0][5]=1;
	Jac[1][0]=x[2] - x[4];
	Jac[1][1]=-1;
	Jac[1][2]=x[0];
	Jac[1][3]=0;
	Jac[1][4]=-x[0];
	Jac[1][5]=1;
	Jac[2][0]=-x[1];
	Jac[2][1]=-x[0] + x[3];
	Jac[2][2]=-1;
	Jac[2][3]=x[1];
	Jac[2][4]=0;
	Jac[2][5]=1;
	Jac[3][0]=0;
	Jac[3][1]=-x[2];
	Jac[3][2]=-x[1] + x[4];
	Jac[3][3]=-1;
	Jac[3][4]=x[2];
	Jac[3][5]=1;
	Jac[4][0]=x[3];
	Jac[4][1]=0;
	Jac[4][2]=-x[3];
	Jac[4][3]=x[0] - x[2];
	Jac[4][4]=-1;
	Jac[4][5]=1;
	Jac[5][0]=0;
	Jac[5][1]=0;
	Jac[5][2]=0;
	Jac[5][3]=0;
	Jac[5][4]=0;
	Jac[5][5]=0;

}
