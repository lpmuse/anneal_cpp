#define NX 6   	     // dim of state variable + number of parameters
#define ND 5         // dim of state variable
#define NS 1
using namespace alglib;
void func_origin(real_1d_array &x, int it, real_1d_array &sti, real_1d_array &func);
void func_DF(real_1d_array &x, int it, real_1d_array &sti, real_2d_array &Jac);

void func_origin(real_1d_array &x, int it, real_1d_array &sti, real_1d_array &func){
	//lorenz96
	//const ae_int_t D = x.length();
	//const double v = 8.17;
	int i;
	for(i=0;i<ND;i++){
		func[i] = x[(i-1+ND)%ND]*(x[(i+1)%ND]-x[(i-2+ND)%ND])-x[i]+x[ND]+sti[it][0];
		//printf("%e\n", sin(it*0.02)-sti[it]);
	}

}

void func_DF(real_1d_array &x, int it, real_1d_array &sti, real_2d_array &Jac){
	//lorenz96 Jacobian
	//const ae_int_t D = x.length();
	//const double v = 8.17;
	int i,j;
	for(i=0;i<ND;i++){
		for(j=0;j<ND;j++)
			Jac[i][j]=0;
		Jac[i][i] = -1;
		Jac[i][(i+1)%ND] = x[(i-1+ND)%ND];
		Jac[i][(i-2+ND)%ND] = - x[(i-1+ND)%ND];
		Jac[i][(i-1+ND)%ND] = x[(i+1)%ND] - x[(i-2+ND)%ND];
		Jac[i][ND]=1;
	}
	for(i=0;i<NX;i++) Jac[ND][i]=0;

}


