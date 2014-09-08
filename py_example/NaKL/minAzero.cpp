//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "func.cpp"

#define NT 160       // number of time steps
#define DT 0.04     // time step size
#define NMEA 1       // number of measurements
#define NPATH 100    // number of paths
#define NBETA 30     // maximal beta
//using namespace alglib;

real_2d_array Ydata;
real_2d_array stimulus;

void readdata(real_2d_array &data){
	FILE *fp;
	fp = fopen("./twin_data.dat","r");
	int i,j;
	for(i=0;i<NT;i++)
		for(j=0;j<NX;j++)
			fscanf(fp,"%lf", &data[i][j]);
}

void readstimulus(real_2d_array &data){
	FILE *fp;
	fp = fopen("./stimulus.dat","r");
	int i,j;
	for(i=0;i<NT;i++)
		for(j=0;j<NS;j++)
			fscanf(fp,"%lf", &data[i][j]);
}

void slice(real_2d_array &matrix, int i, real_1d_array &output){
	int j;
	int c = matrix.cols();
	for(j=0;j<c;j++)
		output[j] = matrix[i][j];
}
	



void action_grad(const real_1d_array &x, double &action, real_1d_array &grad, void *ptr) {
	real_2d_array XX, grad_m, test;
   	XX.setlength(NT,NX);
	grad_m.setlength(NT,NX);
	test.setlength(NT,NX);
	int i,j,k;
	for(i=0;i<NT;i++){
		for(j=0;j<NX;j++){
			XX[i][j] = x[NX*i+j];
		}
	}
	//reshapevector2matrix(x,XX);
	int beta = *(int*)ptr;
	double Rm=2, Rf=0.005*pow(2,beta),tmpd;
	real_1d_array x1,x2,f1,f2,tmpa1,tmpa2;
	x1.setlength(NX);
	x2.setlength(NX);
	f1.setlength(NX);
	f2.setlength(NX);
	tmpa1.setlength(NX);
	tmpa2.setlength(NX);
	real_2d_array J1, J2, tmpm1, tmpm2;
	J1.setlength(NX,NX);
	J2.setlength(NX,NX);
	tmpm1.setlength(NX,NX);
	tmpm2.setlength(NX,NX);

	action = 0;
	for(i=0;i<NT;i++)
		for(j=0;j<NX;j++)
			grad_m[i][j]=0;
	for(i=0;i<NMEA;i++){
		for(j=0;j<NT;j++){
			
			tmpd = (XX[j][i]-Ydata[j][i]);
			
			action = action + Rm*tmpd*tmpd;
			grad_m[j][i] = grad_m[j][i] + 2*Rm*tmpd;
		}
	}
	
	for(j=0;j<(NT-1);j++){
		slice(XX, j, x1);
		slice(XX, j+1, x2);
		func_origin(x1, j, stimulus, f1);
		func_origin(x2, j, stimulus, f2);
		
		func_DF(x1, j, stimulus, J1);
		func_DF(x2, j, stimulus, J2);
		for(i=0;i<NX;i++){
			for(k=0;k<NX;k++){
				tmpm1[i][k]=-0.5*DT*J1[i][k];
				tmpm2[i][k]=-0.5*DT*J2[i][k];
				if(i==k){
					
					tmpm1[i][k] = -1 + tmpm1[i][k];
					tmpm2[i][k] = 1 + tmpm2[i][k];
				}
			}
		}
		for(i=0;i<NX;i++){
			tmpd = x2[i]-x1[i]-0.5*DT*(f1[i]+f2[i]);
			action = action + Rf*tmpd*tmpd;
			for(k=0;k<NX;k++){
				
				grad_m[j][k] = grad_m[j][k] + 2*Rf*tmpd * tmpm1[i][k];
				grad_m[j+1][k] = grad_m[j+1][k] + 2*Rf*tmpd * tmpm2[i][k];
			}
		}
	}
	
	for(i=0;i<NT;i++){
		for(j=0;j<NX;j++){
			grad[NX*i+j] = grad_m[i][j];
		}
	}
}

int main(int argc, char **argv)
{
    //
    // using LBFGS method.
    //
	Ydata.setlength(NT,NX);
	real_1d_array X0,grad_a;
	X0.setlength(NT*NX);
	grad_a.setlength(NT*NX);
	real_2d_array result;
	result.setlength(NBETA,(3+NT*NX));
	double act;
    	int i,j;
	int beta = 1, ipath;
	void *ptr;
    	double epsg = 1e-8;
    	double epsf = 1e-8;
    	double epsx = 1e-8;
	FILE *fp_output;
	char filename[20];

   	ae_int_t maxits = 0;

	readdata(Ydata);
	if(NS!=0){
		stimulus.setlength(NT,NS);
		readstimulus(stimulus);
	}
	minlbfgsstate state;
	minlbfgsreport rep;
	for(ipath=0;ipath<NPATH;ipath++){
		for(i=0;i<NX*NT;i++) X0[i]=20*randomreal()-10;
		sprintf(filename,"path/D%d_M%d_PATH%d.dat", NX,NMEA,ipath);
		fp_output = fopen(filename,"w"); 
		for(beta=0;beta<NBETA;beta++){
			printf("ipath=%d beta=%d\n", ipath, beta);
			ptr = &beta;
			
			minlbfgscreate(1, X0, state);
			minlbfgssetcond(state, epsg, epsf, epsx, maxits);
			minlbfgsoptimize(state, action_grad, NULL, ptr);
			minlbfgsresults(state, X0, rep);
			action_grad(X0, act, grad_a, ptr);
			//printf("run here\n");
			fprintf(fp_output, "%d %d %e ", beta, int(rep.terminationtype), act);
			for(i=0;i<NX*NT;i++)
				fprintf(fp_output,"%e ", X0[i]);
			fprintf(fp_output,"\n");
		}
		fclose(fp_output);
	}
    	return 0;
}
