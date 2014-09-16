// header file for minAzero, containing problem information

#define NT 100  // number of time steps
#define DT 0.02  // time step size
#define NMEA 1  // number of measurements
#define NPATH 1  // number of paths
#define BETA_MAX 10  // maximal beta
#define MULT_BASE 2  // base of the exponential expression MULT_BASE^BETA
#define RM 10000  // Rm
#define RF0 0.01  // initial Rf

#include "func.cpp"

// Set data, stimulus file names below
char data_fname[] = "twin_data.dat";
char stim_fname[] = "stimulus.dat";
