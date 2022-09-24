#include <math.h>
#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <sys/stat.h> //for files
#include <fcntl.h> //for files
#include <string.h>
#include <cmath>

using namespace std;

//Auxiliary Functions
double* Tri(int n, double* a, double* d, double* c, double* b);
double L2Error(double* exact, double* calc);
double LInferror3D(double* exact, double* calc);
double** makeA(int nx, int ny, int a=-4, int b=1, int d=1);

//Methods
double** crankNicolson(double h, double k, int n, int m, double* initVal, double** bounds, double alpha=1.0);
double** backwardDifference(double h, double k, int n, int m, double* initVal, double** bounds, double alpha);
double*** extrapolation(double h, double k, double nx, double ny, double m, double** initVal, double alpha=1);
double*** peacemanRachford(double h, double k, double nx, double ny, double m, double** initVal, double alpha=1);



int main(int argc, char* argv[]){


	//Variables


	//test function
	if(argc>1 and (string(argv[1]) == "-hw")){
		cout << "Hello World\n";
		exit(0);
	}

	return(1);
}


//Auxiliary Functions
double* Tri(int n, double* a, double* d, double* c, double* b){
	double* x;
	double xmult;

	x = (double*)malloc(n * sizeof(double));
	for(int i=0; i<n; i++){
		x[i] = 0;
	}

	for(int i=1; i<n; i++){
		xmult = a[i-1] / d[i-1];
		d[i] -= xmult * c[i-1];
		b[i] -= xmult * b[i-1];
	}

	x[n-1] = b[n-1] / d[n-1];

	for(int i=n-2; n > -1; n--){
		x[i] = (b[i] - c[i] * x[i+1]) / d[i];
	}

	return(x);
	//make sure to free this whereever is gets passed...
}

double L2Error(double* exact, int n, double* calc){
	double sum = 0;

	for(int i=0; i<n; i++){
		sum += pow(calc[i] - exact[i], 2);
	}
	sum = pow(sum/n, 0.5);

	return(sum);

}

double LInferror3D(double** exact, double** calc, int nx, int ny){
	double max = 0;
	double temp;

	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			temp = abs(calc[i][j] - exact[i][j]);

			if(temp > max) max = temp;
		}
	}

	return(max);
}

//double** makeA(int nx, int ny, int a, int b, int d);
//deal with 3d stuff later...
//

//Methods
double** crankNicolson(double h, double k, int n, int m, double* initVal, double** bounds, double alpha){
	//variables
	double r = (k * alpha) / h*h;
	double** U;
	double* d;
	double* c;
	double* v;

	//all the malloc assignments
	U = (double**)malloc(m * sizeof(double*));
	for(int i=0; i<m; i++){
		U[i] = (double*)malloc(n * sizeof(double));
		for(int j=0; j<n; j++){
			U[i][j] = 0;
		}
	}

	d = (double*)malloc((n-2) * sizeof(double));
	c = (double*)malloc((n-2) * sizeof(double));
	v = (double*)malloc((n-2) * sizeof(double));
	for(int i=0; i<(n-2); i++){
		d[i] = 1 + r;
		c[i] = -r / 2;
		v[i] = 0;
	}

	for(int i=0; i<n; i++){
		U[0][i] = initVal[i];
	}

	for(int i=0; i<m; i++){
		U[i][0] = bounds[i][0];
		U[i][n-1] = bounds[i][1];
	}
	
	
	//The method
	for(int j=0; j<(m-1); j++){
		for(int i=0; i<(n-2); i++){
			d[i] = 1 + r;
		}
		v[0] = (1-r) * U[j][1] + (r/2) * (U[j][2] + U[j][0] + U[j+1][0]);
		for(int i=1; i<(n-3); i++){
			v[i] = (r/2) * (U[j][i] + U[j][i+2]) + (1-r) * U[j][i+1];
		}
		v[n-3] = (r/2) * (U[j][n-4] + U[j+1][n-1] + U[j][n-1]) + (1-r) * U[j][n-3];

		v = Tri(n-2, c, d, c, v);

		for(int i=1; i<(n-1); i++){
			U[j+1][i] = v[i-1];
		}
	}

	free(d); free(c); free(v);

	return(U);
	//don't forget to free this up	
}


double** backwardDifference(double h, double k, int n, int m, double* initVal, double** bounds, double alpha){
	double r = (k * alpha) / h*h;
	double** U;
	double* d;
	double* c;
	double* v;

	//all the malloc assignments
	U = (double**)malloc(m * sizeof(double*));
	for(int i=0; i<m; i++){
		U[i] = (double*)malloc(n * sizeof(double));
		for(int j=0; j<n; j++){
			U[i][j] = 0;
		}
	}

	d = (double*)malloc((n-2) * sizeof(double));
	c = (double*)malloc((n-2) * sizeof(double));
	v = (double*)malloc((n-2) * sizeof(double));
	for(int i=0; i<(n-2); i++){
		d[i] = 1 + 2 * r;
		c[i] = -r;
		v[i] = 0;
	}

	for(int i=0; i<n; i++){
		U[0][i] = initVal[i];
	}

	for(int i=0; i<m; i++){
		U[i][0] = bounds[i][0];
		U[i][n-1] = bounds[i][1];
	}

	//The thing
	for(int j=0; j<(m-1); j++){
		for(int i=0; i<(n-2); i++){
			d[i] = 1 + 2 * r;
			v[i] = U[j][i+1];
			if(i==0){
			       	v[i] += r * U[j][i+1];
			}
			else if(i==n-3){
				v[i] += r * U[j+1][n-1];
			}
		}

		v = Tri(n-2, c, d, c, v);

		for(int i=0; i<(n-1); i++){
			U[j+1][i] = v[i-1];
		}
	}

	free(d); free(c); free(v);
	return(U);
}

double** extrapolation(double h, double k, int n, int m, double* initVal, double** bounds, double alpha=1){
	double r = (k * alpha) / h*h;
	double** U;

	if(r > 0.5){
		cout << "ERROR: r too large\n";
		return(U);
		//this is bad. Do better error catching and no extra return statements better
	}

	//all the malloc assignments
	U = (double**)malloc(m * sizeof(double*));
	for(int i=0; i<m; i++){
		U[i] = (double*)malloc(n * sizeof(double));
		for(int j=0; j<n; j++){
			U[i][j] = 0;
		}
	}


	for(int i=0; i<n; i++){
		U[0][i] = initVal[i];
	}

	for(int i=0; i<m; i++){
		U[i][0] = bounds[i][0];
		U[i][n-1] = bounds[i][1];
	}
	
	for(int i=0; i<m; i++){
		for(int j=0; j<(n-1); j++){
			U[i][j] = r * U[i-1][j-1] + (1 - 2*r) * U[i-1][j] + r * U[i-1][j+1];
		}
	}

	return(U);
}

//double*** peacemanRachford(double h, double k, double nx, double ny, double m, double** initVal, double alpha=1);


