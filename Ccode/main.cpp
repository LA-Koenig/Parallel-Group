#include <math.h>
#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <sys/stat.h> //for files
#include <fcntl.h> //for files
#include <string.h>

using namespace std;

//Auxiliary Functions
double* Tri(int n, double* a, double* d, double* c, double* b);
double L2Error(double* exact, double* calc);
double LInferror3D(double* exact, double* calc);
double** makeA(int nx, int ny, int a=-4, int b=1, int d=1);

//Methods
double*** crankNicolson(double h, double k, double nx, double ny, double m, double** initVal, double alpha=1);
double*** backwardDifference(double h, double k, double nx, double ny, double m, double** initVal, double alpha=1);
double*** extrapolation(double h, double k, double nx, double ny, double m, double** initVal, double alpha=1);
double*** peacemanRachford(double h, double k, double nx, double ny, double m, double** initVal, double alpha=1);



int main(int argc, char* argv[]){


	//Variables
	char buffer[100];
	void *handle;
	char *error;
	void (*func)(char *, int);

	//test function
	if(argc>1 and (string(argv[1]) == "-hw")){
		cout << "Hello World\n";
		exit(0);
	}

	return(1);
}
