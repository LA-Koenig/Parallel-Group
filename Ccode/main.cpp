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
