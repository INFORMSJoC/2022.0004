#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include "cg_user.h"
#include "Common.h"
using namespace std;

int main (int argc, char *argv[])
{
    INT i ;
	INT runs ;
	double Results[100] ;
	double Times[100] ;  
	
	double f_best, f_avg, f_worst, sigma, AvgTime ;
	
	srand(time(NULL));
	N          = atoi(argv[1]);
	runs       = atoi(argv[2]);
	Time_limit = atof(argv[3]);
	
	if(argc < 4) { printf("error, too few input parameters, the program needs 3 input parameters !");    return -1;  }
	if(argc > 4) { printf("error, too many input parameters, the program needs 3 input parameters !");   return -1;  }
	
    if(N <= 10) 
    {
    	printf("Please input an integer N > 10 (number of unit spheres) for the first parameter."); 
    	return (0);  
	}
	

    Output_File_Name = "ComputationalResultsCube.txt"; 
	n = 3*N;        // number of variables
	AssignMemery(); // allocate space for solution  

	f_best = 1.0e8;
	f_avg  = 0.0;
	f_worst= 0.0;
	AvgTime= 0.0;
	SR = 0;

	for(i = 0; i < runs; i ++)
	{
        TA_ECP();
        Results[i] = CBestS.CurrentR;
        f_avg += CBestS.CurrentR;
        if(CBestS.CurrentR < f_best - 1.0e-10)
		{
		   SR = 1;
		   f_best = CBestS.CurrentR;
		   for(int j=0;j<n;j++)  GlobalS.x[j] = CBestS.x[j];
	       GlobalS.CurrentR = CBestS.CurrentR;
        }
		else if(fabs(f_best - CBestS.CurrentR) < 1.0e-10)
		{
			SR ++ ;
		}
        if(CBestS.CurrentR > f_worst) f_worst = CBestS.CurrentR;
        Times[i] = final_time; 
        AvgTime += final_time;
	}

	f_avg /= runs;
	AvgTime /= runs;
	sigma = Deviation(Results,runs);

    R= GlobalS.CurrentR;
    Find_Neighbor(GlobalS.x,3.0);
    cg_descent (GlobalS.x, n, NULL, NULL, 1.e-15, myvalueNN, mygradNN, NULL, NULL) ; // Local minimization only considering neighbors
    GlobalS.f = myvalue(GlobalS.x,n); 
    
    Outputing(GlobalS, n); 
    Out_results(f_best, f_avg,  f_worst, SR, AvgTime, sigma, Output_File_Name, N) ;
    return 0; 
}




