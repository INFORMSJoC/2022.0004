#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

void Shake(Solution &S, int n, double Max_delt)
{
    int i;
	for(i=0; i<n; i++)
	{
        S.x[i] += (2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0)*Max_delt;
	}
    
} // Uniformly random perturbation with a perturbation strength 'Max_delt' 

void GeneralPertubationReduce(Solution &S, int Nstep, double MaxDelt, double Sigma, double beta)
{
	int i,j;
	int step;
	double *g;
	double NormG;
	g = new double [n];

	for(step = 0; step < Nstep; step++)
	{

		for(i=0; i< N; i++)
		{
			j = 3*i;
			S.x[j]    +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+1]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+2]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		}
	    for(int m=0; m < 5; m++)
	    {
           mygrad(g, S.x, n);
           NormG = 0 ;
           for(i=0; i<n; i++) NormG += g[i]*g[i];
           NormG = sqrt(NormG);
		   for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG;
	    }

	    MaxDelt = beta* MaxDelt;
	    Sigma   = beta* Sigma;
	   // printf("my=%lf\n",myvalue(S.x,n));
	}
   // printf("perturbation ends\n") ;
	delete [] g;
} // Sequential random perturbation with a perturbation strength 'Max_delt' 

/* The following perturbation operators are only for the test */
void Pertubation(Solution &S, int Nstep, double MaxDelt, double Sigma)
{
	int i,j;
	int step;
    double theta, phi, r;
	double *g;
	double NormG;
	g = new double [n];

	for(step = 0; step < Nstep; step++)
	{
		for(i=0; i< N; i++)
		{
			j = 3*i;
            phi   = 1.0*pi*(rand()%RAND_MAX)/RAND_MAX;
	        theta = 2.0*pi*(rand()%RAND_MAX)/RAND_MAX;
	        r     = MaxDelt;
			S.x[j]    +=   r*sin(phi)*cos(theta);
			S.x[j+1]  +=   r*sin(phi)*sin(theta);
			S.x[j+2]  +=   r*cos(phi);
		}

        mygrad(g, S.x, n);
        NormG = 0 ;
        for(i=0; i<n; i++) NormG += g[i]*g[i];
        NormG = sqrt(NormG);
		for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG ;
	}
	delete [] g;
}


void PertubationMultiple(Solution &S, int Nstep, double MaxDelt, double Sigma)
{
	int i,j;
	int step;
    double theta, phi, r;
	double *g;
	double NormG;
	g = new double [n];

	for(step = 0; step < Nstep; step++)
	{
		for(i=0; i< N; i++)
		{
			j = 3*i;
            phi   = 1.0*pi*(rand()%RAND_MAX)/RAND_MAX;
	        theta = 2.0*pi*(rand()%RAND_MAX)/RAND_MAX;
	        r     = MaxDelt;
			S.x[j]    +=   r*sin(phi)*cos(theta);
			S.x[j+1]  +=   r*sin(phi)*sin(theta);
			S.x[j+2]  +=   r*cos(phi);
		}

		for(int m=0;m < 4;m++)
		{

           mygrad(g, S.x, n);
           NormG = 0 ;
           for(i=0; i<n; i++) NormG += g[i]*g[i];
           NormG = sqrt(NormG);
		   for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG ;
		
	   }
	   
       MaxDelt = 0.99* MaxDelt;
       Sigma   = 0.99* Sigma;
		
	}
	delete [] g;
}

void GeneralPertubation(Solution &S, int Nstep, double MaxDelt, double Sigma)
{
	int i,j;
	int step;
	double *g;
	double NormG;
	g = new double [n];

	for(step = 0; step < Nstep; step++)
	{
		for(i=0; i< N; i++)
		{
			j = 3*i;
			S.x[j]    +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+1]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+2]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		}

        mygrad(g, S.x, n);
        NormG = 0 ;
        for(i=0; i<n; i++) NormG += g[i]*g[i];
        NormG = sqrt(NormG);
		for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG;
	}
	delete [] g;
}


