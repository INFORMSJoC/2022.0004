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
  
} // Uniformly random perturbation in which Max_delt is the perturbation strength


void GeneralPertubationReduce(Solution &S, int Nstep, double MaxDelt, double Sigma, double beta)
{
	int i,j;
	int step;
	double *g;
	double NormG;
	g = new double [n];
  //  printf("perturbation = %lf\n", MaxDelt) ;
	for(step = 0; step < Nstep; step++)
	{

		for(i=0; i< N; i++)
		{
			j = 2*i;
			S.x[j]    +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+1]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		}
	    Find_Neighbor(S.x,4.0);  // find the neighbors of each circle	
	    for(int m=0; m<5; m++)
	    {
           mygradNN(g, S.x, n);
           NormG = 0 ;
           for(i=0; i<n; i++) NormG += g[i]*g[i];
           NormG = sqrt(NormG);
		   for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG;
	    }
	    
	    MaxDelt = beta* MaxDelt;
	    Sigma   = beta* Sigma;
	  //  printf("my=%lf\n",myvalue(S.x,n));
	}
 //   printf("perturbation end = %lf\n", MaxDelt) ;
	delete [] g;
}  // Sequential Random Perturbation


/*  The following perturbation procedures are only for the test  */
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
			j = 2*i;
	        theta = 2.0*pi*(rand()%RAND_MAX)/RAND_MAX;
	        r     = MaxDelt;
			S.x[j]    +=   r*cos(theta);
			S.x[j+1]  +=   r*sin(theta);
		}

        mygrad(g, S.x, n);
        NormG = 0 ;
        for(i=0; i<n; i++) NormG += g[i]*g[i];
        NormG = sqrt(NormG);
		for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG ;
	//	printf("my=%lf\n",myvalue(S.x,n));
	}
	delete [] g;
}// only for the test


void Pertubation1(Solution &S, int Nstep, double MaxDelt, double Sigma)
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
			j = 2*i;
	        theta = 2.0*pi*(rand()%RAND_MAX)/RAND_MAX;
	        r     = MaxDelt;
			S.x[j]    +=   r*cos(theta);
			S.x[j+1]  +=   r*sin(theta);
		}
	   for(int m=0;m<2;m++)
	   {

          mygrad(g, S.x, n);
          NormG = 0 ;
          for(i=0; i<n; i++) NormG += g[i]*g[i];
          NormG = sqrt(NormG);
		  for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG ;
	  }
     // printf("my=%lf\n",myvalue(S.x,n));
	}
  //  printf("perturbation ends\n") ;
	delete [] g;
} // only for the test

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
			j = 2*i;
			S.x[j]    +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+1]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		}
	    for(int m=0;m<3;m++)
	    {
           mygrad(g, S.x, n);
           NormG = 0 ;
           for(i=0; i<n; i++) NormG += g[i]*g[i];
           NormG = sqrt(NormG);
		   for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG;
	    }
	 //   printf("my=%lf\n",myvalue(S.x,n));
	}
   // printf("perturbation ends\n") ;
	delete [] g;
} // only for the test 

void GeneralPertubationReduce1(Solution &S, int Nstep, double MaxDelt, double Sigma, double beta)
{
	int i,j;
	int step;
	double *g;
	double NormG;
	g = new double [n];
 //   printf("perturbation = %lf\n", MaxDelt) ;
	for(step = 0; step < Nstep; step++)
	{
		for(i=0; i< N; i++)
		{
			j = 2*i;
			S.x[j]    +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+1]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		}
	    for(int m=0; m<5; m++)
	    {
           mygrad(g, S.x, n);
           NormG = 0 ;
           for(i=0; i<n; i++) NormG += g[i]*g[i];
           NormG = sqrt(NormG);
		   for(i=0; i<n; i++) S.x[i] -= Sigma*g[i]/NormG;
	    }
	    
	    MaxDelt = beta* MaxDelt;
	    Sigma   = beta* Sigma;
	 //   printf("my=%lf\n",myvalue(S.x,n));
	}
  //  printf("perturbation end = %lf\n", MaxDelt) ;
	delete [] g;
}  // only for the test



void GeneralPertubation1(Solution &S, int Nstep, double MaxDelt)
{
	int i,j;
	int step;
	double *g;
	double NormG;
	double alpha;
	g = new double [n];

	for(int step = 0; step < Nstep; step++)
	{
		for(int i=0; i< N; i++)
		{
			j = 2*i;
			S.x[j]    +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
			S.x[j+1]  +=   MaxDelt*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		}
		
		for(int m=0;m<15;m++)
		{

          mygrad(g, S.x, n);
          NormG = 0 ;
          for(int i=0; i<n; i++) NormG += g[i]*g[i];
          NormG = sqrt(NormG);
	      for(int i=0; i<n; i++) S.x[i] -= 0.01*g[i]/NormG;
	   }
	   // printf("f=%lf\n", myvalue(S.x,n));
	}
	
	delete [] g;
}  // only for the test

