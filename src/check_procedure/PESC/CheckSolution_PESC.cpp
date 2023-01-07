#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

using namespace std; 

typedef struct Configuration
{
	double *x;
	double CurrentR;
	double e;
} Solution;


Solution CS; 
int N, n;     /*  n = 3*N, where N is the number of unit spheres to be packed, and n is the number of continuous variables */  
double L, R;  /*  L is the size of container, R = L/2  */  


/*  Read the coordinate of circles from the solution file  */ 
void Read_Solution(Solution &S)
{
     int i; 
	 int k;
	 double x, y, z;
	 
	 ifstream FIC;
    
	 char buff[80]; 
     sprintf(buff,"%d_CubeSol.txt",N); 
     FIC.open(buff);
     
     if ( FIC.fail() )
     {
           cout << "### Erreur open, File_Name " << buff << endl;
           exit(0);
     } 
	 
     FIC >> k >> L;
     R = L/2.0; 
    
	 printf("N = %d  L = %15.15lf\n", k, L); 
	 
	 i= 0;
     while ( ! FIC.eof() && i<N)
     {
        FIC >> x >> y >> z;
        S.x[3*i]   = x; 
        S.x[3*i+1] = y;
        S.x[3*i+2] = z; 
	    i++;
     }	 
     
     S.CurrentR = R; 
     printf("______________________________________________________________\n"); 
     for(int i=0;i<N;i++) 
     printf("%d  %15.15lf  %15.15lf  %15.15lf\n",i+1, S.x[3*i],S.x[3*i+1],S.x[3*i+2]); 
	 printf("______________________________________________________________\n"); 
     FIC.close();
}  


/*  Calculate the energy of solution  */ 
double CalculateEnergy(double *x, int n)
{
    double f, dist, overlap, SumR;
    double x_i, x_j, y_i, y_j, z_i, z_j;
    int I, J, i, j;
    f = 0.0;
    // The unit sphere c_i is denoted by (x_i,y_i, z_i), and the unit sphere c_j is denoted by (x_j,y_j,z_j), 2R denotes the size of container (i.e., R = L/2).
   	for (I = 0; I < n/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

		overlap =  fabs(x_i)+ 1.0 - R ; // overlap = l_ix = max{ 0, |x_i| + 1 - L/2 }  
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(y_i)+ 1.0 - R ; // overlap = l_iy = max{ 0, |y_i| + 1 - L/2 } 
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(z_i)+ 1.0 - R ; // overlap = l_iz = max{ 0, |z_i| + 1 - L/2 }  
        if(overlap > 0.0) f += (overlap*overlap);
	}

	for (I = 0; I < n/3; I++)
	{
		i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

		for (J = I+1; J < n/3; J++)
		{
			j = 3*J;
			x_j = x[j];
			y_j = x[j+1];
			z_j = x[j+2];
			SumR =  2.0;  // SumR is the diameter of unit spheres 

			if(fabs(x_i - x_j) > SumR || fabs(y_i - y_j) > SumR || fabs(z_i - z_j) > SumR) continue;
			dist = sqrt((x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) + (z_i -z_j)*(z_i - z_j)); // dist = ||c_i - c_j ||
			overlap = SumR - dist ;  // overlap = l_ij = max{0, 2 - ||c_i - c_j ||}
			if(overlap > 0.0) f += (overlap*overlap);
		}

	}

    return (f) ;
}


int main (int argc, char *argv[])
{
    int i ;

	N = atoi(argv[1]); // the number of circles
	n = 3*N;           // number of variables
	CS.x = (double *) malloc (n*sizeof (double)) ; 
	CS.CurrentR = 0;   
	CS.e = 0; 
	
    Read_Solution(CS);  // Read the coordinates of circles from the solution file 
    CS.e = CalculateEnergy(CS.x,n); // calculate the energy of solution
    printf("Energy of Solution  =  %e \n", CS.e);
    
    return 0;
}

