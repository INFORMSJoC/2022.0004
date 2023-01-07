#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

char * Output_File_Name;

INT N;     // number of spheres
INT n;     // number of continous variables n=3*N
INT SR;
double R, R_min, alpha;
double final_time, starting_time, Time_limit;


Solution CS;
Solution NS;
Solution OS;
Solution BestS;
Solution CBestS;
Solution GlobalS;
Fsolution FS;

int **Adjacent;
int *NNeighbors;
double *p;  

void AssignMemery()
{
	 int i;
	 CS.x = (double *) malloc (n*sizeof (double)) ; CS.f = 0;
	 NS.x = (double *) malloc (n*sizeof (double)) ; NS.f = 0;
     OS.x = (double *) malloc (n*sizeof (double)) ; OS.f = 0;
     BestS.x = (double *) malloc (n*sizeof (double)) ; BestS.f = 0;
     CBestS.x = (double *) malloc (n*sizeof (double)) ; CBestS.f = 0;
     GlobalS.x = (double *) malloc (n*sizeof (double)) ; GlobalS.f = 0;
     FS.x = (double *) malloc ((n+1)*sizeof (double)); FS.L = 0;
     Adjacent =  new int *[N];
	 for(i=0; i<N; i++) Adjacent[i] = new int [MaxNN];
	 NNeighbors = new int [N];
	 p = new double [n];
}

void initial(Solution &S, INT n, double IR)
{
	INT I,i;
	double theta, r;
	for (I = 0; I < n/3; I++)
	{
		i = 3*I;
		S.x[i]   = IR*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		S.x[i+1] = IR*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
 	    S.x[i+2] = IR*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
	}

}// generate randomly an initial solution 

void Find_Neighbor(double x[], double Dcut)
{
    int I, J, i, j;
	double x_i, y_i, z_i, x_j, y_j, z_j, dist;
    for(I = 0; I < N; I ++) NNeighbors[I] = 0;

	for(i=0; i<N; i++)
	   for(j=0; j<MaxNN; j++) Adjacent[i][j] = -1;

    for(I = 0; I < N; I ++)
    {
	  i = 3*I;
	  x_i = x[i];
	  y_i = x[i+1];
	  z_i = x[i+2];

	  for(J = I+1; J < N; J ++)
	  {

		j   = 3*J;
		x_j = x[j];
		y_j = x[j+1];
		z_j = x[j+2];

		dist = sqrt((x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) + (z_i-z_j)*(z_i-z_j) );
		if(dist < Dcut)
		{
			Adjacent[I][NNeighbors[I]] = J;
			NNeighbors[I]++;
		}
	  }

	}
}

double Deviation(double arr[],int n)
{
    int i;
    double sum = 0,tmp = 0, x_avg;
    for(i = 0; i < n; ++i) sum += arr[i];
    x_avg = sum / n;
    for(i = 0; i < n; ++i)  tmp += (arr[i] - x_avg)*(arr[i] - x_avg);
    return sqrt(tmp/n);
} // standard deviation

void Outputing(Solution &S, int n)
{
    int i;
    double x,y,z;
	FILE *fp;
	char buff[80];
    sprintf(buff,"%d.txt",N);
    fp=fopen(buff,"a+");
    fprintf(fp,"%d   %.12f\n", n/3, 2.0*S.CurrentR);
    for(i=0;i<n/3;i++)
    {
	   int I = 3*i;
	   fprintf(fp,"%15.12f    %15.12f    %15.12f\n", S.x[I], S.x[I+1], S.x[I+2]);
	}
	fclose(fp);
}


void Out_results1(double RR[], double TT[], double RB, int Nruns, int N)
{
    int i;
	FILE *fp;
	char buff[80];
    sprintf(buff,"%d%s",N,"r.txt");
    fp = fopen(buff,"a+");
    for(i=0;i<Nruns;i++)
    fprintf(fp,"%d    %.10f   %.10f   %.10f\n", N, 2.0*RB, 2.0*RR[i], TT[i]);
	fclose(fp);
}

void Out_results(double best , double ave,  double worst, int sr, double AvgTime, double deviation, char *filename, int N)
{
    int i;
	FILE *fp;
	char buff[80];
    sprintf(buff,"%s",filename);
    fp = fopen(buff,"a+");
    fprintf(fp,"%d   %.12f   %.12f   %.12f   %d   %e   %.2f\n", N, 2.0*best, 2.0*ave, 2.0*worst, SR, deviation, AvgTime);
	fclose(fp);
}
