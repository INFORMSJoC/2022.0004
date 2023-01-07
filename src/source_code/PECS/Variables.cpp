#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"

char * Output_File_Name;

INT N;     // number of circles
INT n;     // number of continous variables n=2*N
INT SR;
double R, R_min, W, alpha;
double final_time, starting_time, Time_limit;

Solution CS;
Solution OS;
Solution NS;
Solution BestS;
Solution CBestS;
Solution GlobalS;
Fsolution FS;
Fsolution FS1;

int **Adjacent;
int *NNeighbors;

void AssignMemery()
{
	 int i;
	 CS.x = (double *) malloc (n*sizeof (double)) ; CS.f = 0;
	 CS.dist =  (double *) malloc (N*sizeof (double)) ;
	 CS.id   =  (int *)malloc (N*sizeof (int)) ;
	 
     OS.x = (double *) malloc (n*sizeof (double)) ; OS.f = 0;
   	 OS.dist =  (double *) malloc (N*sizeof (double)) ;
	 OS.id   =  (int *)malloc (N*sizeof (int)) ;
     
     NS.x    =  (double *) malloc (n*sizeof (double)) ; NS.f = 0;
     NS.dist =  (double *) malloc (N*sizeof (double)) ;
	 NS.id   =  (int *)malloc (N*sizeof (int)) ;
     
     BestS.x = (double *) malloc (n*sizeof (double)) ; BestS.f = 0;
     BestS.dist =  (double *) malloc (N*sizeof (double)) ;
	 BestS.id   =  (int *)malloc (N*sizeof (int)) ;
     
     CBestS.x = (double *) malloc (n*sizeof (double)) ; CBestS.f = 0;
     CBestS.dist =  (double *) malloc (N*sizeof (double)) ;
	 CBestS.id   =  (int *)malloc (N*sizeof (int)) ;
	 
     GlobalS.x = (double *) malloc (n*sizeof (double)) ; GlobalS.f = 0;
     GlobalS.dist =  (double *) malloc (N*sizeof (double)) ;
	 GlobalS.id   =  (int *) malloc (N*sizeof (int)) ;
     
     FS.x  = (double *) malloc ((n+1)*sizeof (double)); FS.L = 0;
     FS1.x = (double *) malloc ((n+1)*sizeof (double)); FS1.L = 0;
     Adjacent =  new int *[N];
	 for(i=0; i<N; i++) Adjacent[i] = new int [MaxNN];
	 NNeighbors = new int [N];
}

void initial(Solution &S, INT n, double IR)
{
	INT I,i;
	double theta, r;
	for (I = 0; I < n/2; I++)
	{
		i = 2*I;
		S.x[i]   = IR*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
		S.x[i+1] = IR*(2.0*(rand()%RAND_MAX)/RAND_MAX - 1.0);
	}

}

void Find_Neighbor(double x[], double Dcut)
{
    int I, J, i, j;
	double x_i, y_i, x_j, y_j, dist;
    for(I = 0; I < N; I ++) NNeighbors[I] = 0;

	for(i=0; i<N; i++)
	   for(j=0; j<MaxNN; j++) Adjacent[i][j] = -1;

    for(I = 0; I < N; I ++)
    {
	  i = 2*I;
	  x_i = x[i];
	  y_i = x[i+1];

	  for(J = I+1; J < N; J ++)
	  {

		j   = 2*J;
		x_j = x[j];
		y_j = x[j+1];

		dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j));
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
    sprintf(buff,"%d.txt",n/2);
    fp=fopen(buff,"a+");
    fprintf(fp,"%d   %.10f\n", n/2, 2.0*S.CurrentR);
    for(i=0;i<n/2;i++)
    {
	   int I = 2*i;
	   fprintf(fp,"%15.15f    %15.15f\n", S.x[I],S.x[I+1]);
	}
	fclose(fp);
}

void Out_results(double best , double ave,  double worst, int sr, double AvgTime, double deviation, char *filename, int N)
{
    int i;
	FILE *fp;
	char buff[80];
    sprintf(buff,"%s",filename);
    fp = fopen(buff,"a+");
    fprintf(fp,"%d   %.10f   %.10f   %.10f  %d  %e   %.2f\n", N, 2.0*best, 2.0*ave, 2.0*worst, SR, deviation, AvgTime);
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
