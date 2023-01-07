
#ifndef MaxNV
#define MaxNV 10000
#define MaxNN 400
#endif // the maximum number of varaibles in the objective function

#ifndef pi
#define pi 3.1415926535898
#endif // the maximum number of varaibles in the objective function

typedef struct Configuration
{
	double f;
	double *x;
	double CurrentR;
}Solution;

typedef struct FlexibleConfiguration
{
	double L;
	double *x;
} Fsolution;



extern INT n,N,SR;
extern Solution CS,NS, OS, BestS, CBestS, GlobalS;
extern Fsolution FS;
extern char * Output_File_Name;
extern double final_time, starting_time, Time_limit;
extern double R, R_min, W, alpha;
extern int **Adjacent;
extern int *NNeighbors;

extern double *p;
double myvalue
(
    double   *x,
    INT       n
) ;

void mygrad
(
    double    *g,
    double    *x,
    INT        n
) ;

double myvalueNN
(
    double   *x,
    INT       n
) ;

void mygradNN
(
    double    *g,
    double    *x,
    INT        n
) ;

double myvalueRR
(
    double   *x,
    INT       n
);
void mygradRR
(
    double    *g,
    double    *x,
    INT        n
);

double Deviation(double arr[],int n);
double AdjustRadius(Fsolution &S);
void Outputing(Solution &S, int n) ;
void Out_results(double best , double ave,  double worst, int sr, double AvgTime, double deviation, char *filename, int N) ;
void Out_results1(double RR[], double TT[], double RB, int Nruns, int N);  
void AssignMemery( );
void initial(Solution &S, INT n, double IR);
void LocalSearch(Solution &S);
void Shake(Solution &S, int n, double Max_delt);
void Pertubation(Solution &S, int Nstep, double MaxDelt, double Sigma);
void GeneralPertubation(Solution &S, int Nstep, double MaxDelt, double Sigma);
void PertubationMultiple(Solution &S, int Nstep, double MaxDelt, double Sigma)  ;
void GeneralPertubationReduce(Solution &S, int Nstep, double MaxDelt, double Sigma, double beta);
void Find_Neighbor(double x[], double Dcut);
void Thredshold_Accepting(double CR, int Flag);
void Thredshold_Search(Solution &S, double CR, int Flag)  ;
void TA_ECP();
