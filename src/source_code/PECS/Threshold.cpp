#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "cg_user.h"
#include "Common.h"


void LocalSearch(Solution &S)
{
    cg_descent (S.x, n, NULL, NULL, 1.0e-2, myvalue, mygrad, NULL, NULL) ; // Local minimization
    S.f = myvalue(S.x,n) ;  // printf("S.f =%e\n",S.f);
    Find_Neighbor(S.x,3.2);  // find the neighbors of each circle
    cg_descent (S.x, n, NULL, NULL, 1.0e-13, myvalueNN, mygradNN, NULL, NULL) ; // Local minimization only considering neighbors
    S.f = myvalue(S.x,n) ;
  //  printf("f1= %e   fn= %e\n", myvalue(S.x,n), myvalueNN(S.x,n) );
} // Two-pahse local search method for solving the subproblems

double AdjustRadius(Fsolution &S)
{
    double pho = 5.0;
    alpha = 1.0e6;
    
    Find_Neighbor(S.x,3.3);
   	cg_descent (S.x, n+1, NULL, NULL, 1.e-12, myvalueRR, mygradRR, NULL, NULL) ;
	for(int i = 0; i < 20; i++) 
	{
        alpha = alpha*pho;
	    Find_Neighbor(S.x,3.0);
        cg_descent (S.x, n+1, NULL, NULL, 1.e-13, myvalueRR, mygradRR, NULL, NULL) ;
	   //printf("alpha =%e\n",alpha);
	}

	return S.x[n];
} //  Container adjustment procedure 

void Thredshold_Accepting(double CR, const int flag)
{
     int i, step, pp;
     INT MaxSteps  = 1000;
     INT MaxNoImprove = 400;
     INT NoImprove;
     double K, K1, K_min, K_min1, K_delt, K_delt1, sigma, T;
     double r ; 
	 double pho = 0.5, alpha = 0.6;
	 double E_diff = 1.0e-4;
	 long int Naccept = 0;
	 long int Nreject = 0;
     double delt_f;

   //printf("flag=%d in the second phase\n",flag); 
     initial(CS, n, CR - 1.0);
     R =  CR;
	 LocalSearch(CS);

     for(i=0;i<n;i++) BestS.x[i]= CS.x[i];
     BestS.f = CS.f;

	 T = 72.0;        // T is a period 
	 sigma = 2.0*pi/T; 
	 K_min = 0.2;
	 K_delt= 0.125; 

     step = 0;
     pp = 0;
     NoImprove = 0 ; 
     while(NoImprove < MaxNoImprove) // the search depth 
     {

       for(i=0;i<n;i++) OS.x[i]= CS.x[i];
       OS.f = CS.f;

       if(flag ==1)
       {
          K = K_min + K_delt*(1.0 + sin(pp*sigma));  // K is the perturbation strength 
          GeneralPertubationReduce(OS,50, K, 4*K, 0.99);  // Sequential Random Perturbation
       }
       else
       {
          Shake(OS, n, 0.8); // K (=0.8) is the perturbation strength, Uniformly random perturbation
	   }
       
       LocalSearch(OS);
       delt_f = OS.f - CS.f ;
       if(delt_f < E_diff && fabs(delt_f) > 1.0e-20 )
	   {
		  for(i=0;i<n;i++) CS.x[i]= OS.x[i];
          CS.f = OS.f;
	      Naccept ++;
	   }
       else
       {
	      Nreject ++;
       }
       
	   if(Naccept > pho*Nreject)
	   {
	      E_diff *= alpha;
	        
	   }
	   else
	   {
		  E_diff /= alpha; 
	   }
	   
	 // printf("%d   %d    %e    %f\n",Naccept, Nreject, E_diff, K);
	   if(OS.f < BestS.f )
	   {
          for(i=0;i<n;i++) BestS.x[i]= OS.x[i];
          BestS.f = OS.f;
          BestS.CurrentR = CR;
          pp = 0.75*T;
          NoImprove = 0;
		 // printf("second stage improved ! f = %12.10e\n",CS.f);
	   }
	   else NoImprove ++;
	   
	   if(BestS.f < 1.0e-25) break;
	   step ++;
	   pp++;
     }
}

void Thredshold_Search(Solution &S, double CR, const int flag)
{
     int i, step, pp;
     INT MaxSteps  = 1000;
     INT MaxNoImprove = 400;
     INT NoImprove;
     double K, K1, K_min, K_min1, K_delt, K_delt1, sigma, T;
	 double pho = 0.5, alpha = 0.6; 
	 double E_diff = 1.0e-4;
	 long int Naccept = 0; // the number of acceptances
	 long int Nreject = 0; // the number of rejections
     double delt_f;

     R =  CR; 
     for(i=0;i<n;i++) CS.x[i] = S.x[i];
	 LocalSearch(CS);

     for(i=0;i<n;i++) BestS.x[i]= CS.x[i];
     BestS.f = CS.f;

	 T = 72.0; // T is a period
	 sigma = 2.0*pi/T;
	 K_min = 0.2;
	 K_delt= 0.125; 

     step = 0;
	 pp = 0;
	 NoImprove = 0 ;
     while(NoImprove < MaxNoImprove) // MaxNoImprove is the search depth 
     {

        for(i=0;i<n;i++) OS.x[i]= CS.x[i];
        OS.f = CS.f;

        if(flag ==1)
        {
        	K = K_min + K_delt*(1.0 + sin(pp*sigma));  // K_max = 0.45
            GeneralPertubationReduce(OS,50, K, 4*K, 0.99);  // Sequential Random Perturbation 
		}
		else 
		{
            Shake(OS, n, 0.8);  // K(=0.8) is the perturbation strength, Uniformly random perturbation
		}
	
       LocalSearch(OS);  
       delt_f = OS.f - CS.f ;
       if(delt_f < E_diff && fabs(delt_f) > 1.0e-20 )
	   {
		  for(i=0;i<n;i++) CS.x[i]= OS.x[i];
          CS.f = OS.f;
	      Naccept ++;
	   }
       else
       {
	      Nreject ++;
       }

      if(Naccept > pho*Nreject)
	   {
	       E_diff *= alpha;
	   }
	   
	   else
	   {
		  E_diff /= alpha;
	   }

	   if(OS.f < BestS.f)
	   {
          for(i=0;i<n;i++) BestS.x[i] = OS.x[i];
          BestS.f = OS.f;
          BestS.CurrentR = CR;
          pp = 0.75*T;
          NoImprove = 0 ;
		  // printf("first stage improved ! f = %12.10e\n",CS.f);
	   }
	   else NoImprove ++;

	   if(BestS.f < 1.0e-25) break;
	   step ++;
	   pp++;
     }
     
     for(i=0;i<n;i++) S.x[i] = BestS.x[i]; 
	 S.CurrentR = R; 
	 S.f = myvalue(S.x,n);  
}

void TA_ECP()
{
    int i;
    int Flag;  
	double ff;
    double CR;
	double SmallerR;
	double beta;
	double delt_beta = 0.001;
    starting_time = clock();
    
    /* The first stage aims to find an intial radius of container */
	beta = 0.7;  // Packing density of initial solution 
    CR = sqrt(pi*N/(4.0*beta) ) ;
    initial(NS, n, CR - 1.0);
    R =  CR; // R= 0.5*L
    LocalSearch(NS);
   	ff = myvalue(NS.x,n); 
    while(ff <= 1.0e-25) 
    {
		beta += (0.5 + 0.5*(rand()%RAND_MAX)/RAND_MAX)*delt_beta;
		CR = sqrt(pi*N/(4.0*beta) ) ; 
 	    Thredshold_Search(NS, CR, 1); // Thresholding Search Procedure
        R =  CR;
		ff = myvalue(NS.x,n);
	}
	
    for(i=0;i<n;i++) FS.x[i] = NS.x[i];
    FS.x[n] = CR;
    AdjustRadius(FS);  // Container Adjustment Procedure
    
    R_min = FS.x[n];  // R_min  = L_min / 2 
    for(i=0; i<n; i++) CBestS.x[i] = BestS.x[i]; 
    CBestS.CurrentR = FS.x[n];
    R = FS.x[n];
    CBestS.f = myvalue(CBestS.x,n);
    final_time = 1.0*(clock()- starting_time)/CLOCKS_PER_SEC ;

   // printf("intial L =%12.15f\n",2.0*R_min) ;
   
   /*  The second stage aims to find an improving packing configuration */
    Flag = 1;
    while(1.0*(clock()- starting_time)/CLOCKS_PER_SEC  < Time_limit)
	{
		SmallerR = R_min;
	    Thredshold_Accepting(SmallerR, Flag);  //Thresholding Search Procedure
		R = SmallerR;
		if(myvalue(BestS.x,n) < 1.0e-25) // if BestS is a feasible solution 
		{
             for(i=0;i<n;i++) FS.x[i] = BestS.x[i];
             FS.x[n] = SmallerR;
             AdjustRadius(FS);   // Container adjustment procedure
             if(FS.x[n] < R_min - 1.0e-11)
             {
                for(i=0;i<n;i++) CBestS.x[i] = FS.x[i];
                CBestS.CurrentR = FS.x[n]; 
                R= CBestS.CurrentR; 
                CBestS.f = myvalue(CBestS.x,n); 
                R_min =  CBestS.CurrentR; 
                final_time = 1.0*(clock()- starting_time)/CLOCKS_PER_SEC ; 
               // printf("Flag =%d  new R = %12.15f time =%f seconds\n",Flag, 2.0*R_min, final_time) ;
             }
             else
			 {
             	Flag = Flag  + 1 ; 
             	Flag = Flag  % 2 ; 
			 } // Change the type of perturbation operator 
	    }
	    else
		{	
             Flag = Flag  + 1 ; 
             Flag = Flag  % 2 ; 
		}  // Change the type of perturbation operator 
    }

}
