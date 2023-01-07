#include <math.h>
#include "cg_user.h"
#include "Common.h"

double myvalue
(
    double   *x,
    INT       n
)
{
     double f, dist, overlap, SumR;
    double x_i, x_j, y_i, y_j;
    INT I, J, i, j;
    f = 0.0;
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), 2*R denotes the size of container.
   	for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		overlap =  fabs(x_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(y_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
	}

	for (I = 0; I < n/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		for (J = I+1; J < n/2; J++)
		{
			j = 2*J;
			x_j = x[j];
			y_j = x[j+1];
			SumR =  2.0;

			if(fabs(x_i - x_j) > SumR || fabs(y_i - y_j) > SumR ) continue;

			dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = SumR - dist ;
			if(overlap > 0.0) f += (overlap*overlap);
		}

	}

    return (f) ;
}

void mygrad
(
    double    *g,
    double    *x,
    INT        n
)
{
    double dist, overlap, SumR ;
    double x_i, x_j, y_i, y_j;
    INT I, J, i, j;

	for (i = 0; i < n; i++) g [i] = 0;
	for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		if(x_i > 0)
		{
		  overlap  = x_i + 1.0 - R;
		  if(overlap > 0)
		  {
		     g[i] += 2.0*overlap;
          }
	    }
		else if(x_i < 0)
		{
		   overlap  = -1.0*x_i + 1.0 - R;
    	   if(overlap > 0)
		   {
		     g[i] -= 2.0*overlap;
           }
        }
		if(y_i > 0)
		{
		  overlap  = y_i + 1.0 - R;
		  if(overlap > 0)
	   	  {
		     g[i+1] += 2.0*overlap;
          }
        }
		else if(y_i < 0)
	    {
		  overlap  = -1.0*y_i + 1.0 - R;
		  if(overlap > 0)
		  {
		    g[i+1] -= 2.0*overlap;
          }
        }
	}  // Calculate the force of each circle c_i

	for (I = 0; I < n/2; I++)
	{
 	    i   = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		for (J = I+1; J < n/2; J++)
		{
   	        j   = 2*J;
			x_j = x[j];
			y_j = x[j+1];

            SumR =  2.0;
			if(fabs(x_i - x_j) > SumR || fabs(y_i - y_j) > SumR ) continue;

	 	    dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
	 	    overlap = SumR - dist;
	 	    if(overlap > 0 )
	 	    {
			   g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			   g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;
			   g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	           g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
			}
		}
	}  // Calculate the forces between the circles
    return ;
}

double myvalueNN
(
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    INT I, J, i, j;
    f = 0.0 ;
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), 2*R denotes the size of container.
    for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		overlap =  fabs(x_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(y_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
	}

	for (I = 0; I < n/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

	//	printf("\n number of neighbors = %d\n",NNeighbors[I]) ;
		for (J = 0; J < NNeighbors[I]; J++)
		{

			j = 2*Adjacent[I][J];  // printf("%d ", Adjacent[I][J]);
			x_j = x[j];
			y_j = x[j+1];

          //  if((x_i - x_j)*(x_i - x_j) > 4.0 || (y_i - y_j)*(y_i - y_j) > 4.0 ) continue;
			dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 - dist ;
			if(overlap > 0.0) f += (overlap*overlap);
		}

	}
    return (f) ;
}

void mygradNN
(
    double    *g,
    double    *x,
    INT        n
)
{
    double dist, overlap ;
    double x_i, x_j, y_i, y_j;
    INT I, J, i, j;

	for (i = 0; i < n; i++) g [i] = 0;
    for (I = 0; I < n/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		if(x_i > 0)
		{
		  overlap  = x_i + 1.0 - R;
		  if(overlap > 0)
		  {
		     g[i] += 2.0*overlap;
          }
	    }
		else if(x_i < 0)
		{
		   overlap  = -1.0*x_i + 1.0 - R;
    	   if(overlap > 0)
		   {
		     g[i] -= 2.0*overlap;
           }
        }
		if(y_i > 0)
		{
		  overlap  = y_i + 1.0 - R;
		  if(overlap > 0)
	   	  {
		     g[i+1] += 2.0*overlap;
          }
        }
		else if(y_i < 0)
	    {
		  overlap  = -1.0*y_i + 1.0 - R;
		  if(overlap > 0)
		  {
		    g[i+1] -= 2.0*overlap;
          }
        }
	}  // Calculate the force of each circle c_i

	for (I = 0; I < n/2; I++)
	{
 	    i   = 2*I;
		x_i = x[i];
		y_i = x[i+1];

 	    for(J = 0; J < NNeighbors[I]; J++)
		{
		   	j = 2*Adjacent[I][J];
			x_j = x[j];
			y_j = x[j+1];

		   // if((x_i - x_j)*(x_i - x_j) > 4.0 || (y_i - y_j)*(y_i - y_j) > 4.0 ) continue;
	 	    dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
	 	    overlap = 2.0 - dist;
	 	    if(overlap > 0 )
	 	    {

			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;

			 g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	         g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
			}
		}
	} //Calculate the forces between the circles
	
    return ;
}

double myvalueRR
(
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j;
    INT I, J, i, j;
    f = 0.0 ;
    // The circle ci is denoted by (x_i,y_i), and the circle cj is denoted by (x_j,y_j), 2*R denotes the size of container.
    for(I = 0; I < (n-1)/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

		overlap =  fabs(x_i)+ 1.0 - x[n-1]  ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(y_i)+ 1.0 - x[n-1] ;
        if(overlap > 0.0) f += (overlap*overlap);
	}

	for (I = 0; I < (n-1)/2; I++)
	{
		i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

	//	printf("\n number of neighbors = %d\n",NNeighbors[I]) ;
		for (J = 0; J < NNeighbors[I]; J++)
		{

			j = 2*Adjacent[I][J];  // printf("%d ", Adjacent[I][J]);
			x_j = x[j];
			y_j = x[j+1];

          //  if((x_i - x_j)*(x_i - x_j) > 4.0 || (y_i - y_j)*(y_i - y_j) > 4.0 ) continue;
			dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
			overlap = 2.0 - dist ;
			if(overlap > 0.0) f += (overlap*overlap);
		}

	}
	
    f += x[n-1]*x[n-1]/alpha;
    
    return (f) ;
}

void mygradRR
(
    double    *g,
    double    *x,
    INT        n
)
{
    double dist, overlap ;
    double x_i, x_j, y_i, y_j;
    INT I, J, i, j;

	for (i = 0; i < n; i++) g [i] = 0;
    for (I = 0; I < (n-1)/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

  	    if(x_i > 0)
		{
		  overlap  = x_i + 1.0 - x[n-1];
		  if(overlap > 0)
		  {
		     g[i]   +=  2.0*overlap;
          }
	    }
		else if(x_i < 0)
		{
		   overlap  = -1.0*x_i + 1.0 - x[n-1];
    	   if(overlap > 0)
		   {
		   g[i] -=  2.0*overlap;
           }
        }
		if(y_i > 0)
		{
		  overlap  = y_i + 1.0 - x[n-1];
		  if(overlap > 0)
	   	  {
		   g[i+1]   +=  2.0*overlap;
          }
        }
		else if(y_i < 0)
	    {
		  overlap  = -1.0*y_i + 1.0 - x[n-1];
		  if(overlap > 0)
		  {
		   g[i+1]   -=  2.0*overlap;
          }
        }
	}  // Calculate the force of each circle c_i
	
	
    for (I = 0; I < (n-1)/2; I++)
	{
	 	i = 2*I;
		x_i = x[i];
		y_i = x[i+1];

        overlap  = fabs(x_i) + 1.0 - x[n-1];
        if(overlap > 0)
		  {
		     g[n-1]   -=  2.0*overlap;
          }

        overlap  = fabs(y_i) + 1.0 - x[n-1];
        if(overlap > 0)
	   	  {
		     g[n-1]   -=  2.0*overlap;
          }

	}   // Calculate the force of each circle c_i

	for (I = 0; I < (n-1)/2; I++)
	{
 	    i   = 2*I;
		x_i = x[i];
		y_i = x[i+1];

 	    for(J = 0; J < NNeighbors[I]; J++)
		{
		   	j = 2*Adjacent[I][J];
			x_j = x[j];
			y_j = x[j+1];

		   // if((x_i - x_j)*(x_i - x_j) > 4.0 || (y_i - y_j)*(y_i - y_j) > 4.0 ) continue;
	 	    dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) );
	 	    overlap = 2.0 - dist;
	 	    if(overlap > 0 )
	 	    {

			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;

			 g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	         g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
	         
			}
			
		}
		
	}
	
	
	g[n-1] += 2*x[n-1]/alpha;
	
    return ;
}

