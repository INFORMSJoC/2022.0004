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
    double x_i, x_j, y_i, y_j, z_i, z_j;
    INT I, J, i, j;
    f = 0.0;
     // The sphere si is denoted by (x_i,y_i,z_i), and the sphere sj is denoted by (x_j,y_j,z_j), 2*R (=L) denotes the size of container.
   	for (I = 0; I < n/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

		overlap =  fabs(x_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(y_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(z_i)+ 1.0 - R ;
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
			SumR =  2.0;

			if(fabs(x_i - x_j) > SumR || fabs(y_i - y_j) > SumR || fabs(z_i - z_j) > SumR) continue;
			dist = sqrt((x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) + (z_i -z_j)*(z_i - z_j));
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
    double x_i, x_j, y_i, y_j, z_i, z_j;
    INT I, J, i, j;

	for (i = 0; i < n; i++) g [i] = 0;
	for (I = 0; I < n/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

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
        
       	if(z_i > 0)
		{
		  overlap  = z_i + 1.0 - R;
		  if(overlap > 0)
	   	  {
		     g[i+2] += 2.0*overlap;
          }
        }
		else if(z_i < 0)
	    {
		  overlap  = -1.0*z_i + 1.0 - R;
		  if(overlap > 0)
		  {
		    g[i+2] -= 2.0*overlap;
          }
        }
        
	}  // Calculate the force of each sphere s_i

	for (I = 0; I < n/3; I++)
	{
 	    i   = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

		for (J = I+1; J < n/3; J++)
		{
   	        j   = 3*J;
			x_j = x[j];
			y_j = x[j+1];
			z_j = x[j+2];

            SumR =  2.0;
			if(fabs(x_i - x_j) > SumR || fabs(y_i - y_j) > SumR || fabs(z_i - z_j) > SumR ) continue;

	 	    dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) + (z_i - z_j)*(z_i - z_j) );
	 	    overlap = SumR - dist;
	 	    if(overlap > 0 )
	 	    {
			   g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			   g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;
			   g[i+2] -=  2.0*overlap*(z_i - z_j) /dist;
			   g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	           g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
	           g[j+2] -=  2.0*overlap*(z_j - z_i) /dist;
			}
		}
	}  // Calculate the forces between the unit spheres
    return ;
}

double myvalueNN
(
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j, z_i, z_j;
    INT I, J, i, j;
    f = 0.0 ;
// The sphere si is denoted by (x_i,y_i,z_i), and the sphere sj is denoted by (x_j,y_j,z_j), 2*R (=L) denotes the size of container.
    for (I = 0; I < n/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];
		
		overlap =  fabs(x_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(y_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(z_i)+ 1.0 - R ;
        if(overlap > 0.0) f += (overlap*overlap);
	}

	for (I = 0; I < n/3; I++)
	{
		i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

	//	printf("\n number of neighbors = %d\n",NNeighbors[I]) ;
		for (J = 0; J < NNeighbors[I]; J++)
		{

			j = 3*Adjacent[I][J];  // printf("%d ", Adjacent[I][J]);
			x_j = x[j];
			y_j = x[j+1];
			z_j = x[j+2];

			dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) + (z_i - z_j)*(z_i -z_j) );
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
    double x_i, x_j, y_i, y_j, z_i, z_j;
    INT I, J, i, j;

	for (i = 0; i < n; i++) g [i] = 0;
    for (I = 0; I < n/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

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
        
       	if(z_i > 0)
		{
		  overlap  = z_i + 1.0 - R;
		  if(overlap > 0)
	   	  {
		     g[i+2] += 2.0*overlap;
          }
        }
		else if(z_i < 0)
	    {
		  overlap  = -1.0*z_i + 1.0 - R;
		  if(overlap > 0)
		  {
		    g[i+2] -= 2.0*overlap;
          }
        }
	}  // Calculate the force of each sphere s_i

	for (I = 0; I < n/3; I++)
	{
 	    i   = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

 	    for(J = 0; J < NNeighbors[I]; J++)
		{
		   	j = 3*Adjacent[I][J];
			x_j = x[j];
			y_j = x[j+1];
			z_j = x[j+2];

	 	    dist = sqrt((x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j) + (z_i - z_j)*(z_i - z_j) );
	 	    overlap = 2.0 - dist;
	 	    if(overlap > 0 )
	 	    {

			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;
	         g[i+2] -=  2.0*overlap*(z_i - z_j) /dist;

			 g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	         g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
	         g[j+2] -=  2.0*overlap*(z_j - z_i) /dist;
			}
		}
	} //Calculate the forces between the spheres
	
    return ;
}

double myvalueRR
(
    double   *x,
    INT       n
)
{
    double f, dist, overlap;
    double x_i, x_j, y_i, y_j, z_i, z_j;
    INT I, J, i, j;
    f = 0.0 ;
    // The sphere si is denoted by (x_i,y_i,z_i), and the sphere sj is denoted by (x_j,y_j,z_j), 2*R (=L) denotes the size of container.
    for(I = 0; I < (n-1)/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];
		
		overlap =  fabs(x_i)+ 1.0 - x[n-1]  ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(y_i)+ 1.0 - x[n-1] ;
        if(overlap > 0.0) f += (overlap*overlap);
		overlap =  fabs(z_i)+ 1.0 - x[n-1] ;
        if(overlap > 0.0) f += (overlap*overlap);
        
	}

	for (I = 0; I < (n-1)/3; I++)
	{
		i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

	//	printf("\n number of neighbors = %d\n",NNeighbors[I]) ;
		for (J = 0; J < NNeighbors[I]; J++)
		{

			j = 3*Adjacent[I][J];  // printf("%d ", Adjacent[I][J]);
			x_j = x[j];
			y_j = x[j+1];
			z_j = x[j+2];

			dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) + (z_i - z_j)*(z_i - z_j) );
			overlap = 2.0 - dist ;
			if(overlap > 0.0) f += (overlap*overlap);
		}

	}

   f  += x[n-1]*x[n-1]/alpha;

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
    double x_i, x_j, y_i, y_j, z_i, z_j;
    INT I, J, i, j;

	for (i = 0; i < n; i++) g [i] = 0;
    for (I = 0; I < (n-1)/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

  	    if(x_i >= 0)
		{
		  overlap  = x_i + 1.0 - x[n-1];
		  if(overlap > 0)
		  {
		     g[i]   +=  2.0*overlap;
          }
	    }
		else
		{
		   overlap  = -1.0*x_i + 1.0 - x[n-1];
    	   if(overlap > 0)
		   {
		   g[i] -=  2.0*overlap;
           }
        }
		if(y_i >= 0)
		{
		  overlap  = y_i + 1.0 - x[n-1];
		  if(overlap > 0)
	   	  {
		   g[i+1]   +=  2.0*overlap;
          }
        }
		else
	    {
		  overlap  = -1.0*y_i + 1.0 - x[n-1];
		  if(overlap > 0)
		  {
		   g[i+1]   -=  2.0*overlap;
          }
        }
        
       	if(z_i >= 0)
		{
		  overlap  = z_i + 1.0 - x[n-1];
		  if(overlap > 0)
	   	  {
		   g[i+2]   +=  2.0*overlap;
          }
        }
		else
	    {
		  overlap  = -1.0*z_i + 1.0 - x[n-1];
		  if(overlap > 0)
		  {
		   g[i+2]   -=  2.0*overlap;
          }
        }
        
	}  // Calculate the force of each sphere s_i
	
    for (I = 0; I < (n-1)/3; I++)
	{
	 	i = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

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
          
        overlap  = fabs(z_i) + 1.0 - x[n-1];
        if(overlap > 0)
	   	  {
		     g[n-1]   -=  2.0*overlap;
          }

	}   // Calculate the force of each sphere s_i

	for (I = 0; I < (n-1)/3; I++)
	{
 	    i   = 3*I;
		x_i = x[i];
		y_i = x[i+1];
		z_i = x[i+2];

 	    for(J = 0; J < NNeighbors[I]; J++)
		{
		   	j = 3*Adjacent[I][J];
			x_j = x[j];
			y_j = x[j+1];
			z_j = x[j+2];

	
	 	    dist = sqrt((x_i - x_j)*(x_i - x_j)+(y_i - y_j)*(y_i - y_j) + (z_i - z_j)*(z_i - z_j) );
	 	    overlap = 2.0 - dist;
	 	    if(overlap > 0 )
	 	    {

			 g[i]   -=  2.0*overlap*(x_i - x_j) /dist;
			 g[i+1] -=  2.0*overlap*(y_i - y_j) /dist;
		     g[i+2] -=  2.0*overlap*(z_i - z_j) /dist;

			 g[j]   -=  2.0*overlap*(x_j - x_i) /dist;
	         g[j+1] -=  2.0*overlap*(y_j - y_i) /dist;
             g[j+2] -=  2.0*overlap*(z_j - z_i) /dist;
	         
			}
			
		}
		
	}
	
	g[n-1] += 2*x[n-1]/alpha;

}

