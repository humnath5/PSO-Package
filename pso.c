#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "matmul2.h"
#include "randomlib.h"
#include "pso.h"
#include "assert.h"
#define SQRT_EPSILON  0.000001
#define CALLOC calloc
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/* defining global variables */
fptr FUN;
double* GBEST;
int ND;
int INDEX;


void pso_optimization(fptr fun, int np, int nd, int ni,double* lb, double* ub, double* value, double* gbest, int ni_stop, double tol){

  int i,j,k;

  double old_value;
  double vb;
 /*double w = 0.7968, c1 = 1.4962 , c2 = 1.4962;*/
 double w = RandomDouble(0.5,1.0) , c1 = RandomDouble(1.4,2.0), c2 = RandomDouble(1.4,2.0);
/*
	xpos = an (np x nd) matrix of positions of particles over dimensions
	rows = particles
	columns = dimensions
	xvel = velocity of particles over dimensions. Similar structure to xpos (np x nd)
	lbest = each particles best achieved position, or its local best. (np x nd)
	gbest = best positi0on found by all particles, or the best lbest position (nd x 1)
*/
	struct matrix * xpos, * xvel, * lbest;
	xpos  = get_matrix(np,nd);
	xvel  = get_matrix(np,nd);
	lbest = get_matrix(np,nd);

/*AssINDEXn xpos and xvel a random number to begin optimization*/
	double * ltemp = (double*)malloc(sizeof(double)*np);
	double * xtemp = (double*)malloc(sizeof(double)*np);

    int np_loc = int(np/4); /* number of particles around global best position*/

    double eps , u,l ;

     double alpha = 0.1;

   // printf("nploc %d\n",np_loc);
	for(i=0; i< np ; i++) {

	   ltemp[i] = 10000000;

	   if( i < np_loc) {

		for(j= 0 ;j< nd ; j++){

		      eps = alpha*(ub[j]- lb[j]);
              l  =  MAX(gbest[j]- eps, lb[j]);
              u =   MIN(gbest[j] + eps, ub[j]);
              vb =  (u-l)/10;

              (xpos->mat)[i][j]  = RandomDouble(l,u);
			  (xvel->mat)[i][j]  = RandomDouble(-vb,vb);
			  (lbest->mat)[i][j] = (xpos->mat)[i][j];
			  }
            }

        else {
            for(j = 0; j< nd ; j++){
             vb = (ub[j]- lb[j])/10;
			(xpos->mat)[i][j] = RandomDouble(lb[j], ub[j]);
			(xvel->mat)[i][j] = RandomDouble(-vb,vb);
			(lbest->mat)[i][j] = (xpos->mat)[i][j];
        }
	  }
   }

    double best  = fun(gbest,nd);

    old_value = best;

	for(j=0; j< nd ;j++){  xpos->mat[0][j] = gbest[j];}

	/*Creates names of files to save intermediate results when str != NULL*/
	double * r1 = (double*)malloc(sizeof(double)* np);
	double * r2 = (double*)malloc(sizeof(double)* np);

	/*For each iteration*/
	for( k= 0; k < ni; k++){

	   /*for each particle*/
		for(i=0; i< np; i++){
			r1[i] = RandomDouble(0,1);
			r2[i] = RandomDouble(0,1);
			}
		for( i=0; i< np ; i++){
		/* Update Velocity */
			double vb;
			for(j= 0;j < nd; j++)
			{
                vb = (ub[j]- lb[j])/10;

                  (xvel->mat)[i][j] = (1*w*(xvel->mat)[i][j]) + (1*c1*r1[i]*((lbest->mat)[i][j] - (xpos->mat)[i][j])) + (1*c2*r2[i]*(gbest[j] - (xpos->mat)[i][j]));

				if((xvel->mat)[i][j] >  vb) { (xvel->mat)[i][j] =  vb; }
				if((xvel->mat)[i][j] < -vb) { (xvel->mat)[i][j] = -vb ; }

				(xpos->mat)[i][j] = (xpos->mat)[i][j] + (xvel->mat)[i][j];

				if((xpos->mat)[i][j] >  ub[j]) { (xpos->mat)[i][j] = ub[j]; }
				if((xpos->mat)[i][j] <  lb[j]) { (xpos->mat)[i][j] = lb[j]; }
			}
           	xtemp[i] = fun((xpos->mat)[i],nd);
		}
    /*First find each best local position and function value,
    then find the best global position and function value.*/
		for( i=0;i< np;i++)
		{
			if(xtemp[i] < ltemp[i])
			{
				update_best((lbest->mat)[i],(xpos->mat)[i],nd);
				ltemp[i] = xtemp[i];
				if(xtemp[i] < best)
				{
					update_best(gbest, (xpos->mat)[i],nd);
					 best = xtemp[i];
				}
			}
		}
         // printf("run %d\t old %f\t new %f\n",k,old_value,pso->value);

          if( (k > 0 ) && (k% ni_stop ==0)){

                        double estimate = fabs((best - old_value)/( best + tol));

                        if(estimate < tol){
                          // printf(" No singnificant improvement!!\t Stopping after %d iteration \n",k);
                          break;
                        }}
         old_value = best;

	}

    best  =  fun(gbest, nd);


    (*value) = best;

	free(r1); free(r2); free(ltemp); free(xtemp);
	free_matrix(xpos); free_matrix(xvel); free_matrix(lbest);

} /* end pso() */



void multi_pso_optimization(fptr fun, int np , int nd, int ni, int maxRun, int* type, double* lb, double* ub, double* value ,double* gbest, int ni_stop, double tol)
{   int i;
    double alpha = 1;
    for( i = 0 ; i < maxRun ; i++) {
       pso_optimization(fun, np , nd, ni, lb, ub, value, gbest, ni_stop,tol);
      // printf("Value after %d run %lf\n ", i, (*value));
       adjust_domain(type, lb , ub, gbest,nd,alpha);
       alpha = alpha/1.1;

    }

} /* end pso() */


void multi_pso_local_optimizer(fptr fun, int np, int nd, int ni, int maxRun, int* type, double* lb, double* ub, double* value, double* gbest, int ni_stop, double tol)
{    int i;

     FUN = fun;
     GBEST = gbest;
     ND = nd;

     double alpha = 1;
     for( i = 0 ; i < maxRun ; i++) {

       pso_optimization(FUN, np , ND, ni, lb, ub, value, GBEST, ni_stop,tol);
       //printf("Value after %d run  %f\n", i, (*value));
       local_optimizer(FUN, ND, lb, ub, value, GBEST, ni_stop, tol);

       //printf("Value after %d run  %f\n", i, (*value));

       adjust_domain(type, lb , ub, GBEST, ND, alpha);
       alpha = alpha/1.1;
    }


} /* end pso() */




void adjust_domain(int*type, double* lb , double* ub, double* gbest, int nd, double alpha){

	int i;
    double epsilon = 0.0000001;
	for ( i= 0; i < nd ; i++) {

            if( type[i] == 0){ lb[i] = lb[i]; ub[i] = ub[i];}/* bounds are unchanged*/

            if( type[i] == -1) {
                  if(lower_bound_check(gbest[i],lb[i],ub[i], epsilon) == 1) { lb[i] =  lb[i] - alpha*(ub[i] - lb[i]);}
                  }/* inner and outer if both close */

		    if( type[i] == 1){
                  if(upper_bound_check(gbest[i], lb[i],ub[i],epsilon) == 1) { ub[i] =  ub[i] + alpha*(ub[i] - lb[i]);}
                  }/* inner and outer if both close */

		} /* inner and outer if both close */

}


void update_best(double *gbest , double * x,int nd){
	int i ;
	double temp;
	for(i=0;i< nd ;i++)
        {
            temp = x[i];
            gbest[i] = temp;
          }
    }



int lower_bound_check(double g, double  l, double u, double epsilon)
  {
         int n = 0;
         if( ( u- l ) > 10*epsilon)
               {
                  if((( g -l) < epsilon))
                   n = 1;
                 }
          return n;
} /* end for */


int upper_bound_check(double g, double  l, double u, double epsilon)
  {
         int n = 0;
         if( ( u - l ) > 10*epsilon)
               {
                  if((( u - g) < epsilon))
                   n = 1;
                 }
          return n;
} /* end for */

void initialize_pso_domain(char * bounds_file, int* type, double* lb, double* ub, double* gbest, int nd){

	/*::::::Choosing Starting  Domain size::::*/
	int i;
	FILE* fp;
	fp = fopen(bounds_file,"r");
	for ( i =0; i < nd; i ++){
         fscanf(fp,"%d %lf  %lf ", &type[i], &lb[i], &ub[i]);

	}
	fclose(fp);

	for (i = 0; i < nd ; i++) { gbest[i] = (ub[i] - lb[i]) / 3; }

}



void  local_optimizer(fptr fun,int nd, double* lb, double* ub, double* value, double* gbest, int ni_stop, double tol) {
		 /* this is extern global variable*/
		int  k;
		double l, u;
		double temp_t, old_value;
		double eps;

         int ntries = 20;
		double best = fun(gbest, nd);
		double alpha = 0.1;

		old_value = best;
		double val;
		//printf("\nbefore local optimizer %f\n",val_old);
		for (k = 0; k< ntries; k++) {

			for (INDEX  = 0; INDEX  < nd; INDEX ++) {

				eps = alpha*(ub[INDEX ] - lb[INDEX ]);

				l = MAX(gbest[INDEX ] - eps, lb[INDEX ]);

				u = MIN(gbest[INDEX ] + eps, ub[INDEX ]);

				temp_t = gbest[INDEX ];

				gbest[INDEX ] = brent_algorithm(l, u, fun_1D, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[INDEX ] = temp_t; }
				else { best = val; }


			}

			for ( INDEX= nd - 1; INDEX >= 0; INDEX--) {


				eps = alpha*(ub[INDEX ] - lb[INDEX ]);

				l = MAX(gbest[INDEX ] - eps, lb[INDEX ]);

				u = MIN(gbest[INDEX ] + eps, ub[INDEX ]);

				temp_t = gbest[INDEX ];

				gbest[INDEX ] = brent_algorithm(l, u,fun_1D, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[INDEX ] = temp_t; }
				else { best = val; }


			}

			for (INDEX= nd - 1; INDEX >= 0; INDEX--) {


				eps = alpha*(ub[INDEX ] - lb[INDEX ]);

				l = MAX(gbest[INDEX ] - eps, lb[INDEX ]);

				u = MIN(gbest[INDEX ] + eps, ub[INDEX ]);

				temp_t = gbest[INDEX ];

				gbest[INDEX ] = brent_algorithm(l, u,fun_1D, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[INDEX ] = temp_t; }
				else { best = val; }

			}

			for (INDEX = 0; INDEX < nd; INDEX++) {


				eps = alpha*(ub[INDEX ] - lb[INDEX ]);

				l = MAX(gbest[INDEX ] - eps, lb[INDEX ]);

				u = MIN(gbest[INDEX ] + eps, ub[INDEX ]);

				temp_t = gbest[INDEX ];

				gbest[INDEX ] = brent_algorithm(l, u, fun_1D, tol);

				val = fun(gbest, nd);

				if (val > best) { gbest[INDEX ] = temp_t; }
				else { best = val; }

			}

			best = fun(gbest, nd);
			// printf("\n val_old  %f\t  val_best%f\n",val_old,val_best);
			if ((k > 0) && (k % ni_stop == 0)) {
				double estimate = fabs((best - old_value) / (best + tol));
				if (estimate < tol) break;
			}
			old_value = best;
		}
		best = fun(gbest, nd);

		(*value) = best;
	}



/*
************************************************************************
*	    		    C math library
* function brent_algorithm - one-dimensional search for a function minimum
*			  over the given range
*
* Input
*	double brent_algorithm(a,b,f,tol)
*	double a; 			Minimum will be seeked for over
*	double b;  			a range [a,b], a being < b.
*	double (*f)(double x);		Name of the function whose minimum
*					will be seeked for
*	double tol;			Acceptable tolerance for the minimum
*					location. It have to be positive
*					(e.g. may be specified as EPSILON)
*
* Output
*	brent_algorithm returns an estimate for the minimum location with accuracy
*	3*SQRT_EPSILON*abs(x) + tol.
*	The function always obtains a local minimum which coincides with
*	the global one only if a function under investINDEXation being
*	unimodular.
*	If a function being examined possesses no local minimum within
*	the given range, brent_algorithm returns 'a' (if f(a) < f(b)), otherwise
*	it returns the rINDEXht range boundary value b.
*
* Algorithm
*	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
*	computations. M., Mir, 1980, p.202 of the Russian edition
*
*	The function makes use of the "gold section" procedure combined with
*	the parabolic interpolation.
*	At every step program operates three abscissae - x,v, and w.
*	x - the last and the best approximation to the minimum location,
*	    i.e. f(x) <= f(a) or/and f(x) <= f(b)
*	    (if the function f has a local minimum in (a,b), then the both
*	    conditions are fulfiled after one or two steps).
*	v,w are previous approximations to the minimum location. They may
*	coincide with a, b, or x (although the algorithm tries to make all
*	u, v, and w distinct). Points x, v, and w are used to construct
*	interpolating parabola whose minimum will be treated as a new
*	approximation to the minimum location if the former falls within
*	[a,b] and reduces the range enveloping minimum more efficient than
*	the gold section procedure.
*	When f(x) has a second derivative positive at the minimum location
*	(not coinciding with a or b) the procedure converges superlinearly
*	at a rate order about 1.324
*
************************************************************************


#include "assert.h"
#include "math.h"
#include "stdio.h"


#define SQRT_EPSILON  .000001



/* An estimate to the min location*/

/*
double a;				/* Left border | of the range
double b;  				/* RINDEXht border| the min is seeked
double (*f)(double x);			/* Function under investINDEXation
double tol;				/* Acceptable tolerance
*/
double brent_algorithm(double a, double b, double(*f)(double x), double tol)
{
	double x, v, w;				/* Abscissae, descr. see above	*/
	double fx;				/* f(x)				*/
	double fv;				/* f(v)				*/
	double fw;				/* f(w)				*/
	const double r = (3. - sqrt(5.0)) / 2;	/* Gold section ratio		*/

	assert(tol > 0 && b > a);

	v = a + r*(b - a);  fv = (*f)(v);       /* First step - always gold section*/
	x = v;  w = v;
	fx = fv;  fw = fv;

	for (;;)		/* Main iteration loop	*/
	{
		double range = b - a;			/* Range over which the minimum */
										/* is seeked for		*/
		double middle_range = (a + b) / 2;
		double tol_act =			/* Actual tolerance		*/
			SQRT_EPSILON*fabs(x) + tol / 3;
		double new_step;      		/* Step at this iteration       */



		if (fabs(x - middle_range) + range / 2 <= 2 * tol_act)
			return x;				/* Acceptable approx. is found	*/

									/* Obtain the gold section step	*/
		new_step = r * (x<middle_range ? b - x : a - x);


		/* Decide if the interpolation can be tried	*/
		if (fabs(x - w) >= tol_act)		/* If x and w are distinct      */
		{					/* interpolatiom may be tried	*/
			register double p; 		/* Interpolation step is calcula-*/
			register double q;              /* ted as p/q; division operation*/
											/* is delayed until last moment	*/
			register double t;

			t = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v)*q - (x - w)*t;
			q = 2 * (q - t);

			if (q>(double)0)		/* q was calculated with the op-*/
				p = -p;			/* posite sINDEXn; make q positive	*/
			else				/* and assINDEXn possible minus to	*/
				q = -q;			/* p				*/

			if (fabs(p) < fabs(new_step*q) &&	/* If x+p/q falls in [a,b]*/
				p > q*(a - x + 2 * tol_act) &&		/* not too close to a and */
				p < q*(b - x - 2 * tol_act))            /* b, and isn't too large */
				new_step = p / q;			/* it is accepted         */
											/* If p/q is too large then the	*/
											/* gold section procedure can 	*/
											/* reduce [a,b] range to more	*/
											/* extent			*/
		}

		if (fabs(new_step) < tol_act)	/* Adjust the step to be not less*/
			if (new_step > (double)0)	/* than tolerance		*/
				new_step = tol_act;
			else
				new_step = -tol_act;

		/* Obtain the next approximation to min	*/
		{				/* and reduce the enveloping range	*/
			register double t = x + new_step;	/* Tentative point for the min	*/
			register double ft = (*f)(t);
			if (ft <= fx)
			{                                 /* t is a better approximation	*/
				if (t < x)			/* Reduce the range so that	*/
					b = x;                        /* t would fall within it	*/
				else
					a = x;

				v = w;  w = x;  x = t;		/* AssINDEXn the best approx to x	*/
				fv = fw;  fw = fx;  fx = ft;
			}
			else                              /* x remains the better approx  */
			{
				if (t < x)			/* Reduce the range enclosing x	*/
					a = t;
				else
					b = t;

				if (ft <= fw || w == x)
				{
					v = w;  w = t;
					fv = fw;  fw = ft;
				}
				else if (ft <= fv || v == x || v == w)
				{
					v = t;
					fv = ft;
				}
			}

		}			/* ----- end-of-block ----- */
	}		/* ===== End of loop ===== */

}



/* The following  one dimensional objective function is used by brent_algorithm() and CCD local optimizer*/
double fun_1D(double t){
	double sum;
	GBEST[INDEX] = t;
	sum = FUN(GBEST, ND);
}

