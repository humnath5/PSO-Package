#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "randomlib.h"
#include "matmul2.h"
#include "pso.h"

double Sphere(double * x,int nd);
double Rastrigin(double * x, int nd);
double Rosenbrock(double*x);
double Griewank(double * x, int nd);
double Ackley(double * x, int nd);
double Levy(double * x, int nd);
double DropWave(double * x, int nd);
double Testfun(double* x, int nd);
double Schwef(double*x, int nd);



int main(int argc,char*argv[]){
int i,j,k;
int seed1,seed2;
srand((unsigned)time(NULL));
seed1= rand()%30000;
seed2= rand()%30000;

/* Reading command line inputs */
char func[80];
strcpy(func, argv[1]);
fptr f  = functionLookup(func);

int np = atoi(argv[2]);
int nd = atoi(argv[3]);
int ni = atoi(argv[4]);
int maxRun =atoi(argv[5]);
double l  = atof(argv[6]);
double u =  atof(argv[7]);
int option =atoi(argv[8]);


/* Alocationg memory spaces */
int* type = (int*)malloc(sizeof(int)*nd);
double*lb = (double*)malloc(sizeof(double)*nd);
double*ub = (double*)malloc(sizeof(double)*nd);
double * gbest = (double*)malloc(sizeof(double)*nd);
for(i = 0;i < nd ;i++){
    type[i] = 0;
    lb[i] = l; ub[i] = u;
    gbest[i] = (ub[i]-lb[i])/3;
}

RandomInitialise(seed1, seed2);
double tol = 0.000001;
double value = 9999999;
int ni_stop = 100;
/*Start Optimization */
double start= time(NULL);
if(option==0){
pso_optimization(f,np,nd,ni,lb,ub,&value,gbest,ni_stop,tol);}
else if(option==1){
multi_pso_optimization(f,np,nd,ni,maxRun,type,lb,ub,&value,gbest,ni_stop,tol);}
else if(option==2){
multi_pso_local_optimizer(f,np,nd,ni,maxRun,type,lb,ub,&value,gbest,ni_stop,tol);}
else{printf("\n No option is found !!!\n");}
double finish= time(NULL);
double time=difftime(finish, start);

	FILE* fp;
	fp = fopen("outputpara.txt", "w");
	fprintf(fp, "\n PSO parameters used:\n");
	fprintf(fp, "seed1:%d\n", seed1);
	fprintf(fp, "seed2:%d\n", seed2);
	fprintf(fp, "nd:%d\n", nd);
	fprintf(fp, "np:%d\n", np);
	fprintf(fp, "ni:%d\n", ni);
	fprintf(fp, "maxRun:%d\n", maxRun);
	fprintf(fp, "Elapsed time:%lf\n", time);
	fprintf(fp, "best value: %lf\n", value);
	fprintf(fp, "\n Optimal parameters(gbest)\n");
	for (j = 0; j < nd; j++) { fprintf(fp, "%lf\n", gbest[j]); }
	fclose(fp);

free(lb);free(ub);free(gbest);free(type);

return 0;
}



double Sphere(double * x,int nd){
	double sum = 0;
	int i;
	for(i = 0; i < nd; i++){
	sum += pow((x)[i],2);
	}
	return(sum);
}

double Rosenbrock(double * x, int nd){
	double sum = 0;
	int i;
	for(i=0;i<(nd-1);i++){
		sum += 100*pow(x[i+1]- pow(x[i],2),2)+pow(x[i]-1,2);
	}
	return(sum);
}

double Rastrigin(double * x, int nd)
{

	double sum = 0, pi = 3.14159265359, temp;
	int i;
	for(i=0;i<nd;i++)
	{
		temp = (double)i/5.0;
		sum += pow((x)[i] - temp,2) - (10*cos(2*pi*((x)[i] - temp)));
	}
	return(10*nd + sum);
}

double Griewank(double * x, int nd){
	double sum = 0, prod = 1;
	int i;
	for(i=0;i<nd;i++)
	{
		sum += pow((x)[i]-1,2);
		prod *= cos(((x)[i]-1)/sqrt(i+1));
	}
	return(1 + sum/4000-prod);
}

double Ackley(double * x, int nd){
   double sum1 = 0, sum2 = 0, pi = 3.14159265359, temp;
   int i;
   for ( i = 0; i < nd; i++){
      temp = (double)i/5.0;
      sum1+= pow((x)[i]- temp,2);
      sum2+= cos(2*pi*((x)[i] - temp));
   }

   return -20*exp(-0.2*sqrt(sum1/nd)) - exp(sum2/nd) + 20 + exp(1);

}

double Levy(double * x, int nd)
{
	double sum = 0,pi = 3.14159265359;
	int i;
	for(i=0;i<(nd-1);i++)
	{
		sum += pow((x)[i]-1,2)*(1 + 10*pow(sin(pi * (x)[i] + 1),2));
	}
	double a = pow(sin(pi*(x)[0]),2);
	double b = (pow((x)[nd-1]-1,2)*(1 + pow(sin(2*pi*(x)[nd-1]),2)));
	return( a + b + sum);
}


double DropWave(double * x, int nd)
{

	double sum = 0.0;
	int i;

	for(i=0;i< nd ;i++)
	{
		sum += pow((x)[i],2);
	}

	return -(1 + cos(12*sqrt(sum)))/(2 + (0.5*(sum)));
}


double Testfun(double* x, int nd){
            double sum;
            sum =  pow(x[0] - 0.2, 2)* pow(x[1] - 0.3,2)*( 1 + 20*pow(x[0] + x[1] - 1.5, 2)*pow( x[0] + 2*(x)[1] - 2.3, 2));
            return sum;
}

double Schwef(double* x, int nd){

    int i;
    double sum;
    for ( i = 0; i < nd; i++){
        sum += x[i]*sin(sqrt(abs(x[i])));
    }

   sum = 418.9829*nd - sum;
  return sum;
}




fptr functionLookup(char* fName){
if (strcmp(fName,"Sphere")==0){return &Sphere;}
else if (strcmp(fName,"Rosenbrock")==0){return &Rosenbrock;}
else if (strcmp(fName,"Rastrigin")==0){return &Rastrigin;}
else if (strcmp(fName,"Griewank")==0){return &Griewank;}
else if (strcmp(fName,"Ackley")==0){return &Ackley;}
else if (strcmp(fName,"Levy")==0) {return &Levy;}
else if (strcmp(fName,"DropWave")==0){return &DropWave;}
else if (strcmp(fName,"Schwef")==0){return &Schwef;}
else if (strcmp(fName,"Testfun")==0){return &Testfun;}
else {printf("Function Name %s not recognized\n",fName);return NULL;}

}


