#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include "randomlib.h"
#include "matmul2.h"
#include "pso.h"
#include "assert.h"


//#define MAX  5000
/* global variable */

#define numofgeom 282 // change this number if number of data points are different


/* Global Variables*/
double V[numofgeom],rH1I[numofgeom],rH2I[numofgeom],rOI[numofgeom],rGI[numofgeom];
double xH1[numofgeom],xH2[numofgeom],xO[numofgeom],xI[numofgeom],xG[numofgeom];
double yH1[numofgeom],yH2[numofgeom],yO[numofgeom],yI[numofgeom],yG[numofgeom];
double zH1[numofgeom],zH2[numofgeom],zO[numofgeom],zI[numofgeom],zG[numofgeom];
double approxV[numofgeom];

double fit_withghost(double , int);
double poten_energy(double *, double, double, double, double);
void read_coord_data(void);

int  main (int argc, char *argv[]){
int i,j,k;
int seed1,seed2;
srand((unsigned)time(NULL));
seed1= rand()%30000;
seed2= rand()%30000;
//seed1= 17795;
//seed2= 1795;

//seed1 = 17373;
//seed2 = 17496;

/* Reading command line inputs */
char func[80];
strcpy(func, argv[1]);
fptr f  = functionLookup(func);

int np     =  atoi(argv[2]);
int nd      =  atoi(argv[3]);
int ni     = atoi(argv[4]);
int maxRun = atoi(argv[5]);
int option = atoi(argv[6]);

read_coord_data();


/* If inter-atomic  distances are provided, we can use the following.
char file[20];
strcpy(file,"datafile.txt");
data = read_matrix(file,row,col);
*/

/* Alocationg memory spaces */
int* type = (int*)malloc(sizeof(int)*nd);
double* lb = (double*)malloc(sizeof(double)*nd);
double* ub = (double*)malloc(sizeof(double)*nd);
double * gbest = (double*)malloc(sizeof(double)*nd);/*gbest is global variable */

char bounds_file[40];
strcpy(bounds_file,"bounds_withghost.txt");
initialize_pso_domain(bounds_file,type,lb,ub,gbest,nd);
/*for( i = 0; i < nd; i ++){

   printf("%d\t %f\t %f\t %f\n",type[i],lb[i], ub[i], gbest[i]);
}
*/

 /*============Setting up MPI environment=============*/
    int myrank, comsize;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Barrier(MPI_COMM_WORLD);

    RandomInitialise(myrank + seed1, seed2);

/*Start Optimization */
double tol = 0.000001;
double value = 9999999;
int ni_stop = 100;
int ntries = 10;
/*Start Optimization */
double  start = MPI_Wtime(NULL);
if(option==0){
mpipso_optimization(f,np,nd,ni,lb,ub,&value,gbest,ni_stop,tol, myrank);}
else if(option==1){
mpimulti_pso_optimization(f,np,nd,ni,maxRun,type,lb,ub,&value,gbest,ni_stop,tol,myrank);}
else if(option==2){
mpimulti_pso_local_optimizer(f,np,nd,ni,maxRun,ntries,type,lb,ub,&value,gbest,ni_stop,tol,myrank);}
else{printf("\n No option is found !!!\n");}
double finish = MPI_Wtime(NULL);

double time=difftime(finish, start);

//double approxV[MAX];
for(i= 0;i< numofgeom;i++){
approxV[i]= poten_energy(gbest,rH1I[i],rH2I[i],rOI[i],rGI[i]);
}
//fprintf(fp,"%lf\t%lf\n",V[i],approxV[i]);}



gbest[7]=round(gbest[6]+gbest[7]),  gbest[6]=  round(gbest[6]);
gbest[15]=round(gbest[14]+gbest[15]); gbest[14]= round(gbest[14]);
gbest[23]=round(gbest[22]+gbest[23]);gbest[22]= round(gbest[22]);



/* Write outputs */
if(myrank == 0){
FILE* fp;
fp=fopen("output_withghost.txt","w");
fprintf(fp,"\n**================== PSO Parameters Used ==================**\n");
fprintf(fp,"seed1:%d\n",seed1);
fprintf(fp,"seed2:%d\n",seed2);
fprintf(fp,"dimension : %d\n",nd);
fprintf(fp,"num_particles: %d\n",np);
fprintf(fp,"max_num_iterations per run allowed: %d\n",ni);
fprintf(fp,"maxRun for multi-run PSO: %d\n",maxRun);
fprintf(fp,"total wall-clock time: %lf\n",time);
fprintf(fp,"best rmse : %lf\n",value);

fprintf(fp,"\n**================= Optimal Model Parameters =================**\n");
for(j=0;j<nd;j++){fprintf(fp,"%lf\n",gbest[j]);}

fprintf(fp,"\n**================= DFT Energy and Fitted Energy =================**\n");
//double approxV[MAX];
for(i= 0;i< numofgeom;i++){

 fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",rH1I[i],rH2I[i],rOI[i],rGI[i],V[i],approxV[i]);
}
fclose(fp);
}
MPI_Finalize();
free(lb);free(ub);free(gbest);free(type);
//free_matrix(data);

return 0;
}


void read_coord_data(void) {
 int i;
 FILE* input,*inputV;

    input  = fopen("geom.txt","r");
    inputV = fopen("energy.txt","r");
    char tempblank[1];

    for( i=0;i<numofgeom;i++)
    {
        fscanf(input,"%lf\t %lf\t %lf",&xO[i],&yO[i],&zO[i]);
        fscanf(input,"%lf\t%lf\t%lf",&xH1[i],&yH1[i],&zH1[i]);
        fscanf(input,"%lf\t%lf\t%lf",&xH2[i],&yH2[i],&zH2[i]);
        fscanf(input,"%lf\t%lf\t%lf",&xI[i],&yI[i],&zI[i]);
        fscanf(input,"%c",tempblank);
        fscanf(inputV,"%lf",&V[i]);


  /* calculating distances: */
        xG[i] = 0.00000; yG[i]= 0.00000; zG[i]= 0.00000;
        rH1I[i]= sqrt(pow((xH1[i]-xI[i]),2) +pow((yH1[i]-yI[i]),2)+pow((zH1[i]-zI[i]),2));
        rH2I[i]= sqrt(pow((xH2[i]-xI[i]),2)+pow((yH2[i]-yI[i]),2)+pow((zH2[i]-zI[i]),2));
        rOI[i] = sqrt(pow((xO[i] - xI[i]),2) +pow((yO[i]-yI[i]),2)+pow((zO[i]-zI[i]),2));
        rGI[i]= 0.0;

    }
  fclose(input);
  fclose(inputV);


}

double fit_withghost(double * x, int nd){/* nd is not necessary in this function but kept to be consistent for PSO input */
  int i;
     double rmse = 0;
     double sse = 0;

     for( i= 0; i< numofgeom; i++) {

          xG[i] =0.00000; yG[i]  = x[24]; zG[i] = 0.0000; // x[24] is alpha parameter

          rGI[i] = sqrt(pow((xG[i] - xI[i]),2) + pow((yG[i]-yI[i]),2)+pow((zG[i]-zI[i]),2));

          approxV[i] = poten_energy(x,rH1I[i],rH2I[i],rOI[i],rGI[i]);
          sse = sse + pow(approxV[i] - V[i],2);
       }

      rmse = sqrt(sse/numofgeom);
      return rmse;

}



double poten_energy(double * x, double r1, double r2, double r3, double r4){

     /* x -> parameter vector to be optimized
        r -> vector of inter-atomic distance */

    double  m1,m2,m3; /* integer parameters */
    double  n1,n2,n3; /* integer parameters */

    double pe = 0;
  //Calculating Nitrogen-Carbon potential energy

    m1  = round(x[6]);    n1 = round(x[6] + x[7]);

    m2  = round(x[14]);   n2 = round(x[14] + x[15]);

    m3  = round(x[22]);   n3 = round(x[22] + x[23]);


									// I-H1
            pe = pe + x[0] * exp(-x[1]*r1) + x[2]/ pow((r1+ x[4]), m1) + x[3]/ pow((r1+x[5]),n1);
                                    // I-H2
            pe = pe + x[0] * exp(-x[1]*r2) + x[2]/ pow((r2+ x[4]), m1) + x[3]/ pow((r2+x[5]),n1);
                                    // I-O
            pe = pe + x[8] * exp(-x[9]*r3) + x[10]/ pow((r3+ x[12]), m2) + x[11]/ pow((r3+x[13]),n2);
                                   // I-G
            pe = pe + x[16] * exp(-x[17]*r4) + x[18]/ pow((r4+ x[20]), m3) + x[19]/ pow((r4+x[21]),n3);

    return pe;

}



fptr functionLookup(char* fName){
/*if (strcmp(fName,"rosenbrock")==0){return &rosenbrock;}*/
if (strcmp(fName,"fit_withghost")==0){return &fit_withghost;}
else {printf("Function Name %s not recognized\n",fName);return NULL;}
}

