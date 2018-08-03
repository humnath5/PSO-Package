
typedef double (*fptr)(double*,int);
fptr functionLookup(char* fName);

extern double * GBEST;
extern int ND;
extern int INDEX;
extern fptr FUN;


double fun_1D(double t);

void update_best(double * g, double * x, int nd);

void initialize_pso_domain( char * file, int* type, double* lb, double* ub ,double * gbest, int nd);

int lower_bound_check(double g, double l, double u, double epsilon);
int upper_bound_check(double g, double l, double u, double epsilon);

void adjust_domain(int* type, double* lb, double* ub, double * gbest, int nd, double alpha);

double brent_algorithm(double a, double b, double (*f) (double t), double tol);

void  local_optimizer(fptr fun, int nd, double* lb,
                       double* ub, double* value, double* gbest, int ni_stop, double tol, int myrank);

void mpipso_optimization(fptr fun, int np, int nd, int ni, double* lb, double* ub,
                      double* value, double* gbest, int ni_stop, double tol, int myrank);

void mpimulti_pso_optimization( fptr fun, int np, int nd, int ni, int maxRun, int* type, double* lb,
                          double* ub, double* value, double* gbest, int ni_stop, double tol, int myrank);

void mpimulti_pso_local_optimizer(fptr fun, int np, int nd, int ni, int maxRun,
                    int* type, double* lb, double* ub, double* value, double* gbest, int ni_stop, double tol, int myrank);
