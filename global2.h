// error processing
#define ERROR(x) fprintf(stderr, x), fprintf(stderr, "\n"), exit(1)

// set the reference point as externer variable 
extern const double point_ref[];
extern int dimension;

void Swap(double *front[], int i, int j);
void copy(double *temp[], double *front[], int i, int j);
void Deallocate_memrory(double *front[], int noPoints);
void allocate_memory(double **front[], int noPoints);
void Merge_Sort(double *front[], int p, int r, int objectives);
void Sort(double *front[], int noPoints, int objective);
static void Merge(double *front[], int p, int q, int r, int objectives);
int ReduceFront(double *front[], int noPoints, int noObjectives, double *slice_depth, int slice_number);
int Dominates(double *point1, double *point2, int noObjectives);
int FilterNondominatedSet(double *front[], int noPoints, int noObjectives);
int Check_Dominate(double *front[], int noPoints, int objective);
int Is_In_Set(double *temp[], int temp_length, double *point);
int ReadFile(double **front[], FILE *fpt, int noObjectives);
double CalculateHypervolume(double *front[], int noPoints, int noObjectives);
double CalculateArea(double *front[], int noPoints, int noObjectives);
double getBound(double *front[], int noPoints, int i, int j);
int IsRepeaded(int x, int array[], int size);
void Extract(double **front[], FILE *fpt, int noPoints, int noObjectives);
double CalculateArea_2D(double *front[], int noPoints, int noObjectives);
