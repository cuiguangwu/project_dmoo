/* This file contains the variable and function declarations */
#pragma once
#include "emo/IndividualMO.h"
# define INF 1.0e14
# define EPS 1.0e-14
# define E  2.71828182845905
# define PI 3.14159265358979
# define epsilon 1.0e-7
# define MYSIGN(x) ((x)>0?1.0:-1.0)

#define	REP_SIZE 1		//number of representives in each subpopulation
#define CLONE_NUM 1		//number of clone copies
#define VAR_NUMII 5
#define COMPETE "1"
#define SAMPLE 5
#define CR1 0.5
#define RETENT 5
#define ADIST 0.7 

#define OBJ_NUM 2
#define VAR_NUM 30
#define VAR_BIT 30

//EA PARAMETER
//#define SEED 6764421
#define SEED 100
//1670 360
#define GENS 100
#define	PC 0.8
//#define	PM 1.0/VAR_BIT
#define	PM1 1.0/VAR_BIT
#define POP1_SIZE 10
#define	POP2_SIZE 100
#define RADIUS 0.01

typedef struct
{
    int rank;						//rank值
    double constr_violation;		//约束违反
    double *xreal;					//决策（实数编码）
    int **gene;						//决策（二进制编码二进制串）
    double *xbin;					//决策（二进制编码）
    double *obj;					//目标向量
    double *constr;					//(带)约束向量
    double crowd_dist;				//聚集距离
	int boundary_ind;				//初始化为0，若该个体的某一目标值在该维上最小，并参与构造极端个体，则该值为1
} individual;

typedef struct
{
    individual *ind;
} population;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
} list;										//list节点类型，包含双向指针，分别指向上一个节点和下一个节点



typedef struct
{
    int area_exist;							//表示单元域是否存在（个体）
	int area_indcount;
	list *area_ind;							//单元域中个体队列
	int area_repreasent_ind;
	list *area_non_domianted;				//单元域中的非支配单元域
	int area_nondominated_number;
} unit_area;

extern int kkk;
extern int lll;
extern int nreal;							//实数编码决策维数
extern int nbin;							//二进制编码决策维数
extern int nobj;							//目标维数
extern int ncon;							//（带）约束维数
extern int popsize;
extern double pcross_real;
extern double pcross_bin;
extern double pmut_real;
extern double pmut_bin;
extern double eta_c;
extern double eta_m;
extern int ngen;
extern int nbinmut;
extern int nrealmut;
extern int nbincross;
extern int nrealcross;
extern int *nbits;
extern double *min_realvar;
extern double *max_realvar;
extern double *min_binvar;
extern double *max_binvar;
extern int bitlength;
extern int environment_div;				//每一维单元域数目
extern int environment_indcount;		//环境中个体数
extern int environment_areacount;		//环境中单元域数
extern int environment_capacity;
extern double *environment_min;			//环境（各维）下界
extern double *environment_max;			//环境（各维）上界
extern double *environment_distance;	//各单元域（各维）宽度
extern list *environment_list;			//环境单元域链表（指向链表头结点的指针）
extern list *environment_dealwitharea_list;
extern list *beyond_ind;				//越界个体链表（指向链表头结点的指针）
extern list *empty_ind;
extern list *beyond_boundary;
extern list *ultimate_list;
extern unit_area *area;
extern int beyond_judge;
extern int *extreme_ind;
extern int eval_time;
extern int cur_indcount;				//DE_accelate
extern list *cur_ind;					//DE_accelate
extern double CR;
extern double F;
extern double CR_improve;
extern double F_improve;

extern double hypervolume;

//allocate.cpp
//Memory allocation and deallocation routines 
void allocate_memory_pop (population *pop, int size);
void allocate_memory_ind (individual *ind);
void deallocate_memory_pop (population *pop, int size);
void deallocate_memory_ind (individual *ind);
void allocate_memory_area (unit_area *area);
void deallocate_memory_area (unit_area *area);


//auxiliary.cpp
double maximum (double a, double b);
double minimum (double a, double b);


//crossover.cpp
void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2);
void realcross (individual *parent1, individual *parent2, individual *child1, individual *child2);
void bincross (individual *parent1, individual *parent2, individual *child1, individual *child2);
void bincross1 (individual *parent1, individual *parent2, individual *child1, individual *child2);
void crossover1 (individual *parent1, individual *parent2, individual *child1, individual *child2);
void realcross1 (individual *parent1, individual *parent2, individual *child1, individual *child2);
void crossover_arithmetic (individual *parent1, individual *parent2, individual *child1, individual *child2);
void realcross_arithmetic (individual *parent1, individual *parent2, individual *child1, individual *child2);
void realcross_heuristic(individual *parent1, individual *parent2, individual *child1, individual *child2, double para);


//??
void assign_crowding_distance_list (population *pop, list *lst, int front_size);
void assign_crowding_distance_indices (population *pop, int c1, int c2);
void assign_crowding_distance (population *pop, int *dist, int **obj_array, int front_size);
//??
void decode_pop (population *pop);
void decode_ind (individual *ind);


//dominance.cpp
int check_dominance (individual *a, individual *b);
int check_dominance1 (individual *a, individual *b);
int check_area_dominance (int a, int b);


//??
void evaluate_pop (population *pop);
void evaluate_child (population *pop, int num);
//eval.cpp
void evaluate_ind (individual *ind, az::mea::CIndividualMO& eind);


//??
void fill_nondominated_sort (population *mixed_pop, population *new_pop, double *time1);
void crowding_fill (population *mixed_pop, population *new_pop, int count, int front_size, list *cur);


//initialize.cpp
void initialize_pop (population *pop);
void initialize_ind (individual *ind);
void initialize_area (unit_area *area);


//list.cpp
void insert (list *node, int x);
list* del (list *node);


//merge.cpp
void merge(population *pop1, population *pop2, population *pop3);
void copy_ind (individual *ind1, individual *ind2);


//mutation.cpp
void mutation_pop (population *pop);
void mutation_ind (individual *ind);
void bin_mutate_ind (individual *ind);
void real_mutate_ind (individual *ind);


//注释
void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr);


//??
void assign_rank_and_crowding_distance (population *new_pop);


//??
void report_pop (population *pop, FILE *fpt, int size);
void report_feasible (population *pop, FILE *fpt);
void report_ind (individual *ind, FILE *fpt);
void report_best(population *pop, FILE *fpt);
void report_time(double time, FILE *fpt);
void report_parent(individual *ind, FILE *fpt);
void report_child(population *pop, FILE *fpt, int num);
void report_best_list(population *pop, list *node, FILE *fpt);
void report_all(population *pop, FILE *fpt);
void report_best_reduce(population *pop, FILE *fpt, int obj_num);
void report_best_extent(population *pop, FILE *fpt, int count);


//sort.cpp
void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size);
void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right);
void quicksort_dist(int *dist, int front_size);
void q_sort_dist(int *dist, int left, int right);


//??
void selection (population *old_pop, population *new_pop);
//score.cpp
void score_pop(population *pop_final, int final_size);

individual* tournament (individual *ind1, individual *ind2);
individual* tournament1 (individual *ind1, individual *ind2);
void get_igd(population *pop, int final_size);
void get_test_points(double** obj,int num);

//score.cpp
void score_extent(population *pop_final, int final_size);
double calhypervolume(population *pop_final, int ind_num, int obj_num);
int filter_nondomi(population *pop_final, int ind_num, int obj_num);
int dominates(individual *ind1, individual *ind2, int obj_num) ;
void  swap(population *pop_final, int  i, int  j);
double  surfaceUnchangedTo(population *pop_final, int ind_num, int obj_num);
int  reduceNondominatedSet(population *pop_final, int ind_num, int obj_num, double threshold);

//env_set.cpp
void construct_environment (population *pop);
void enter_environment (population *pop, int index);
void enter_directly_area(int index, int number);
int judge_area(population *pop, int index, int number);
int judge_environment(population *pop, int index, int number);
list* findnode(list *node, int mark);
void change (list *node, int x);
void del_area(int number);
double cal_dist_origin(population *pop, int index, int number);
void accelerate(population *pop, int index, int number,az::mea::CIndividualMO& eind);
double cal_fitness(int a, int b);
void mutation (population* pop, int index,az::mea::CIndividualMO& eind);
void clear_environment(population *pop);
int choose_ind ();
void lead_environment(population *pop,az::mea::CIndividualMO& eind);
void maintain_environment(population *pop);
int compare_area(list *node, list *node1);
//??????????????????
void report_nondomi_set(population *pop);

void DE_accelerate(population *pop, int index, int number,az::mea::CIndividualMO& eind);
void DE_accelerate_improve(population *pop, int index, int number,az::mea::CIndividualMO& eind);

//diversity_score.cpp
void newscore_diversity(population *pop_new, int final_size);

//score.cpp
void find_extent(population *pop_final, int final_size);

//allocate.cpp
void allocate_memory_pop_rediv (population *pop, int size, int obj_num);
void allocate_memory_ind_rediv (individual *ind, int obj_num);

//?????????????????????
void report_nondomi_set_represent(population *pop, list *ultimate_list);
void report_best_final(population *pop, FILE *fpt, int final_size);
void get_gd(population *pop, int front_size);


int check_repeat (individual *a, individual *b);

//?????????????????????
void report_best_final_hv (population *pop, FILE *fpt, int final_size);
void get_test_points_gd(double** obj,int num);
void gd_igd_uf(population *pop, int front_size);
void gd_igd_wfg(population *pop, int front_size);
void gd_igd_cf(population *pop, int front_size);
void init_problem(int flag);

