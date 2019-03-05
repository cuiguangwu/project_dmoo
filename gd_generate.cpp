# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>

# include "gd_generate.h"

namespace az
{
//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
	//!\brief namespace of dynamic evolutionary algoirhtm
namespace dea
{

void gd:: generater (unsigned int sizenew, az::mea::CPopulationMO& popnew, az::mea::CPopulationMO& pop)
{
	 int i, j, k;
	population *parent_pop;
	list *temp, *temp1, *printf_env_list;

	int random, mark;
	double mux;
	double time1,start1,end1;

	nobj=pop.P().FSize();
	nreal=pop.P().XSize();
	popsize = sizenew;
	min_realvar = (double *)malloc(nreal*sizeof(double));
	max_realvar = (double *)malloc(nreal*sizeof(double));
	for (int i=0;i<nreal;i++)
	{
		min_realvar[i]=pop.P().XLow(i);
		max_realvar[i]=pop.P().XUpp(i);
	}
	//  随机种子
	srand((unsigned)time(NULL));
	random = rand()%1000;
	seed=(float) random/1000.0;

	pcross_real=0.9;

	mux = 1.0/nreal;
	pmut_real=mux;
	ncon=0;
	bitlength = 0;
	eta_c=15.0;
	eta_m=20.0;

	CR = 1.0;
	F = 0.2;
	CR_improve = 1.0;
	F_improve = 0.2;

	nrealmut = 0;
    nrealcross = 0;

	// 设置环境属性
	environment_div=10;                                                   // 环境每维单元域的数目
	environment_indcount=0;                                              // 环境个体总数
	environment_capacity=200;                                              // 环境容量
	environment_areacount=0;                                             // 环境中域总数

	
	environment_min = (double *)malloc(nobj*sizeof(double));             // 环境最小边界
	environment_max = (double *)malloc(nobj*sizeof(double));             // 环境最大边界

	extreme_ind = (int *)malloc(nobj*sizeof(int));                       // 每维上极端个体
	


	environment_distance = (double *)malloc(nobj*sizeof(double));        // 环境单位域大小

	environment_list = (list *)malloc(sizeof(list));                     // 环境中有个体的域
	environment_list->index = -1;
    environment_list->parent = NULL;
    environment_list->child = NULL;

	environment_dealwitharea_list = (list *)malloc(sizeof(list));        // 环境中待导向的域
    environment_dealwitharea_list->index = -1;
    environment_dealwitharea_list->parent = NULL;
    environment_dealwitharea_list->child = NULL;

	printf_env_list = (list *)malloc(sizeof(list));                      // 输出环境中个体
    printf_env_list->index = -1;
    printf_env_list->parent = NULL;
    printf_env_list->child = NULL;

	ultimate_list = (list *)malloc(sizeof(list));                         // 为变异存储临时个体链表
    ultimate_list->index = -1;
    ultimate_list->parent = NULL;
    ultimate_list->child = NULL;
	


	area = (unit_area *)malloc(pow((double)environment_div, nobj)*sizeof(unit_area));         // 分配域空间
    for (i=0; i<pow((double)environment_div, nobj); i++)                                      // 初始化域属性
	{
		allocate_memory_area (&area[i]);
		initialize_area (&area[i]);
	}
	


	// 设置外部种群
	
	parent_pop = (population *)malloc(sizeof(population));               // 分配种群空间
    allocate_memory_pop (parent_pop, 100*popsize);
	
	beyond_ind = (list *)malloc(sizeof(list));                           // 分配超个体链表空间
	beyond_ind->index = -1;
    beyond_ind->parent = NULL;
    beyond_ind->child = NULL;

	empty_ind = (list *)malloc(sizeof(list));                            // 分配空位个体空间
	empty_ind->index = -1;
    empty_ind->parent = NULL;
    empty_ind->child = NULL;

    beyond_boundary = (list *)malloc(sizeof(list));                      // 分配超界个体链表空间
	beyond_boundary->index = -1;
    beyond_boundary->parent = NULL;
    beyond_boundary->child = NULL;

	cur_ind = (list *)malloc(sizeof(list));                              // 记录当代个体
	cur_ind->index = -1;
    cur_ind->parent = NULL;
    cur_ind->child = NULL;
  

	randomize();
	
	for (k=0; k<1; k++)
	{
	    ngen = 0;
		eval_time = 0;

		for (i=0;i<popsize;i++)
		{
			for (j=0;j<nreal;j++)
			{
				parent_pop->ind[i].xreal[j]=pop[i][j];
			}
			parent_pop->ind[i].boundary_ind = 0;
			parent_pop->ind[i].constr_violation = 0.0;
		}
		for (i=0;i<popsize;i++)
		{
			for (j=0;j<nobj;j++)
			{
				parent_pop->ind[i].obj[j]=pop[i].F(j);
			}
		}
		//initialize_pop (parent_pop);                                         // 初始化外部种群
		//evaluate_pop (parent_pop);                                           // 评价外部种群
		//	
		//clear_environment(parent_pop);                                       // 清空环境
		
		temp = empty_ind->child;
		while (temp!=NULL)                                                   
		{
			temp = del (temp);
			temp = temp->child;
		}

		temp = cur_ind->child;
		while (temp!=NULL)                                                         // 清空临时域内个体
		{
			temp = del (temp);
			temp = temp->child;
		}		

		temp = beyond_ind->child;
		while (temp!=NULL)                                                   
		{
			temp = del (temp);
			temp = temp->child;
		}
		
		temp = empty_ind;                                                    // 初始化空位个体链表
		for (i=popsize; i<100*popsize; i++)                                    
		{
			insert (temp, i);
			temp = temp->child;
		}
		
		temp = beyond_ind;                                                   // 初始化超个体链表
		for (i=0; i<popsize; i++)
		{
			insert(temp, i);
			temp = temp->child;
		}
				
		construct_environment (parent_pop);                                  // 构造环境		
		
		temp = beyond_ind->child;                                            // 个体进入环境

		temp1 = cur_ind;
		while(temp!=NULL)
		{
			enter_environment (parent_pop, temp->index);
			insert (temp1, temp->index);
			temp = temp->child;
			temp1 = temp1->child;
		}
	    cur_indcount = popsize;

/*		
		fflush(fpt13);
		report_all(parent_pop, fpt13);
		fflush(fpt13);
*/		
		start1=clock();

		while(ngen<25)
//	    while(eval_time<50000)
		{				
		//	printf("\n\n !!! generation is %d!\n", ngen);
			for (i=0;i<popsize;i++)
			{
				for (j=0;j<nreal;j++)
				{
					pop[i][j]=parent_pop->ind[i].xreal[j];
				}
			}
			for (i=0;i<popsize;i++)
			{
				for (j=0;j<nobj;j++)
				{
					pop[i].F(j)=parent_pop->ind[i].obj[j];
				}
			}
			temp = beyond_ind->child;
			while (temp!=NULL)                                                   // 清空超域内个体
			{
				temp = del (temp);
				temp = temp->child;
			}
			
			temp = beyond_boundary->child;
			while (temp!=NULL)                                                   // 清空超界域
			{
				temp = del (temp);
				temp = temp->child;
			}
			
			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			temp = environment_list->child;                                      // SBX促进
			while (temp!=NULL)
			{
				temp1 = area[temp->index].area_ind->child;
				while (temp1!=NULL)
				{
					accelerate(parent_pop, temp1->index, temp->index, pop[0]);
					temp1 = temp1->child;
				}
				temp = temp->child;
			}
			
			temp = environment_list->child;                                      // 差分促进
			while (temp!=NULL)
			{
				temp1 = area[temp->index].area_ind->child;
				while (temp1!=NULL)
				{
					DE_accelerate(parent_pop, temp1->index, temp->index, pop[0]);
					temp1 = temp1->child;
				}
				temp = temp->child;
			}
			
			temp = environment_list->child;                                      // 环境中个体变异操作
			while (temp!=NULL)
			{
				temp1 = area[temp->index].area_ind->child;
				while(temp1!=NULL)
				{
					mutation(parent_pop, temp1->index, pop[0]); 
					temp1=temp1->child;
				}
				temp = temp->child;
			}
/*	
			temp = beyond_ind->child;                                            // 超个体链表复制到临时个体链表
			while (temp!=NULL)
			{
			insert(tempind_list, temp->index);
			temp = temp->child;
			}
			
			temp = tempind_list->child;                                        // 超链表中个体变异
			while(temp!=NULL)
			{
				mutation(parent_pop, temp->index);
				temp = temp->child;
			}
*/	

//		    clear_environment(parent_pop);                                       // 清空环境
			
/*
			if (beyond_judge==1)                                                 // 是否要重新构造环境
			{
				construct_environment (parent_pop);                                               
			}			  
*/      
	
			
			temp = beyond_ind->child;                                            // 个体进入环境
			while(temp!=NULL)
			{
				enter_environment (parent_pop, temp->index);
				temp = temp->child;
			}
      
			if (environment_indcount>environment_capacity)						// 可维护可不维护  
			{
				maintain_environment (parent_pop);                               
			}
			
/*		    printf("the number of individual in environment is %d!!!\n", environment_indcount); 
			printf("the number of area in environment is %d!!!\n", environment_areacount);
			temp = environment_list->child;
			i = 0;
			while (temp!=NULL&&i<30)
			{	
			printf("->%d", temp->index);
			temp = temp->child;
			i++;
			}
			
			  printf("the number of area in environment is %d!!!\n", environment_areacount);
			  if (i==30)
			  exit(1);
*/
/*
			fflush(fpt13);
			report_all(parent_pop, fpt13);
			fflush(fpt13);
*/		
			
			temp = beyond_ind->child;
			while (temp!=NULL)                                                  // 清空超域内个体
			{
				temp = del (temp);
				temp = temp->child;
			}
			
			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			lead_environment(parent_pop, pop[0]);							    // 环境导向操作

/*			if (beyond_judge == 1)
			{
				temp = environment_list->child;                                      
				while (temp!=NULL)
				{
					temp1 = area[temp->index].area_ind->child;
					while(temp1!=NULL)
					{
						insert(beyond_ind, temp1->index); 
						temp1=temp1->child;
					}
					temp = temp->child;
				}
				clear_environment(parent_pop);									// 清空环境
				construct_environment (parent_pop);
			}
*/			
			
			temp = beyond_ind->child;											// 个体进入环境
			while(temp!=NULL)
			{
				enter_environment (parent_pop, temp->index);
				temp = temp->child;
			}
/*
			fflush(fpt13);
			report_all(parent_pop, fpt13);
			fflush(fpt13);
*/
			
			if (environment_indcount>environment_capacity)
			{
				maintain_environment (parent_pop);								// 环境维护 
			}
/*			
			fflush(fpt13);
			report_all(parent_pop, fpt13);
			fflush(fpt13);
*/			
			
			temp = beyond_ind->child;
			while (temp!=NULL)													// 清空超域内个体
			{
				temp = del (temp);
				temp = temp->child;
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////			
			temp = ultimate_list->child;										// 清空最终链表内个体
			while (temp!=NULL)
			{
				temp = del (temp);
				temp = temp->child;
			}
			
			temp = environment_list->child;                                      
			while (temp!=NULL)
			{
				insert(ultimate_list, area[temp->index].area_repreasent_ind); 
				temp1 = area[temp->index].area_ind->child;
				while(temp1!=NULL)
				{
					insert(beyond_ind, temp1->index); 
					temp1=temp1->child;
				}
				temp = temp->child;
			}
			for (i=0; i<nobj; i++)
			{
				temp = ultimate_list->child;
				while(temp!=NULL&&extreme_ind[i]!=temp->index)
				{
					temp = temp->child;
				}
				if (temp==NULL)
				{
					insert(ultimate_list, extreme_ind[i]);
				}
			}
			
			
			temp = beyond_boundary->child;
			while (temp!=NULL)
			{
				insert(beyond_ind, temp->index);
				temp = temp->child;
			}

			popnew=pop;
			for (i=0;i<popsize;i++)
			{
				for (j=0;j<nreal;j++)
				{
					popnew[i][j]=parent_pop->ind[i].xreal[j];
				}
			}
			for (i=0;i<popsize;i++)
			{
				for (j=0;j<nobj;j++)
				{
					popnew[i].F(j)=parent_pop->ind[i].obj[j];
				}
			}

			//evaluate new solutions
			popnew.Evaluate();

			//environmental select
			pop.Combine(popnew);
	
			az::mea::sel::SCrowd2 sel;
			sel.Select(pop, 100);
			
			for (i=0;i<popsize;i++)
			{
				for (j=0;j<nreal;j++)
				{
					parent_pop->ind[i].xreal[j]=pop[i][j];
				}
				parent_pop->ind[i].boundary_ind = 0;
				parent_pop->ind[i].constr_violation = 0.0;
			}
			for (i=0;i<popsize;i++)
			{
				for (j=0;j<nobj;j++)
				{
					parent_pop->ind[i].obj[j]=pop[i].F(j);
				}
			}

			clear_environment(parent_pop);                                             // 清空环境
            
			construct_environment (parent_pop);                                        // 一代结束后重新构造环境
			
																					   
/*			temp = beyond_ind->child;                                                  // 删除被极端个体支配的个体、个体进入环境
			while(temp!=NULL)
			{
				mark = 0;
				for (i=0; i<nobj; i++)
				{
					flag = check_dominance(&parent_pop->ind[temp->index], &parent_pop->ind[extreme_ind[i]]);
					if (flag == -1)
					{
						mark = 1;
						break;
					}
				}
				if (mark != 1)
				{
					enter_environment (parent_pop, temp->index);
				}
				else
				{
					temp1 = empty_ind;                                                 // 进入空位链表
					insert(temp1, temp->index);
				}
				temp = temp->child;	
			}
*/			
			temp = cur_ind->child;
			while (temp!=NULL)                                                         // 清空临时域内个体
			{
				temp = del (temp);
				temp = temp->child;
			}		
			cur_indcount = 0;
			
			temp = beyond_ind->child;                                                  // 个体进入环境
			temp1 = cur_ind;
			while(temp!=NULL)
			{
				enter_environment (parent_pop, temp->index);
				insert(temp1, temp->index);
				temp = temp->child;
				temp1 = temp1->child;
				cur_indcount ++;
			}
//		    printf("cur_indcount is %d", cur_indcount);
/*		    fflush(fpt13);
		    report_all(parent_pop, fpt13);
		    fflush(fpt13);
*/		
			ngen++;
		}
		
		end1=clock();
		
		time1=(end1-start1)/CLOCKS_PER_SEC;
		
		printf("\n time is %f\n",time1);
		

	
	
/*	    temp = environment_list->child;                                                                               
		while (temp!=NULL)
		{
		insert(printf_env_list, area[temp->index].area_repreasent_ind); 
		temp = temp->child;
		}
*/	
//		report_nondomi_set (parent_pop);
//      report_nondomi_set_represent (parent_pop, ultimate_list);


		clear_environment(parent_pop);                                               // 清空环境
	}
	if (nreal!=0)
	{
		free (min_realvar);
		free (max_realvar);
	}
	free (environment_min);
	free (environment_max);
	free (extreme_ind);

	while (environment_list!=NULL)
	{
		temp = environment_list;
		environment_list = environment_list->child;
		free (temp);
	}
	while (printf_env_list!=NULL)
	{
		temp = printf_env_list;
		printf_env_list = printf_env_list->child;
		free (temp);
	}
	while (environment_dealwitharea_list!=NULL)
	{
		temp = environment_dealwitharea_list;
		environment_dealwitharea_list = environment_dealwitharea_list->child;
		free (temp);
	}
	while (beyond_ind!=NULL)
	{
		temp = beyond_ind;
		beyond_ind = beyond_ind->child;
		free(temp);
	}
	while (empty_ind!=NULL)
	{
		temp = empty_ind;
		empty_ind = empty_ind->child;
		free(temp);
	}
	while (cur_ind!=NULL)
	{
		temp = cur_ind;
		cur_ind = cur_ind->child;
		free(temp);
	}
	while (beyond_boundary!=NULL)
	{
		temp = beyond_boundary;
		beyond_boundary = beyond_boundary->child;
		free(temp);
	}

	for (i=0; i<pow((double)environment_div, nobj); i++)                                                                      
	{
		deallocate_memory_area (&area[i]);
	}
	free (area);

	//free (div_distance);

	deallocate_memory_pop (parent_pop, 100*popsize);
	//deallocate_memory_pop (child_pop, popsize);
	//deallocate_memory_pop (mixed_pop, 2*popsize);
	free (parent_pop);
	}
	}
	}
	}