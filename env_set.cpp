// environment setting
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "global.h"
# include "rand.h"


//����Խ��������¹�������
void construct_environment (population *pop)
{
	int j, flag;
	list *temp;
	double *ind_max, *ind_min;
	ind_max = (double *)malloc(nobj*sizeof(double));
	ind_min = (double *)malloc(nobj*sizeof(double));

    
	// find boundary ind
	for (j=0; j<nobj; j++)
	{
		ind_min[j]=INF;
		ind_max[j]=EPS;
		temp = beyond_ind->child;
		while (temp!=NULL)
		{
			if (pop->ind[temp->index].obj[j]<ind_min[j])
			{
				ind_min[j] = pop->ind[temp->index].obj[j];
				extreme_ind[j] = temp->index;
			}
			if (pop->ind[temp->index].obj[j]==ind_min[j])
			{
				flag = check_dominance(&pop->ind[temp->index], &pop->ind[extreme_ind[j]]);
				if (flag == 1)//������1֧�����2
				{
					extreme_ind[j] = temp->index;
				}
			} 
			if (pop->ind[temp->index].obj[j]>ind_max[j])
				ind_max[j] = pop->ind[temp->index].obj[j];
			temp = temp->child;
		}
		pop->ind[extreme_ind[j]].boundary_ind = 1;
		environment_distance[j]=(ind_max[j]-ind_min[j])*(environment_div+1)/(environment_div*environment_div);//area_size(j)
		environment_min[j] = ind_min[j]-(environment_distance[j]*environment_div-(ind_max[j]-ind_min[j]))/2;
        environment_max[j] = ind_max[j]+(environment_distance[j]*environment_div-(ind_max[j]-ind_min[j]))/2;
 
	}

	free(ind_max);
	free(ind_min);
	return;
	
}

//��������еĵ�Ԫ��
void clear_environment(population *pop)
{
	list *temp;
	int i;
	temp = environment_list->child;
	while(temp!=NULL)//���ɾ�������еĵ�Ԫ��
	{
		del_area(temp->index);//ɾ����Ԫ��
		temp = del(temp);     //�ڻ���������ɾ����Ӧ���
		temp = temp->child;
	}
	temp = environment_dealwitharea_list->child;//?????????????????????????????
	while(temp!=NULL)
	{
		temp = del(temp);
	//	temp = temp->index;
	}
	environment_indcount = 0;
	environment_areacount = 0;
	for (i=0; i<nobj; i++)
	{
		pop->ind[extreme_ind[i]].boundary_ind = 0;//��ʼ�����˵�
	}

	return;
}


//��indcount>capacityʱ���޳�����ֱ������Ҫ��
void maintain_environment(population *pop)
{
	int *dist;
	list *temp, *temp1, *temp2;
	int i, mark, count, random_int, delete_ind, index;
//	printf("maintain_environment\n");							 
	dist = (int *)malloc(environment_areacount*sizeof(int));
	temp = environment_list->child;
	for (i=0; i<environment_areacount; i++)
	{
		dist[i] = temp->index;
		temp = temp->child;
	}
	//�Ի�������ݾۼ���������,����dist��˳�����environment_listָ��ĵ�ַ
    quicksort_dist (dist, environment_areacount);
	temp = environment_list->child;
	while(temp!=NULL)
	{
		temp = del(temp);
		temp = temp->child;
	}
	temp = environment_list;
	for (i=0; i<environment_areacount; i++)
	{
		insert(temp, dist[i]);
	}

/*	temp = environment_list->child;
	while (temp!=NULL)
	{	
		printf("->%d\t", temp->index);
		temp = temp->child;
	}
	printf("\n");
*/	
                                                                                      // �޳������и���
	while (environment_indcount>environment_capacity)
	{
		temp = environment_list->child;
		if (area[temp->index].area_indcount==1)                                      // һ�㲻������������
		{
			printf("�����������������ã���\n");
			temp1 = area[temp->index].area_ind->child;
			insert(empty_ind, temp1->index);		
			del_area(temp->index);
			environment_areacount--;		
			temp = del(temp);
			environment_indcount--;
		}
		else
		{
			                                                                         
			temp1 = area[temp->index].area_ind->child;
			while (temp1!=NULL&&pop->ind[temp1->index].boundary_ind==1)		      //�߽���岻ɾ��
			{
				temp1=temp1->child;
			}
			if (temp1==NULL)                                                         // ����ȫΪ���˸������������򽻻�����environment_list�������
			{
				printf("���������ĳ����ˣ���");
				temp1 = temp;                                                       
				while (temp1->child!=NULL)
				{
					temp1=temp1->child;
				}
				change(temp, temp1->index);			
			}
			else
			{
				do 
				{
					random_int = rnd(1, area[temp->index].area_indcount);
					temp1 = area[temp->index].area_ind->child;
					count = 1;
					while (temp1!=NULL&&count<random_int)
					{
						count++;
						temp1=temp1->child;
					}
					if (pop->ind[temp1->index].boundary_ind == 1)                         // ����ɾ�����˸���
					{
						temp1 = area[temp->index].area_ind->child;
						while (temp1!=NULL&&pop->ind[temp1->index].boundary_ind==1)
						{
							temp1=temp1->child;
						}
						break;
					}
				} 
				while(area[temp->index].area_repreasent_ind==temp1->index);                     
				insert(empty_ind, temp1->index);
				temp1 = del(temp1);
				area[temp->index].area_indcount--;
				environment_indcount--;
				temp1 = temp;
				mark = 0;
				// if 1, represent temp better; else return 0
				while (temp1->child!=NULL&&compare_area(temp, temp1->child)==1)
				{
					temp1 = temp1->child;
					mark = 1;
				}
				if (mark == 1)
				{
					change(temp, temp1->index);							
				}
/*				index = 0;
				temp2 = environment_list->child;
				while (temp2!=NULL&&index!=20)
				{	
					printf("->%d\t", temp2->index);
					temp2 = temp2->child;
					index++;
				}
				printf("index is %d\n ", index);
*/				
			}
												
		}

	}
	free (dist);
//	printf("exit maintain_environment\n"); 
	return;
}


//������뻷��
void enter_environment (population *pop, int index)
{
	int *ind_coordinate;
	int i, j;
	list *temp, *temp1, *temp2;
	int number, judgement, end, flag, flag1, extreme, extreme1;
	double distance1, distance2;
//  printf("enter_environment\n");
	for (j=0; j<nobj; j++)                                                      // �ж��뼫�˸���Ĺ�ϵ
	{
		flag = check_dominance (&pop->ind[index], &pop->ind[extreme_ind[j]]);
		if (flag == -1)//�����˸���֧��
		{
			temp = empty_ind;                                                   // �����λ����
			insert(temp, index);
			return;
		}
	}

                                                                                // �жϸ����Ƿ�Խ��
    for (j=0; j<nobj; j++)
	{
		if (pop->ind[index].obj[j]<environment_min[j]||pop->ind[index].obj[j]>environment_max[j])
		{
			printf("����Խ�磡����");
			exit(1);
		}		
	}

    ind_coordinate = (int *)malloc(nobj*sizeof(int));
	for (j=0; j<nobj; j++)
	{
		ind_coordinate[j] = (int)floor((pop->ind[index].obj[j]-environment_min[j])/environment_distance[j]);
//		printf("ind_coordinate[%d] is %d\t", j, ind_coordinate[j]);
	}
/*	printf("\n");
    for (j=0; j<nobj; j++)
	{
		printf("%e\t", pop->ind[index].obj[j]);
	}
	printf("\n");
*/
	
	number = 0;
	for (j=0; j<nobj; j++)                                                      // ��������ڻ�����λ��
	{
		number += ind_coordinate[j]*pow((double)environment_div, (nobj-1-j));
	}
	judgement = 0;
	if (environment_indcount==0)                                                // �������в����ڸ���ʱ
	{
		
		enter_directly_area(index, number);                                     // ֱ�ӽ�����

		temp = environment_list;                                                // ����������Ϣ
		insert(temp, number);
		environment_indcount++;
		environment_areacount++;
	}
	else                                                                        // �������д��ڸ���ʱ
	{
		// judgement is 1 represents ind should enter environment, -1 represents non-enter
		judgement = judge_environment(pop, index, number);                      // �жϸ����Ƿ������
		if (judgement == -1)
		{
			temp = empty_ind;                                                   // �����λ����
			insert(temp, index);
			free(ind_coordinate);
			return;
		}
		else
		{
			                                                                    // ������뻷��
			temp = environment_list->child;
			end = 0;
			do
			{	                                                              
				flag = check_area_dominance (number, temp->index);             
			    if (flag == 3)													//�������Ѿ����ڣ��򽫸ø�����õ�Ԫ���еĸ�����бȽϣ���������֮�󣬽��ø������õ�Ԫ��
				{
					temp1 = area[number].area_ind->child;
					while (temp1!=NULL)
					{
						flag1 = check_dominance (&pop->ind[index], &pop->ind[temp1->index]);
		                switch(flag1)															
						{
		                    case 1:																//1�����õ�Ԫ���еĸ��壨1�����ø��壨2��֧�䣬�򽫸��壨1��������ɾ�����������λ����
							{
				                temp2 = empty_ind;                              
								insert(temp2, temp1->index);
								temp1 = del(temp1);
								area[number].area_indcount--;
								environment_indcount--;
								temp1 = temp1->child;
								break;
							}
		                    case 0:																//0:���õ�Ԫ���еĸ��壨1����ø��壨2������֧�䣬������
							{
				                temp1 = temp1->child;
				                break;
							}
							case -1:
							{
								printf("this situation should not be happen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
								temp1 = temp1->child;
							}
							case 2:
							{
								printf("this situation should not be happen too!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
								temp1 = temp1->child;									
							}
						}			
					}
					area[number].area_indcount++;
					temp1 = area[number].area_ind;
					insert(temp1, index);
					                                                                					
					flag1 = check_dominance(&pop->ind[index], &pop->ind[area[number].area_repreasent_ind]);       // �Ƚϸ�����������      
					if(flag1 == 1)
						area[number].area_repreasent_ind = index;
					else
					{
						distance1 = cal_dist_origin(pop, index, number);
						distance2 = cal_dist_origin(pop, area[number].area_repreasent_ind, number);
						if (distance1 <= distance2)
							area[number].area_repreasent_ind = index;
					}
						
					
					end = 1;				
				}
				temp = temp->child;
					
			}
			while (temp!=NULL&&end!=1);
			if (end == 0)                                                                   // ���������µ�����
			{
				area[number].area_exist = 1;
				area[number].area_indcount++;
				temp = area[number].area_ind;
				insert (temp, index);
				area[number].area_repreasent_ind = index;
				temp = environment_list->child;
				do
				{			                                                              
					flag = check_area_dominance (number, temp->index);              // ������뻷�����򡢻�����Ϣ����
					switch(flag)
					{
					case 1:																	//���µĵ�Ԫ��֧�����еĵ�Ԫ���򽫱�֧��ĵ�Ԫ��ɾ�����������������λ����
						{
							environment_indcount=environment_indcount-area[temp->index].area_indcount;
							temp1 = area[temp->index].area_ind->child;
							while(temp1!=NULL)
							{
								temp2 = empty_ind;
								insert(temp2, temp1->index);
								temp1 = temp1->child;
							}
							del_area(temp->index);
							environment_areacount--;
							temp = del(temp);
							temp =temp->child;
							break;
						}
					case 2:																	//���µĵ�Ԫ��*��֧�����еĵ�Ԫ��ִ����ͬ�Ĳ���
						{
							environment_indcount=environment_indcount-area[temp->index].area_indcount;
							temp1 = area[temp->index].area_ind->child;
							while(temp1!=NULL)
							{
							temp2 = empty_ind;
							insert(temp2, temp1->index);
							temp1 = temp1->child;
							}
							del_area(temp->index);
							environment_areacount--;
							temp = del(temp);
							temp =temp->child;
							break;
						}
					case -2:																//���µĵ�Ԫ�����еĵ�Ԫ��*��֧�䣬��ֱ������
						{
							temp = temp->child;
							break;
						}				
					case 0:																	//���µĵ�Ԫ�������еĵ�Ԫ�򻥲�֧�䣬��������Ԫ���ڸ��Եķ�֧�䵥Ԫ�������м���Է�
						{
							temp1 = area[temp->index].area_non_domianted;
							insert(temp1, number);
							area[temp->index].area_nondominated_number++; 
							temp1 = area[number].area_non_domianted;
							insert(temp1, temp->index);
							area[number].area_nondominated_number++;
							temp = temp->child;
							break;
						}
					}
				}
				while (temp!=NULL);
				temp = environment_list;
			    insert(temp, number);											//����µĵ�Ԫ������environment_list
				environment_areacount++;
			}
			for (j=0; j<nobj; j++)												//�뼫�˸���Ƚ�
			{
				if (pop->ind[index].obj[j]<pop->ind[extreme_ind[j]].obj[j])
				{
					pop->ind[extreme_ind[j]].boundary_ind = 0;
					extreme_ind[j] = index;			
					pop->ind[index].boundary_ind = 1;					
				}
				if (pop->ind[index].obj[j]==pop->ind[extreme_ind[j]].obj[j])
				{
					flag = check_dominance(&pop->ind[index], &pop->ind[extreme_ind[j]]);
					if (flag == 1)
					{
						pop->ind[extreme_ind[j]].boundary_ind = 0;
						extreme_ind[j] = index;
					    pop->ind[index].boundary_ind = 1;
					}
/*					if (flag == 0)
					{
						if (randomperc() <= 0.5)
						{
							pop->ind[extreme_ind[j]].boundary_ind = 0;
							extreme_ind[j] = index;
							pop->ind[index].boundary_ind = 1;
						}
					}				
*/				}
			}
			environment_indcount ++;
/*			temp = environment_list->child;
			while(temp!=NULL)
			{
				printf("->%d(%d, %d)", temp->index, area[temp->index].area_indcount, area[temp->index].area_nondominated_number);
				
				temp=temp->child;
			}
			temp = environment_list->child;
			while(temp!=NULL)
			{
				
				printf("\narea[%d] non-dominated is\t", temp->index);
				temp1 = area[temp->index].area_non_domianted->child;
				while (temp1!=NULL)
				{
					printf("->%d", temp1->index);
					temp1 = temp1->child;
				}
				temp=temp->child;
			}
			printf("\n");
*/		
		}

	}
//	printf("individual in environment is %d\n", environment_indcount);
//	printf("area in environment is %d\n", environment_areacount);
//    printf("exit enter_environment\n");
    free(ind_coordinate);
}


//�жϸ����Ƿ���Խ��뻷��
int judge_environment(population *pop, int index, int number)
{
	list *temp;
	int flag, end;
	
	temp = environment_list->child;
	end = 0;
	do
	{   
		// a dominate b return 1; a weak dominate b return 2; b dominate a return -1; b weak dominate a return -2; a = b return 3; a non-dominate b return 0
		flag = check_area_dominance (number, temp->index);
		switch(flag)
		{
		case 1:
			{
				return 1;
			}
		case -1:
			{
				return -1;
			}
		case -2:
			{
			//	return -1;
				temp = temp->child;
				break;
			}
		case 3:
			{
				end = 1;
				break;
			}
		case 0:
			{
				temp = temp->child;
				break;
			}
		case 2:
			{
				temp = temp->child;
				break;
			}
		}

	}
	while (end!= 1 && temp!=NULL);
	if (end == 1)
	{
		// �ж��Ƿ������ 1 represent enter -1 represent non-enter
		flag = judge_area(pop, index, number);
		if (flag == 1)
			return 1;
		else 
			return -1;
		
	}
	else
		return 1;
}

//crossover
void accelerate(population *pop, int index, int number, az::mea::CIndividualMO& eind)
{
	list *temp;
	int i, j, k, chosen_area, random_int, cross_ind, flag, flag1, flag2, end;
	double sum, *area_fitness, prob;
    population *child;
	int mark[2];
	
	child = (population *)malloc(sizeof(population));               // ������Ⱥ�ռ�
    allocate_memory_pop (child, 2);
    
//	printf("accelerate\n");

	if (area[number].area_nondominated_number==0)                   // �������ڷ�֧����ʱ ���޸ģ����ѡ��
	{
		if (environment_indcount>1)
		{
			do 
			{
				cross_ind = choose_ind();
			} 
			while(index==cross_ind);
			printf("this area has no non-dominated area......\n");
		}
		else
		{
			printf("this environment has only one individual!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			DE_accelerate_improve(pop, index, number,eind);
			deallocate_memory_pop(child, 2);
        	free (child);
//			printf("exit accelerate\n");
			return;
//			exit(1);			
		}
	}
	else
	{
	    area_fitness = (double *)malloc(area[number].area_nondominated_number*sizeof(double));
	    i = 0;
	    sum = 0;
	    temp = area[number].area_non_domianted->child;				//	�������fitness
	    while (temp!=NULL)
		{
		    area_fitness[i] = cal_fitness(number, temp->index);
		    sum += area_fitness[i];
		    i++;
		    temp = temp->child;
		}															//�������fitness��ʹ�����̶ķ�����non_dominated list��ѡ��һ����Ԫ��>chosen_area
	    for (i=0; i<area[number].area_nondominated_number; i++)
		{
		    area_fitness[i] = area_fitness[i]/sum;
		}
	    for (i=1; i<area[number].area_nondominated_number; i++)
		{
		    area_fitness[i] = area_fitness[i-1]+area_fitness[i];
		}
        prob = randomperc();

        temp = area[number].area_non_domianted->child;
	
	    i=0;
	    while(temp!=NULL&&prob>area_fitness[i])
		{
		    i++;
		    temp=temp->child;
		}
	    chosen_area = temp->index;
	


	    random_int = rnd(1, area[chosen_area].area_indcount);		//����ѡ��ķ�֧������ѡ��һ������
	    temp = area[chosen_area].area_ind->child;
	    i = 1;
	    while (temp!=NULL&&i<random_int)
		{
		    i++;
		    temp=temp->child;
		}
	    cross_ind = temp->index;
		free(area_fitness);
	}

	crossover (&pop->ind[index], &pop->ind[cross_ind], &child->ind[0], &child->ind[1]);
	evaluate_ind(&child->ind[0], eind);
	evaluate_ind(&child->ind[1], eind);
	
	mark[0]=1;
	mark[1]=1;
	//ֻ����֧����߻���֧��ĸ��壨�����������Ӵ����壩
	flag = check_dominance(&child->ind[0], &child->ind[1]);
	if (flag == 1)
		mark[1] = 0;
	if (flag == -1||flag == 2)
		mark[0] = 0;

	for (i=0; i<2; i++)
	{
		if (mark[i]==1)
		{
			flag1 = check_dominance(&child->ind[i], &pop->ind[cross_ind]);
            flag2 = check_dominance(&child->ind[i], &pop->ind[index]);
			if (flag1==-1||flag2==-1||flag2==2||flag1==2)
				mark[i] = 0;
		}
	}


	for (i=0; i<2; i++)                                                       
	{
		if (mark[i]==1)
		{
			temp = empty_ind->child;
	        copy_ind (&child->ind[i], &pop->ind[temp->index]);                // ���λ��
			end = 0;
			for (j=0; j<nobj; j++)
			{
				if (child->ind[i].obj[j]>environment_max[j]||child->ind[i].obj[j]<environment_min[j])         
				    end = 1;
			}
			if (end == 1)
			{
				insert (beyond_boundary, temp->index);  }                      // �����峬�������߽磬����beyond_boundary��
			else
				insert (beyond_ind, temp->index);                              // ������δ���������߽磬����beyond_ind��
			temp = del(temp);                                                  
		}
	}
    
	deallocate_memory_pop(child, 2);
	free (child);
//	printf("exit accelerate\n");
	return;
}


//DE����
void DE_accelerate(population *pop, int index, int number, az::mea::CIndividualMO& eind)
{
	int i, j, k;
    int random_int, chosen_area, mark1, mark2, j_rand, flag, end, cross_ind;
	double pc;
	list *temp;
	individual* temp_ind;
//  printf("DE_accelerate\n");
	if (environment_indcount==1)
	{
		DE_accelerate_improve(pop, index, number,eind);
//		printf("exit DE_accelerate\n");
		return;
	}
	temp_ind = (individual*)malloc(sizeof(individual));
	allocate_memory_ind (temp_ind);
    
	do 
	{
        random_int = rnd(1, environment_areacount);
		temp = environment_list->child;
	    i = 1;
	    while (temp!=NULL&&i<random_int)
		{
		    i++;
		    temp=temp->child;
		}
	    chosen_area = temp->index;
		random_int = rnd(1, area[chosen_area].area_indcount);
	    temp = area[chosen_area].area_ind->child;
	    i = 1;
	    while (temp!=NULL&&i<random_int)
		{
		    i++;
		    temp=temp->child;
		}
	    cross_ind = temp->index;

	}
	while(index==cross_ind);

	if (cur_indcount<4)
	{
		do 
		{
			mark1 = rnd(0, popsize-1);
		} 
		while(mark1 == cross_ind || mark1 == index);
		
		do 
		{
			mark2 = rnd(0, popsize-1);
		} 
		while(mark2 == cross_ind || mark2 == index || mark2 == mark1);
	}
	else
	{
		do 
		{
			random_int = rnd(1, cur_indcount);
			i = 1;
			temp = cur_ind->child;
			while (temp!=NULL&&i<random_int)
			{
				i++;
				temp=temp->child;
			}
			mark1 = temp->index;
		} 
		while(mark1 == cross_ind || mark1 == index);
		
		do 
		{
			random_int = rnd(1, cur_indcount);
			i = 1;
			temp = cur_ind->child;
			while (temp!=NULL&&i<random_int)
			{
				i++;
				temp=temp->child;
			}
			mark2 = temp->index;
		} 
		while(mark2 == cross_ind || mark2 == index || mark2 == mark1);

	}

	

	j_rand = rnd(0, nreal - 1);

	for (k = 0; k < nreal; k ++)
	{
		pc = rndreal (0, 1);
		if (pc < CR || k == j_rand)			
		{
			temp_ind->xreal[k] = pop->ind[cross_ind].xreal[k] + F * (pop->ind[mark1].xreal[k] - pop->ind[mark2].xreal[k]);
		}
		else
			temp_ind->xreal[k] = pop->ind[index].xreal[k];
		
		if (temp_ind->xreal[k] > max_realvar[k])
		{
			temp_ind->xreal[k] = max_realvar[k];
		}
		else
		{
			if (temp_ind->xreal[k] < min_realvar[k])
			{
				temp_ind->xreal[k] = min_realvar[k];
			}
		}
	}

	evaluate_ind(temp_ind, eind);
	flag = check_dominance(temp_ind, &(pop->ind[index]));
	if (flag == 1 || flag == 0)
	{
		temp = empty_ind->child;
		copy_ind (temp_ind, &pop->ind[temp->index]);                          // ���λ��
		end = 0;
		for (j=0; j<nobj; j++)
		{
			if (temp_ind->obj[j]>environment_max[j]||temp_ind->obj[j]<environment_min[j])         
				end = 1;
		}
		if (end == 1)
			insert (beyond_boundary, temp->index);                            // �����峬�������߽� ���볬��������
		else
			insert (beyond_ind, temp->index);                                 // ������δ���������߽� ���볬����
		temp = del(temp);   
	}

    deallocate_memory_ind(temp_ind);
	free(temp_ind);
//	printf("exit DE_accelerate\n");
	return;
}


//DE�������������壩
void DE_accelerate_improve(population *pop, int index, int number, az::mea::CIndividualMO& eind)
{
	int i, j, k;
    int random_int, chosen_area, mark1, mark2, j_rand, flag, end, cross_ind;
	double pc;
	list *temp;
	individual* temp_ind;
//  printf("DE_accelerate_improve\n");
	temp_ind = (individual*)malloc(sizeof(individual));
	allocate_memory_ind (temp_ind);

  if (cur_indcount<4)
	{
		do 
		{
			cross_ind = rnd(0, popsize-1);	
		}		
		while(index==cross_ind);
		
		do 
		{
			mark1 = rnd(0, popsize-1);
		} 
		while(mark1 == cross_ind || mark1 == index);
		
		do 
		{
			mark2 = rnd(0, popsize-1);
		} 
		while(mark2 == cross_ind || mark2 == index || mark2 == mark1);

	}
	else
	{
		do 
		{
			random_int = rnd(1, cur_indcount);
			i = 1;
			temp = cur_ind->child;
			while (temp!=NULL&&i<random_int)
			{
				i++;
				temp=temp->child;
			}
			cross_ind = temp->index;
			
		}
		while(index==cross_ind);
		
		do 
		{
			random_int = rnd(1, cur_indcount);
			i = 1;
			temp = cur_ind->child;
			while (temp!=NULL&&i<random_int)
			{
				i++;
				temp=temp->child;
			}
			mark1 = temp->index;
		} 
		while(mark1 == cross_ind || mark1 == index);
		
		do 
		{
			random_int = rnd(1, cur_indcount);
			i = 1;
			temp = cur_ind->child;
			while (temp!=NULL&&i<random_int)
			{
				i++;
				temp=temp->child;
			}
			mark2 = temp->index;
		} 
		while(mark2 == cross_ind || mark2 == index || mark2 == mark1);
	}

	j_rand = rnd(0, nreal - 1);

	for (k = 0; k < nreal; k ++)
	{
		pc = rndreal (0, 1);
		if (pc < CR_improve || k == j_rand)			
		{
			temp_ind->xreal[k] = pop->ind[index].xreal[k] + F_improve * (pop->ind[mark1].xreal[k] - pop->ind[mark2].xreal[k]);
		}
		else
			temp_ind->xreal[k] = pop->ind[cross_ind].xreal[k];
		
		if (temp_ind->xreal[k] > max_realvar[k])
		{
			temp_ind->xreal[k] = max_realvar[k];
		}
		else
		{
			if (temp_ind->xreal[k] < min_realvar[k])
			{
				temp_ind->xreal[k] = min_realvar[k];
			}
		}
	}

	evaluate_ind(temp_ind, eind);
	flag = check_dominance(temp_ind, &(pop->ind[index]));
	if (flag == 1 || flag == 0)
	{
		temp = empty_ind->child;
		copy_ind (temp_ind, &pop->ind[temp->index]);                          // ���λ��
		end = 0;
		for (j=0; j<nobj; j++)
		{
			if (temp_ind->obj[j]>environment_max[j]||temp_ind->obj[j]<environment_min[j])         
				end = 1;
		}
		if (end == 1)
			insert (beyond_boundary, temp->index);                            // �����峬�������߽� ���볬��������
		else
			insert (beyond_ind, temp->index);                                 // ������δ���������߽� ���볬����
		temp = del(temp);   
	}

    deallocate_memory_ind(temp_ind);
	free(temp_ind);
//	printf("exit DE_accelerate_improve\n");
	return;
}


//�������
void mutation (population* pop, int index, az::mea::CIndividualMO& eind)
{
	individual *ind;
	int flag, flag1;
	int j, k, end;
	list *temp;
	end=0;
	ind = (individual*)malloc(sizeof(individual));
	allocate_memory_ind(ind);

	copy_ind (&pop->ind[index], ind); 
	mutation_ind(ind);
	evaluate_ind(ind, eind);

    flag = check_dominance(ind, &pop->ind[index]);
	if (flag!=-1&&flag!=2)
	{
		temp = empty_ind->child;
		copy_ind (ind, &pop->ind[temp->index]);                               // ���λ��
		for (j=0; j<nobj; j++)
		{
			if (ind->obj[j]>environment_max[j]||ind->obj[j]<environment_min[j])         
				end = 1;
		}
		if (end == 1)
			insert (beyond_boundary, temp->index);                            // �����峬�������߽� ���볬��������
		else
			insert (beyond_ind, temp->index);                                 // ������δ���������߽� ���볬����
		temp = del(temp);
	}
	
	deallocate_memory_ind(ind);
	free (ind);
	return;

}


//�����������
void lead_environment(population *pop, az::mea::CIndividualMO& eind)
{
	int i, j, k, m, n, end;
	int **record;
    list *temp, *temp1, *temp2;
	int **a_coordinate, *ind_coordinate;
	int remainder, mark1, mark2, flag1, flag2, cross_count, nearest_number, chosen_fitness, flag;
	double *low_fitness, *high_fitness;
	int *low_fitness_int, *high_fitness_int;
	list *low_list, *high_list, *nearest_list;
	int count_low, count_high, nearest_distance;
	double sum_low, sum_high, prob, cross_para;
	int chosen_area1, chosen_area2, random_int1, random_int2, cross_ind1, cross_ind2, number;
	individual *ind1, *ind2;

//	printf("lead_environment\n");

	if (environment_areacount==1)                   // ��������ֻ����һ���и������ʱ
	{	
		printf("this environment has only one area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		for (i=0; i<nobj*environment_div; i++)
		{		
			if (environment_indcount==1)
			{
				temp = environment_list->child;
				temp1 = area[temp->index].area_ind->child;
				DE_accelerate_improve(pop, temp1->index, temp->index,eind);
			}
			else
			{
				temp = environment_list->child;
				random_int1 = rnd(1, area[temp->index].area_indcount);
				temp1 = area[temp->index].area_ind->child;
				n = 1;
				while (temp1!=NULL&&n<random_int1)
				{
					n++;
					temp1=temp1->child;
				}
				DE_accelerate_improve(pop, temp1->index, temp->index,eind);
			}
		}
		return;
	}
	
	ind_coordinate = (int *)malloc(nobj*sizeof(int));

	ind1 = (individual*)malloc(sizeof(individual));
	allocate_memory_ind(ind1);

	ind2 = (individual*)malloc(sizeof(individual));
	allocate_memory_ind(ind2);

	low_list = (list *)malloc(sizeof(list));                     
	low_list->index = -1;
    low_list->parent = NULL;
    low_list->child = NULL;

	high_list = (list *)malloc(sizeof(list));                     
	high_list->index = -1;
    high_list->parent = NULL;
    high_list->child = NULL;

	nearest_list = (list *)malloc(sizeof(list));                     
	nearest_list->index = -1;
    nearest_list->parent = NULL;
    nearest_list->child = NULL;

	low_fitness = (double *)malloc(environment_areacount*sizeof(double));
    high_fitness = (double *)malloc(environment_areacount*sizeof(double));

	low_fitness_int = (int *)malloc(environment_areacount*sizeof(int));
    high_fitness_int = (int *)malloc(environment_areacount*sizeof(int));
	
	record = (int**)malloc(nobj*sizeof(int));
    a_coordinate = (int **)malloc(environment_areacount*sizeof(int));
	for (j=0; j<environment_areacount; j++)
	{
		a_coordinate[j] = (int *)malloc(nobj*sizeof(int));
	}

    for (j=0; j<nobj; j++)
	{
		record[j]=(int*)malloc(environment_div*sizeof(int));
	}
    for (j=0; j<nobj; j++)
	{
		for (i=0; i<environment_div; i++)
		{
			record[j][i]=0;
		}
	}
    
	k = 0;
	temp = environment_list->child;                                                       // �ҳ�ÿά�д��������
	while (temp!=NULL)
	{
	    remainder = temp->index;
	    for (j=0; j<nobj; j++)
		{
		    if (j == nobj-1)
			{
				a_coordinate[k][j] = remainder;
                record[j][a_coordinate[k][j]] = 1;
			}
			    
		    else
			{
			    a_coordinate[k][j] = (int)floor(remainder/pow((double)environment_div, nobj-1-j));
				record[j][a_coordinate[k][j]] = 1;
		        remainder = remainder%(int)pow((double)environment_div, nobj-1-j);
			}			
		}
		k++;
		temp = temp->child;
	}

    

	for (j=0; j<nobj; j++)                                                               // �������
	{
		for (i=0; i<environment_div; i++)
		{
			if (record[j][i]==0)
			{
				k = 0;
				count_low = 0;
				count_high = 0;
				sum_low = 0;
				sum_high = 0;
				
				temp = low_list->child;
				while(temp!=NULL)
				{
					temp = del(temp);
					temp = temp->child;
				}
				temp = high_list->child;
				while(temp!=NULL)
				{
					temp = del(temp);
					temp = temp->child;
				}

				temp = nearest_list->child;
				while(temp!=NULL)
				{
					temp = del(temp);
					temp = temp->child;
				}
			
	            nearest_distance = environment_div;
				temp = environment_list->child;
				temp1 = low_list;
				temp2 = high_list;
				while(temp!=NULL)
				{
					
					if (a_coordinate[k][j]<i)
					{
						insert(temp1, temp->index);
						low_fitness[count_low] = 1.0/(i-a_coordinate[k][j]);
						sum_low += low_fitness[count_low];
						low_fitness_int[count_low] = i-a_coordinate[k][j];
						if ((i-a_coordinate[k][j])<nearest_distance)
						{
							nearest_distance = i-a_coordinate[k][j];							
						}
						temp1 = temp1->child;
						count_low++;
					}
					else
					{
						insert(temp2, temp->index);
                        high_fitness[count_high] = 1.0/(a_coordinate[k][j]-i);
						sum_high += high_fitness[count_high];
						high_fitness_int[count_high] = a_coordinate[k][j]-i;
						if ((a_coordinate[k][j]-i)<nearest_distance)
						{
							nearest_distance = a_coordinate[k][j]-i;							
						}
						temp2 = temp2->child;
					
						count_high++;
					}
					k++;
					temp = temp->child;
				}
				
				if (count_low==0||count_high==0)                                       // ����������߻����ұ���û�и���
				{
				//	printf("�㷺�Դ���������");
					if (count_high==0)                                                   
					{
						n=0;                                                           
						nearest_number = 0;
						temp = low_list->child;
						while(temp!=NULL)
						{
							if (low_fitness_int[n]==nearest_distance)
							{
								insert(nearest_list, temp->index);
								nearest_number++;
							}
							n++;
							temp = temp->child;								
						}
						cross_count = 0;
						do 
						{
						    cross_count++;                                                // ѡ���һ����
							random_int1 = rnd(1, nearest_number);
							temp = nearest_list->child;
							n = 1;
							while (temp!=NULL&&n<random_int1)
							{
								n++;
								temp=temp->child;
							}
							chosen_area1 = temp->index;
																						  // ѡ��ڶ�����  (��������ֻ��һ����ʱ����д)
							do 
							{
								random_int2 = rnd(1, count_low);
								temp = low_list->child;
								n = 1;
								while (temp!=NULL&&n<random_int2)
								{
									n++;
									temp=temp->child;
								}
								chosen_area2 = temp->index;
								chosen_fitness = low_fitness_int[n-1];
                        		
							} 
							while(chosen_area2==chosen_area1);
																						  // ��������ʽ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							if (chosen_fitness==nearest_distance)
								cross_para = nearest_distance;
							else
								cross_para = (double)nearest_distance/(double)(chosen_fitness-nearest_distance);
/*							
							if ((chosen_fitness-nearest_distance)<=nearest_distance)
								cross_para = 1.0;
							else
								cross_para = (double)nearest_distance/(double)(chosen_fitness-nearest_distance);
*/																						  // ѡ������ĸ���
							random_int1 = rnd(1, area[chosen_area1].area_indcount);
							temp = area[chosen_area1].area_ind->child;
							n = 1;
							while (temp!=NULL&&n<random_int1)
							{
								n++;
								temp=temp->child;
							}
							cross_ind1 = temp->index;

							random_int2 = rnd(1, area[chosen_area2].area_indcount);
							temp = area[chosen_area2].area_ind->child;
							n = 1;
							while (temp!=NULL&&n<random_int2)
							{
								n++;
								temp=temp->child;
							}
							cross_ind2 = temp->index;

							realcross_heuristic(&pop->ind[cross_ind1], &pop->ind[cross_ind2], ind1, ind2, cross_para);			
							evaluate_ind(ind1, eind);
							evaluate_ind(ind2, eind);
																										
							mark1 = 0;
							for (m=0; m<nobj; m++)
							{
								if (ind1->obj[m]>environment_max[m]||ind1->obj[m]<environment_min[m])
								{
									mark1 = -1;
									break;
								}
							}
							if (mark1 == -1)
							{
								flag1 = check_dominance(ind1, &pop->ind[cross_ind1]);
								flag2 = check_dominance(ind1, &pop->ind[cross_ind2]);
								if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
								{
									temp = empty_ind->child;
									copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_boundary, temp->index);                            
									temp = del(temp);
								}
							}

							if (mark1!=-1)
							{
								for (m=0; m<nobj; m++)
								{
									ind_coordinate[m] = (int)floor((ind1->obj[m]-environment_min[m])/environment_distance[m]);       // ��������ڻ�����λ��
								}
								number = 0;
								for (m=0; m<nobj; m++)                                                      
								{
									number += ind_coordinate[m]*pow((double)environment_div, (nobj-1-m));
								}
								flag1 = check_area_dominance (number, chosen_area1);
								flag2 = check_area_dominance (number, chosen_area2);
								
								if (ind_coordinate[j]==i&&flag1!=-1&&flag2!=-1)
								{
									mark1 = 1;
									temp = empty_ind->child;
									copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_ind, temp->index);                                      // �������
									temp = del(temp);                                                      // �������
								}
								else
								{
									flag1 = check_dominance(ind1, &pop->ind[cross_ind1]);
									flag2 = check_dominance(ind1, &pop->ind[cross_ind2]);
									if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
									{
										temp = empty_ind->child;
										copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
										insert (beyond_ind, temp->index);                                      // �������
										temp = del(temp);                                                      // �������
									}

								}

							}

							mark2 = 0;
							for (m=0; m<nobj; m++)
							{
								if (ind2->obj[m]>environment_max[m]||ind2->obj[m]<environment_min[m])
								{
									mark2 = -1;
									break;
								}
							}
							if (mark2 == -1)
							{
								flag1 = check_dominance(ind2, &pop->ind[cross_ind1]);
								flag2 = check_dominance(ind2, &pop->ind[cross_ind2]);
								if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
								{
									temp = empty_ind->child;
									copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_boundary, temp->index);                            
									temp = del(temp);
								}
							}

							if (mark2!=-1)
							{
								for (m=0; m<nobj; m++)
								{
									ind_coordinate[m] = (int)floor((ind2->obj[m]-environment_min[m])/environment_distance[m]);
								}
								number = 0;
								for (m=0; m<nobj; m++)                                                      
								{
									number += ind_coordinate[m]*pow((double)environment_div, (nobj-1-m));
								}
								flag1 = check_area_dominance (number, chosen_area1);
								flag2 = check_area_dominance (number, chosen_area2);
								
								if (ind_coordinate[j]==i&&flag1!=-1&&flag2!=-1)
								{
									mark2 = 1;
									temp = empty_ind->child;
									copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_ind, temp->index);                                      // �������
									temp = del(temp);                                                      // �������
								}
								else
								{
									flag1 = check_dominance(ind2, &pop->ind[cross_ind1]);
									flag2 = check_dominance(ind2, &pop->ind[cross_ind2]);
									if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
									{
										temp = empty_ind->child;
										copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
										insert (beyond_ind, temp->index);                                      // �������
										temp = del(temp);                                                      // �������
									}
								}
							}
						}
						while(cross_count<environment_div&&mark1!=1&&mark2!=1);


					
					
					}
                    else                                                   
					{
						n=0;                                                           
						nearest_number = 0;
						temp = high_list->child;
						while(temp!=NULL)
						{
							if (high_fitness_int[n]==nearest_distance)
							{
								insert(nearest_list, temp->index);
								nearest_number++;
							}
							n++;
							temp = temp->child;								
						}
                        cross_count = 0;
                        do
						{
							cross_count++;                                                 // ѡ���һ����
							random_int1 = rnd(1, nearest_number);
							temp = nearest_list->child;
							n = 1;
							while (temp!=NULL&&n<random_int1)
							{
								n++;
								temp=temp->child;
							}
							chosen_area1 = temp->index;
																						   // ѡ��ڶ�����
							do 
							{
								random_int2 = rnd(1, count_high);
								temp = high_list->child;
								n = 1;
								while (temp!=NULL&&n<random_int2)
								{
									n++;
									temp=temp->child;
								}
								chosen_area2 = temp->index;
								chosen_fitness = high_fitness_int[n-1];
                        		
							} 
							while(chosen_area2==chosen_area1);
																						  // ��������ʽ�������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							if ((chosen_fitness-nearest_distance)<=nearest_distance)
								cross_para = 1.0;
							else
								cross_para = (double)nearest_distance/(double)(chosen_fitness-nearest_distance);
																						  // ѡ������ĸ���
							random_int1 = rnd(1, area[chosen_area1].area_indcount);
							temp = area[chosen_area1].area_ind->child;
							n = 1;
							while (temp!=NULL&&n<random_int1)
							{
								n++;
								temp=temp->child;
							}
							cross_ind1 = temp->index;

							random_int2 = rnd(1, area[chosen_area2].area_indcount);
							temp = area[chosen_area2].area_ind->child;
							n = 1;
							while (temp!=NULL&&n<random_int2)
							{
								n++;
								temp=temp->child;
							}
							cross_ind2 = temp->index;

							realcross_heuristic(&pop->ind[cross_ind1], &pop->ind[cross_ind2], ind1, ind2, cross_para);			
							evaluate_ind(ind1, eind);
							evaluate_ind(ind2, eind);
																										
							mark1 = 0;
							for (m=0; m<nobj; m++)
							{
								if (ind1->obj[m]>environment_max[m]||ind1->obj[m]<environment_min[m])
								{
									mark1 = -1;
									break;
								}
							}
							if (mark1 == -1)
							{
								flag1 = check_dominance(ind1, &pop->ind[cross_ind1]);
								flag2 = check_dominance(ind1, &pop->ind[cross_ind2]);
								if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
								{
									temp = empty_ind->child;
									copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_boundary, temp->index);                            
									temp = del(temp);
								}
							}

							if (mark1!=-1)
							{
								for (m=0; m<nobj; m++)
								{
									ind_coordinate[m] = (int)floor((ind1->obj[m]-environment_min[m])/environment_distance[m]);
								}
								number = 0;
								for (m=0; m<nobj; m++)                                                      
								{
									number += ind_coordinate[m]*pow((double)environment_div, (nobj-1-m));
								}
								flag1 = check_area_dominance (number, chosen_area1);
								flag2 = check_area_dominance (number, chosen_area2);
								
								if (ind_coordinate[j]==i&&flag1!=-1&&flag2!=-1)
								{
									mark1 = 1;
									temp = empty_ind->child;
									copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_ind, temp->index);                                      // �������
									temp = del(temp);                                                      // �������
								}
								else
								{
									flag1 = check_dominance(ind1, &pop->ind[cross_ind1]);
									flag2 = check_dominance(ind1, &pop->ind[cross_ind2]);
									if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
									{
										temp = empty_ind->child;
										copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
										insert (beyond_ind, temp->index);                                      // �������
										temp = del(temp);                                                      // �������
									}

								}

							}

							mark2 = 0;
							for (m=0; m<nobj; m++)
							{
								if (ind2->obj[m]>environment_max[m]||ind2->obj[m]<environment_min[m])
								{
									mark2 = -1;
									break;
								}
							}
							if (mark2 == -1)
							{
								flag1 = check_dominance(ind2, &pop->ind[cross_ind1]);
								flag2 = check_dominance(ind2, &pop->ind[cross_ind2]);
								if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
								{
									temp = empty_ind->child;
									copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_boundary, temp->index);                            
									temp = del(temp);
								}
							}

							if (mark2!=-1)
							{
								for (m=0; m<nobj; m++)
								{
									ind_coordinate[m] = (int)floor((ind2->obj[m]-environment_min[m])/environment_distance[m]);
								}
								number = 0;
								for (m=0; m<nobj; m++)                                                      
								{
									number += ind_coordinate[m]*pow((double)environment_div, (nobj-1-m));
								}
								flag1 = check_area_dominance (number, chosen_area1);
								flag2 = check_area_dominance (number, chosen_area2);
								
								if (ind_coordinate[j]==i&&flag1!=-1&&flag2!=-1)
								{
									mark2 = 1;
									temp = empty_ind->child;
									copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_ind, temp->index);                                      // �������
									temp = del(temp);                                                      // �������
								}
								else
								{
									flag1 = check_dominance(ind2, &pop->ind[cross_ind1]);
									flag2 = check_dominance(ind2, &pop->ind[cross_ind2]);
									if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
									{
										temp = empty_ind->child;
										copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
										insert (beyond_ind, temp->index);                                      // �������
										temp = del(temp);                                                      // �������
									}
								}
							}
						}
						while(cross_count<environment_div&&mark1!=1&&mark2!=1);			
					}
				}                
				else                                                                   // �����������ұ߶��и���
				{
               //     printf("�����Դ���������");
					for (m=0; m<count_low; m++)
					{
						low_fitness[m] = low_fitness[m]/sum_low;
					}
					for (m=1; m<count_low; m++)
					{
						low_fitness[m] = low_fitness[m-1]+low_fitness[m];
					}
					for (m=0; m<count_high; m++)
					{
						high_fitness[m] = high_fitness[m]/sum_high;
					}
					for (m=1; m<count_high; m++)
					{
						high_fitness[m] = high_fitness[m-1]+high_fitness[m];
					}


					cross_count = 0;
					do
					{
						cross_count++;
						prob = randomperc();        
						temp1 = low_list->child;	
						m=0;
						while(temp1!=NULL&&prob>low_fitness[m])
						{
							m++;
							temp1=temp1->child;
						}
						chosen_area1 = temp1->index;
						
						prob = randomperc();        
						temp2 = high_list->child;	
						m=0;
						while(temp2!=NULL&&prob>high_fitness[m])
						{
							m++;
							temp2=temp2->child;
						}
						chosen_area2 = temp2->index;
						
						random_int1 = rnd(1, area[chosen_area1].area_indcount);
						temp = area[chosen_area1].area_ind->child;
						n = 1;
						while (temp!=NULL&&n<random_int1)
						{
							n++;
							temp=temp->child;
						}
						cross_ind1 = temp->index;

						random_int2 = rnd(1, area[chosen_area2].area_indcount);
						temp = area[chosen_area2].area_ind->child;
						n = 1;
						while (temp!=NULL&&n<random_int2)
						{
							n++;
							temp=temp->child;
						}
						cross_ind2 = temp->index;

						realcross_arithmetic(&pop->ind[cross_ind1], &pop->ind[cross_ind2], ind1, ind2);	         // ��������		
						evaluate_ind(ind1, eind);
						evaluate_ind(ind2, eind);
						
						mark1 = 0;
                        for (m=0; m<nobj; m++)
						{
							if (ind1->obj[m]>environment_max[m]||ind1->obj[m]<environment_min[m])
							{
								mark1 = -1;
								break;
							}
						}
						if (mark1 == -1)
						{
							flag1 = check_dominance(ind1, &pop->ind[cross_ind1]);
							flag2 = check_dominance(ind1, &pop->ind[cross_ind2]);
							if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
							{
								temp = empty_ind->child;
								copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
								insert (beyond_boundary, temp->index);                            
								temp = del(temp);
							}
						}

						if (mark1!=-1)
						{
							for (m=0; m<nobj; m++)
							{
								ind_coordinate[m] = (int)floor((ind1->obj[m]-environment_min[m])/environment_distance[m]);
							}
							number = 0;
							
							for (m=0; m<nobj; m++)                                                      
							{
								number += ind_coordinate[m]*pow((double)environment_div, (nobj-1-m));
							}
							flag1 = check_area_dominance (number, chosen_area1);
							flag2 = check_area_dominance (number, chosen_area2);
							if (ind_coordinate[j]==i&&flag1!=-1&&flag2!=-1)
							{
								mark1 =1;
								temp = empty_ind->child;
								copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
								insert (beyond_ind, temp->index);                                      // �������
								temp = del(temp);                                                      // �������
							}
							else
							{
								flag1 = check_dominance(ind1, &pop->ind[cross_ind1]);
								flag2 = check_dominance(ind1, &pop->ind[cross_ind2]);
								if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
								{
									temp = empty_ind->child;
									copy_ind (ind1, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_ind, temp->index);                                      // �������
									temp = del(temp);                                                      // �������
								}

							}

						}

						mark2 = 0;
						for (m=0; m<nobj; m++)
						{
							if (ind2->obj[m]>environment_max[m]||ind2->obj[m]<environment_min[m])
							{
								mark2 = -1;
								break;
							}
						}
						if (mark2 == -1)
						{
							flag1 = check_dominance(ind2, &pop->ind[cross_ind1]);
							flag2 = check_dominance(ind2, &pop->ind[cross_ind2]);
							if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
							{
								temp = empty_ind->child;
								copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
								insert (beyond_boundary, temp->index);                            
								temp = del(temp);
							}
						}
						if (mark2!=-1)
						{
							for (m=0; m<nobj; m++)
							{
								ind_coordinate[m] = (int)floor((ind2->obj[m]-environment_min[m])/environment_distance[m]);
							}
							number = 0;
								
							for (m=0; m<nobj; m++)                                                      
							{
								number += ind_coordinate[m]*pow((double)environment_div, (nobj-1-m));
							}
							flag1 = check_area_dominance (number, chosen_area1);
							flag2 = check_area_dominance (number, chosen_area2);
							
							if (ind_coordinate[j]==i&&flag1!=-1&&flag2!=-1)
							{
								mark2 =1;
								temp = empty_ind->child;
								copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
								insert (beyond_ind, temp->index);                                      // �������
								temp = del(temp);                                                      // �������
							}
							else
							{
								flag1 = check_dominance(ind2, &pop->ind[cross_ind1]);
								flag2 = check_dominance(ind2, &pop->ind[cross_ind2]);
								if (flag1 != -1 && flag2 != -1 && flag1 != 2 && flag2 != 2)
								{
									temp = empty_ind->child;
									copy_ind (ind2, &pop->ind[temp->index]);                               // ���λ��
									insert (beyond_ind, temp->index);                                      // �������
									temp = del(temp);                                                      // �������
								}

							}

						}
					}
					while(cross_count<environment_div&&mark1!=1&&mark2!=1);
				}			
			}
		}
	}
	free(ind_coordinate); 
    deallocate_memory_ind(ind1);
	deallocate_memory_ind(ind2);
	while (low_list!=NULL)
	{
		temp = low_list;
		low_list = low_list->child;
		free(temp);
	}
	while (high_list!=NULL)
	{
		temp = high_list;
		high_list = high_list->child;
		free(temp);
	}
	while (nearest_list!=NULL)
	{
		temp = nearest_list;
		nearest_list = nearest_list->child;
		free(temp);
	}
	free(low_fitness);
	free(high_fitness);
	free(low_fitness_int);
	free(high_fitness_int);
	

	for (j=0; j<environment_areacount; j++)
	{
		free(a_coordinate[j]);
	}
	free(a_coordinate);
	for (j=0; j<nobj; j++)
	{
		free(record[j]);
	}
	free(record);

//	printf("exit lead_environment\n");

	return;

}
