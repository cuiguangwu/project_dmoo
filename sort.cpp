/* Routines for randomized recursive quick-sort */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Randomized quick sort routine to sort a population based on a particular objective chosen */
/* 对某一维目标值进行快速排序 */
void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size)
{
	q_sort_front_obj (pop, objcount, obj_array, 0, obj_array_size-1);
	return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
/* 对某一位目标值进行快速排序函数 */
void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right)
{
	int index;
	int temp;
	int i, j;
	double pivot;
	if (left<right)
	{
		index = rnd (left, right);
		temp = obj_array[right];
		obj_array[right] = obj_array[index];
		obj_array[index] = temp;
		pivot = pop->ind[obj_array[right]].obj[objcount];
		i = left-1;
		for (j=left; j<right; j++)
		{
			if (pop->ind[obj_array[j]].obj[objcount] <= pivot)
			{
				i+=1;
				temp = obj_array[j];
				obj_array[j] = obj_array[i];
				obj_array[i] = temp;
			}
		}
		index=i+1;
		temp = obj_array[index];
		obj_array[index] = obj_array[right];
		obj_array[right] = temp;
		q_sort_front_obj (pop, objcount, obj_array, left, index-1);
		q_sort_front_obj (pop, objcount, obj_array, index+1, right);
	}
	return;
}

/* Randomized quick sort routine to sort a population based on crowding distance */
/* 种群个体按聚集距离排序 */
void quicksort_dist(int *dist, int front_size)
{
	q_sort_dist (dist, 0, front_size-1);
	return;
}






/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
/* 种群个体按聚集距离排序函数 */
void q_sort_dist(int *dist, int left, int right)
{
	int index;
	int temp;
	int i, j;
	double pivot, pivot2;
	if (left<right)
	{
		index = rnd (left, right);
		temp = dist[right];
		dist[right] = dist[index];
		dist[index] = temp;
		pivot = area[dist[right]].area_indcount;
		pivot2 = area[dist[right]].area_nondominated_number;
		i = left-1;
		for (j=left; j<right; j++)
		{
			if (area[dist[j]].area_indcount < pivot)
			{
				i+=1;
				temp = dist[j];
				dist[j] = dist[i];
				dist[i] = temp;
			}
			if (area[dist[j]].area_indcount == pivot && area[dist[j]].area_nondominated_number >= pivot2)
			{
				i+=1;
				temp = dist[j];
				dist[j] = dist[i];
				dist[i] = temp;
			}
		}
		index=i+1;
		temp = dist[index];
		dist[index] = dist[right];
		dist[right] = temp;
		q_sort_dist (dist, left, index-1);
		q_sort_dist (dist, index+1, right);
	}
	return;
}
