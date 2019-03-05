/* A custom doubly linked list implemenation */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Insert an element X into the list at location specified by NODE */
/* 将x插入到node之后 */
void insert (list *node, int x)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = (list *)malloc(sizeof(list));
    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != NULL)
    {
        node->child->parent = temp;
    }
    node->child = temp;
    return;
}
/* 将node后移至x之后 */
void change (list *node, int x)
{
	list *temp, *temp1;
	if (node->child==NULL)
    {
        printf("\n Error!!  NULL child pointer, hence exiting \n");
        exit(1);
    }
	do
	{		
		temp = node->child;
		temp1 = temp->child;
		node->parent->child = temp;
		temp->parent = node->parent;
		temp->child = node;		
		node->parent = temp;
		if (temp1!=NULL)
		{
			node->child = temp1;
			temp1->parent = node;
		}
		else 
		{
			node->child = NULL;
		}
	}
	while(temp->index!=x);
	return;
	
	

}

/* 返回index=mark的个体 */
list* findnode(list *node, int mark)
{
	list *temp;
    if (node==NULL)
    {
        printf("\n Error!!  NULL pointer, hence exiting \n");
        exit(1);
    }
	temp=node;
	while (temp!=NULL && temp->index!=mark)
	{
		temp=temp->child;
	}
	return temp;
}

/* 删除node结点，返回node结点上一个结点 */
list* del (list *node)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child!=NULL)
    {
        temp->child->parent = temp;
    }
    free (node);
    return (temp);
}
/*

void adjust_coordinate(int coordinate_x, int coordinate_y)          //找出代表个体
{
	list temp;
	double best;
	double distance[2];
	int mark;
	best = 100000;
	temp=area[coordinate_x][coordinate_y]->ind_list->child;
	while (temp!=NULL)
	{
		if (temp->index!=area[coordinate_x][coordinate_y].rep_index)
		{
		    distance[0] = env_set[temp->index]->obj[0]-(coordinate_x*div_distance[0]+min[0]);
		    distance[1] = env_set[temp->index]->obj[1]-(coordinate_y*div_distance[1]+min[1]);
            if (best>pow((distance[0]*distance[0]+distacne[1]*distance[1]), 0.5))
			{
				best=pow((distance[0]*distance[0]+distacne[1]*distance[1]), 0.5);
                mark=temp->index;
			}
		}
		temp=temp->child;
	}
	area[coordinate_x][coordinate_y].rep_index = mark;

    return;
}

int compare_rep(int coordinate_x, int coordinate_y, int mark)       //比较新个体与代表个体
{
	double distance[2];
	double rep_dis[2];
	int index;
	distance[0] = env_set[mark]->obj[0]-(coordinate_x*div_distance[0]+min[0]);
    distance[1] = env_set[mark]->obj[1]-(coordinate_y*div_distance[1]+min[1]);

	rep_dis[0]=env_set[area[coordinate_x][coordinate_y].rep_index]->obj[0]-(coordinate_x*div_distance[0]+min[0]);  
    rep_dis[1]=env_set[area[coordinate_x][coordinate_y].rep_index]->obj[1]-(coordinate_y*div_distance[1]+min[1]);
    if (pow((distance[0]*distance[0]+distance[1]*distance[1]), 0.5)<pow((rep_dis[0]*rep_dis[0]+rep_dis[1]*rep_dis[1]), 0.5))
        index=mark;
	else
		index=area[coordinate_x][coordinate_y].rep_index;
	
	return index;

}

int find_oldest(int coordinate_x, int coordinate_y)        //年龄最大的非代表个体
{
	list temp;
	int max_age=0;
	int mark_oldest;
	temp=area[coordinate_x][coordinate_y]->ind_list->child;
	while (temp!=NULL)
	{
        if (env_set[temp->index].age>max_age && temp->index!=area[coordinate_x][coordinate_y].rep_index)
		{
			max_age=env_set[temp->index].age;
			mark_oldest=temp->index;
		}
		temp=temp->child;
	}
	return mark_oldest;

}
*/