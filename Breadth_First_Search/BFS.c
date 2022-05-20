//
//  This program runs a Breadth First Search algorithm on a graph in order to find the size of the Giant component,
//  how many clusters there are in the graph, the minimal path from each node to the root node, the distance and who's the partent of each node
//  it takes into account the possibility that the graph can have more than one cluster obviously. 
//
//  Created by gino on 3/21/16.
//
// Compile: ./a.aout N initial_vertex File_data_graph 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#define directory "../Graphs" //directory containing the graph data
#define initial 0
#define waiting 1
#define visited 2

/* --- list to store the variable node ---- */
struct variable{
    
    int parent; // to store who's the parents
    int status; // visited, not visited
    int distance; // distance from the starting node
    int cluster;   // to which islands belongs to
    int degree;
    int *neigh;
};

/* --- list to store the queue ---- */
struct list{

    int val;
    struct list *next;

};

/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */

// "front" and "p_queue" are indexes of the same array queue[] where the queue is stored.
// "front" is used for the current node in the graph which is being visited "p_queue" refers to the last element added in the queue

// All the nodes from 0 to N-1 are stored in a linked list 0<-1<-2<-...< N-1 by using a list of structure called "curr"
// each node which is visited is deleted from this list. When a cluster of nodes in the graph has been
//completely visited, the function check_delete_list tells BFS from which node to restart in order to continue the searching in other clusters

// creates the linked list of nodes //
void create_list(struct list **head, struct list **curr, int N){
    
    int i,j;
    
    *head = NULL;
    
    for(i=0;i<N;i++){
        *curr = (struct list *)malloc(sizeof(struct list));
        (*curr)->val = i;
        (*curr)->next = *head;
        *head = *curr;
    }
    *curr = *head;
    
    return;
    
}


/***************************************************/

// Insert a new element in the queue array //
void insert_queue(int ver, int *queue, int *p_queue, int N){
    
    if(*p_queue > N-1){  // If the queue exceeds the space dedicated to it
        fprintf(stderr,"Queue overflow ----------- \n");
        exit(errno);
    }
    else{
        (*p_queue)++; // shift the index of the queue array of one position ahead 
        queue[*p_queue] = ver; // add a new element in the queue in the next position
    
        //    printf("   node %d in the queue[%d]\n",ver,*p_queue);
    }
    
    return;
}

/***************************************************/

// First part: remove the node from the list, second part: check if the current cluster is completed and if yes, start
// from a node in a new cluster
void check_delete_list(struct list **curr, int ver, int *front, int *p_queue, int *queue, int *cluster, int N){
    
    struct list *currP, *prevP;
    
    prevP = NULL;           // the preceding pointer to first element of the list, i.e. N, is NULL
    
    for(currP = *curr; currP != NULL; prevP = currP, currP = currP->next){
        if(currP->val == ver){       // if the node is the one to be removed
            if(prevP==NULL){         // if the node is the first of the list
                *curr = currP->next; // shift the pointer to next position
            }
            else                           // if the node is every other node and not the first of the list
                prevP->next = currP->next; // then his predecessor points to his descendant
            
            
            
            
            // Second part: If one cluster has been visited but there are still nodes to be visited //
            if(*front > *p_queue && *front < N){
                
                if(currP->next ==NULL)   // if the node is the last in the list
                    queue[*front]=prevP->val;
                else                    // if the node is every other node and not the first in the list
                    queue[*front]=(currP->next)->val;
                
                //    printf("cluster %d finished, moving to node %d\n",*cluster,queue[*front]);
                insert_queue(queue[*front],queue,p_queue,N); //put vertex ver into the queue
                (*cluster)++;
            }
            
                free(currP);    // free space of the item removed from the list
            return;
        }
    }
}


/***************************************************/

// Breadth first search function. Traverse the graph, visit nodes and mark them //
void BFS(struct variable *site, struct list *curr, int *queue, int ver, int N){

    int j,nn, count, front, p_queue, cluster;
    
    p_queue = -1;
    front = 0;
    count=0;
    cluster= 1;
    
    insert_queue(ver,queue,&p_queue,N); //put first vertex into the queue
    
    while(count < N){ // until you don't visit all the nodes
        
        site[ver].status = visited;  // mark the vertex as visited
        site[ver].cluster = cluster;
    
        for(j=0;j<site[ver].degree;j++){ // for all the neighbours of that node
            nn = site[ver].neigh[j];
         
            if(site[nn].status == initial){ // if the neighbours has not been visited yet, mark and add them to the queue
                
                insert_queue(nn,queue,&p_queue,N);    // put the neighbour of ver into the queue
                
                site[nn].status = waiting;                     // mark the vertex as visited
                site[nn].parent = ver;                         // assign the parents to nn
                site[nn].distance = site[ver].distance + 1;   // assign distance from original node to nn
                site[nn].cluster = cluster;                 // mark nn belongin to the cluster X
            }
        }
        
        front++;  // one node has been visited, shift ahead the front queue pointer
        check_delete_list(&curr,ver,&front,&p_queue,queue,&cluster,N); // remove the node which has been visited from the list and check if the cluster is complete
        ver = queue[front];   // next vertex to visit is the one at the front of the queue
        
        count++;                // count how many vertex have been visited
    }
    
    
    return;
}


/***************************************************/

// Read data graph from a file. Here J1 and J2 are useless. It's just that we have this data
void read_and_allocate(struct variable **site, int **queue, int N, char **argv){
    
    FILE *fp_file;
    char filename[101];
    int i,j,degree,nn;
    double J1,J2;
    
    *queue = (int *)malloc(N*sizeof(int));
    *site=(struct variable *)calloc(N,sizeof(struct variable));
    
    
    sprintf(filename,"%s_MC_graph.dat",argv[3]);
    fp_file = fopen(filename,"r");
    
    if(fp_file==NULL){
        fprintf(stderr,"Problem opening file with graph data MC_graph.dat \n");
        exit(errno);
    }
    
    for(i=0;i<N;i++){
        fscanf(fp_file,"%d",&degree);
        //     printf("degree %d = %d\n",i,degree);
        (*site)[i].neigh=(int *)malloc(N*sizeof(int));
        (*site)[i].degree = degree;
        (*site)[i].status = initial;    // initial status of the node = not visited
        (*site)[i].distance = 0;        //initial distance from the root node
        (*site)[i].cluster = 1;
        
        for(j=0;j<degree;j++){
            fscanf(fp_file,"%d%d%lf%lf",&i,&nn,&J1,&J2);
            (*site)[i].neigh[j]=nn;
            
            //        printf("%d ---> %d \t J[%d][%d] = %lf \t J[%d][%d] = %lf\n",i,nn,i,nn,J1,nn,i,J2);
        }
        
    }
    
    
    
    return;
}

/***************************************************/

// Print results after the whole graph has been traversed.
void print_structure(struct variable *site, int N){

    int i,n_cluster;
    n_cluster = 0;
    
    printf("\n");
    for(i=0;i<N;i++){
        
        if(site[i].cluster > 1){
            printf("i = %d, par_i = %d, status_i = %d, distance = %d, cluster = %d\n",i,site[i].parent,site[i].status,site[i].distance,site[i].cluster);
            n_cluster = 1;
        }
    }
        if(n_cluster == 0) printf("There is only one unique cluster of size %d\n\n",N);

}





/* =========================================================================== */
/*                                     MAIN                                    */
/* =========================================================================== */

int main(int argc, char *argv[]){

    chdir(directory); // Move to the directory
    
    int front,p_queue,N,ver;
    struct variable *site;
    struct list *curr, *head, *temp;
    int *queue,i,j,degree;
    
    N = atoi(argv[1]);   // numbers of nodes
    ver = atoi(argv[2]); // node from which to start traversing the graph
    
    read_and_allocate(&site,&queue,N,argv); 
    create_list(&head,&curr,N);
    
    BFS(site,curr,queue,ver,N);
    print_structure(site,N);

    return;

}





