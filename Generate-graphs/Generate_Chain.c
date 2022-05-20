//
//  Generate_ERRG.c
//  
//
//  Created by gino on 1/20/15.
//
//  Compile: ./a.out N  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define directory "DMP_data" //directory containing the data

struct variable{ //Structure for MC simulations
    int degree; // number of neighbours of a given site
    int *neigh; //neighbours of a given site
    int *J;     // couplings J_{ij}
};

void convert_parameter(int argc, char **argv);
void allocate_memory(struct variable **site);
void generate_open_chain(struct variable *site,int J1, int J11, int J2, int J22); // generate a CHAIN with OPEN boundary conditions
void generate_periodic_chain(struct variable *site,int J1, int J11, int J2, int J22); // generate a CHAIN with PERIODIC boundary conditions
void print_graph_on_file(struct variable *site, char **argv);   // print graph on 2 files

// Rand numb generator
double get_rand(void);

int N;


/* =========================================================================== */
/*                                  MAIN                                       */
/* =========================================================================== */

int main(int argc, char **argv){
    
    struct variable *site;
    int i,j;
  
    
    chdir(directory); // Move to the directory
    srand(time(NULL)); //initialize get_rand(), random numb generator
    
    convert_parameter(argc,argv);
    allocate_memory(&site);
   // generate_open_chain(site,1,0,1,0);
    generate_periodic_chain(site,1,1,1,1);
    print_graph_on_file(site,argv);
    

    return 0;
}




/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */


void convert_parameter(int argc, char **argv){
    
    if(argc>1){
        
        N=atoi(argv[1]);
        
    }

    return;

}

/***************************************************/

void allocate_memory(struct variable **site){
    
    int i;
    
    /* allocate memory for the array of structures */
    *site=(struct variable *)malloc(N*sizeof(struct variable));        //all.mem. for an array of structures
    for(i=0;i<N;i++){
        
        (*site)[i].neigh=(int *)malloc(N*sizeof(int));           // For each element of the array, all.mem. for each spin's neighbours
        (*site)[i].J=(int *)malloc(N*sizeof(int));              //For each element of the array, all. mem. for its couplings
        (*site)[i].degree=0;                                    //set the initial degree equal to zero
    }
    
    return;
    
}

/***************************************************/

// Generate a chain, 1D network, OPEN BOUNDARY CONDITIONS

void generate_open_chain(struct variable *site,int J1, int J11, int J2, int J22){
    
    int i;
    
    printf("Chain, OPEN boundary conditions, N = %d\n\n",N);
    
    site[0].neigh[0]=1; //left leaf
    site[0].degree=1;
    
    site[1].neigh[0]=0;
    
    // coupling of the left leaf and second node
    if(get_rand()<0.5){
        site[0].J[1]=J1; // J_{0,1}
        site[1].J[0]=J11;  // J_{1,0}
    }
    else{
        site[0].J[1]=J2;  // J_{0,1}
        site[1].J[0]=J22;  // J_{1,0}
    }
    
    // all other nodes
    for(i=1;i<N-1;i++){
        
        site[i].neigh[1]=i+1;  // to the right i ---> i+1
        site[i+1].neigh[0]=i;  // to the left i+1 ---> i
    
        // Couplings J_{ij} //
        if(get_rand()<0.5){
            site[i].J[i+1]=J1;  // J_{i,i+1}
            site[i+1].J[i]=J11;  // J_{i+1,i}
        }
        else{
            site[i].J[i+1]=J2;  // J_{i,i+1}
            site[i+1].J[i]=J22;  // J_{i+1,i}
        }
        
        // degree of any node 
        site[i].degree=2;
    }
    
    site[N-1].degree=1;
    
    return;
}

/***************************************************/

// Generate a chain, 1D network, PERIODIC BOUNDARY CONDITIONS

void generate_periodic_chain(struct variable *site, int J1, int J11, int J2, int J22){
    
    int i;
    
    printf("Chain, PERIODIC boundary conditions, N = %d\n\n",N);
    
    site[0].neigh[0]=N-1; //left leaf
    
    if(get_rand()<0.5){
        
        site[N-1].J[0]=J1;  // J_{i,i+1}
        site[0].J[N-1]=J11;   // J_{i+1,i}
          
    }
    else{
        site[N-1].J[0]=J2;  // J_{i,i+1}
        site[0].J[N-1]=J22;   // J_{i+1,i}
    }
   
    
    // all other nodes
    for(i=0;i<N-1;i++){
        
        site[i].neigh[1]=i+1;  // to the right i ---> i+1
        site[i+1].neigh[0]=i;  // to the left i+1 ---> i
        
        // Couplings J_{ij} //
        if(get_rand()<0.5){
            site[i].J[i+1]=J1;  // J_{i,i+1}
            site[i+1].J[i]=J11;  // J_{i+1,i}
        }
        else{
            site[i].J[i+1]=J2;  // J_{i,i+1}
            site[i+1].J[i]=J22;  // J_{i+1,i}
        }
        
        // degree of any node
        site[i].degree=2;
    }
    
    site[N-1].neigh[1]=0;
    site[N-1].degree=2;
    
    return;
}

/***************************************************/

// Print graph on the file ERRgraph.dat to be plotted using graphics and on MC_graph.dat
// to be read by the C program MC_DMP.c which runs MC simulations based on the same graph

void print_graph_on_file(struct variable *site, char **argv){
    
    int i,j;
    
    FILE *fp_graph,*fp_MC;
    
    char filename[101];
    
    sprintf(filename,"%s_ERR_graph.dat",argv[2]);
    fp_graph=fopen(filename,"w");
    sprintf(filename,"%s_MC_graph.dat",argv[2]);
    fp_MC=fopen(filename,"w");
    
    fprintf(fp_graph,"digraph {\n        graph [rankdir=LR];\n");
    for(i=0;i<N;i++){
        fprintf(fp_MC,"%d\n",site[i].degree);
        if(site[i].degree==0){
            fprintf(fp_graph,"        %d;\n",i);
        }
        for(j=0;j<site[i].degree;j++){
            if(site[i].J[site[i].neigh[j]]>0)
                fprintf(fp_graph,"        %d -> %d[color=\"red\"];\n",i,site[i].neigh[j]);
            if(site[i].J[site[i].neigh[j]]<0)
                fprintf(fp_graph,"        %d -> %d[color=\"blue\"];\n",i,site[i].neigh[j]);
                
            fprintf(fp_MC,"%d\t%d\t%d\t%d\n",i,site[i].neigh[j],site[i].J[site[i].neigh[j]],site[site[i].neigh[j]].J[i]); // print i --> i-neigh, J[i][nn]
        }
    }
    fprintf(fp_graph,"}\n");
    
    
    fclose(fp_graph);
    fclose(fp_MC);
    return;
}



/***************************************************/

double get_rand(void){
    
    return(1.0 * rand()/(1.0 + RAND_MAX));
}










