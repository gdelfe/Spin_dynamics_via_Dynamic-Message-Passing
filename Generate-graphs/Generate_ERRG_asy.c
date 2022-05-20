//
//  Generate_ERRG.c
// To run the code on linux: $ gcc Generate_ERRG_comm_line.c -lm, then: $./a.out N c J1 J11 J2 J22, for instance,
// $./a.out 1000 3 1 1 -1 -1
//
//  Created by gino on 1/20/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define directory "DMP_data" //directory containing the data


struct variable{ //Structure for MC simulations
    
    int degree;   // number of generic neighbours
    int degree_i; // number of incoming neighbours of a given site
    int degree_o; //number of oucoming neighbours of a given site
    
    int *neigh;    // generic neighbours
    int *neigh_i; //incoming neighbours of a given site
    int *neigh_o;  // outcoming neighbours of a given site
    
    double *J;     // couplings J_{ij}
};

// get parameters from stdout ------
void get_parameter(int argc, char **argv);

void allocate_memory(struct variable **site); // allocate memory
void generate_ERRG(struct variable *site); // generate ERRG WITHOUT isolated nodes

void print_graph_on_file(struct variable *site, char **);   // print graph on 2 files
// Rand numb generator
double get_rand(void);

int N,c;  // N spin and average connectivity c
double J1,J11,J2,J22; // couplings


/* =========================================================================== */
/*                                  MAIN                                       */
/* =========================================================================== */

int main(int argc, char *argv[]){
    
    struct variable *site;
    int i,j;

    chdir(directory); // Move to the directory
    
    srand(time(NULL)); //initialize get_rand(), random numb generator
    
    get_parameter(argc,argv);
    allocate_memory(&site);
    generate_ERRG(site);
    print_graph_on_file(site,argv);
    

    return 0;
}




/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */


void get_parameter(int argc, char **argv){
    
    if(argc>1){
        
        N=atoi(argv[1]);
        c=atoi(argv[2]);
        J1=atof(argv[3]);  
        J2=atof(argv[4]);
        
    }

    return;

}

/***************************************************/

void allocate_memory(struct variable **site){
    
    int i;
    
    /* allocate memory for the array of structures */
    *site=(struct variable *)malloc(N*sizeof(struct variable));        //all.mem. for an array of structures
    for(i=0;i<N;i++){
        
        (*site)[i].neigh=(int *)malloc(N*sizeof(int));
        (*site)[i].neigh_i=(int *)malloc(N*sizeof(int));           // For each element of the array, all.mem. for each spin's neighbours_incoming
        (*site)[i].neigh_o=(int *)malloc(N*sizeof(int));           // For each element of the array, all.mem. for each spin's neighbours_outcoming
        (*site)[i].J=(double *)malloc(N*sizeof(double));              //For each element of the array, all. mem. for its couplings
        (*site)[i].degree=0;
        (*site)[i].degree_i=0;                                    //set the initial degree_incoming equal to zero
        (*site)[i].degree_o=0;                                    //set the initial degree_outcoming equal to zero
    }
    
    return;
    
}

/***************************************************/

void generate_ERRG(struct variable *site){
    
    int i,j,l;
    double z;
    
    
    printf("\nER random graph, N = %d, c = %d, ",N,c);
    
    z=((double)c/(double)N); //average connectivity of the graph
    
    for(i=0;i<N;i++){               //note: when calling this function, site[i].degree=0 initially
        for(j=0;j<N;j++){
            
            if(i!=j){ //esclude auto-connections
                if(get_rand()<z){  // get_rand() gives a random numb between [0,1]
                
                    if(site[j].J[i]==0.){  //check if j is already connect to i, if not then proceed, if yes move on. The network is fully asymmetric
                        if(get_rand()<0.5){
                            site[i].J[j]=(double)J1/c; // J_{ij}
                        }
                        else{
                            site[i].J[j]=(double)J2/c;  // J_{ij}
                        }
                        // the connection J[j][i]=0 is assigned when the network is printed on the file. See function print_graph_on_file
            
                        site[i].neigh_o[site[i].degree_o]=j; //connect i_out to j_in
                        site[i].degree_o++;
                        site[j].neigh_i[site[j].degree_i]=i; // connect j_in to i_out
                        site[j].degree_i++;
                    }
       
    
                    if((site[i].degree_i + site[i].degree_o) ==0){  // go back and recompute the nodes if it was isolated
                        i--;
                    }
                }
            }
        }
    }
    return;
}



/***************************************************/

// Print graph on the file ERRgraph.dat to be plotted using graphics and on MC_graph.dat
// to be read by the C program MC_DMP.c which runs MC simulations based on the same graph

void print_graph_on_file(struct variable *site, char **argv){
    
    int i,j,nn;
    
    FILE *fp_graph; // to be plotted graphically
    FILE *fp_MC;    // to be used by MC
    char filename[101];
    
    sprintf(filename,"%s_ERRgraph.dat",argv[5]);
    fp_graph=fopen(filename,"w");
    sprintf(filename,"%s_MC_graph.dat",argv[5]);
    fp_MC=fopen(filename,"w");
    
    fprintf(fp_graph,"digraph {\n        graph [rankdir=LR];\n");
    for(i=0;i<N;i++){
        fprintf(fp_MC,"%d\t%d\n",site[i].degree_i,site[i].degree_o);
        if(site[i].degree_i==0 && site[i].degree_o==0){ //if the node is isolated
            fprintf(fp_graph,"        %d;\n",i);
        }
        for(j=0;j<site[i].degree_o;j++){
            nn=site[i].neigh_o[j];
            
            if(site[i].J[nn]>0.)
                fprintf(fp_graph,"        %d -> %d[color=\"red\"];\n",i,site[i].neigh_o[j]);
            if(site[i].J[nn]<0.)
                fprintf(fp_graph,"        %d -> %d[color=\"blue\"];\n",i,site[i].neigh_o[j]);
            
            fprintf(fp_MC,"%d\t%d\t%lf\t0.\n",i,nn,site[i].J[nn]); // print i --> i-neigh, J[i][nn], J[nn][i]=0
        }
    }
    fprintf(fp_graph,"}\n");
    
    
    fclose(fp_graph);
    fclose(fp_MC);
    return;
}



/***************************************************/
/*
void generate_gaussian(){

    for(i=0;i<2)

}

/***************************************************/



double get_rand(void){
    
    return(1.0 * rand()/(1.0 + RAND_MAX));
}










