// This can be used in two ways: using the function print_links the adjacency matrix is read from the file of the networks
//and a file with a list of random links is written. This file is then used to change some link of the FERRO into SG.

// By using the function generate_rand_config a number of argv[2] files is created, each one contain a list of random number
// between 0 and N-1 with no repetition. These files are then used to generate the initial condition of the SG in order to have
//a number k of sample starting with the same initial conditions

// compile ./a.out Name_file_only_links to generate a random list of links
// compile ./a.out N Num_files to generate a number Num_files of files containing random list of
//numbers between 0 and N-1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#define directory "DMP_data" //directory containing the data


void print_links(char *argv[]);   // print links in a random order
double get_rand(void); // Rand numb generator

void generate_rand_config(char *argv[]); // generate a random list of number between 0 and N-1 with no repetitions
// used to generate random initial configuration of spins for the SG histogram program




/* =========================================================================== */
/*                                  MAIN                                       */
/* =========================================================================== */

int main(int argc, char *argv[]){
    
    struct variable *site;
    int i,j;
    
    chdir(directory); // Move to the directory
    
    srand(time(NULL)); //initialize get_rand(), random numb generator
    
    print_links(argv);
 //   generate_rand_config(argv);
    
    
    return 0;
}




/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */



/***************************************************/

// Print graph on the file ERRgraph.dat to be plotted using graphics and on MC_graph.dat
// to be read by the C program MC_DMP.c which runs MC simulations based on the same graph

void print_links(char *argv[]){
    
    int i, count, temp, randomIndex, ch, tot_change, l;
    int sitei,sitej;
    double rate,J1,J2;
    
    FILE *fp_MC, *fp_SG, *fp_link;
    char filename[101];
    
    fp_MC= fopen("link_couplings.dat","r");  // file from which I read the links. This file is generated with couplings.awk from the original file containing the structure of the graph, usually ..._MC_graph
    fp_link = fopen("list_links.dat","w"); // file where I write the random list of links which will be turned into negative
    
    if(fp_MC==NULL){
        fprintf(stderr,"PROBLEM OPENING FILE %s\n\n" ,"link_couplings.dat");
        exit(errno);
    }
    
    
    count=0; // count the total number of link in the file
    while(!feof(fp_MC))
    {
        ch = fgetc(fp_MC);
        if(ch == '\n')
        {
            count++;
        }
    }
    
    
    
    int nlink[3][count];  // Generate a random sequence of number with elements in [0,tot_link-1]
    double b,v;
    for (l = 0; l < count; l++){
        fscanf(fp_MC,"%d%d%lf%lf",&sitei,&sitej,&b,&v);
        nlink[0][l]= sitei;
        nlink[0][l]= sitej;
        nlink[2][l]=l;
        printf("%d\t%d\t%d\n",nlink[0][l],nlink[1][l],nlink[2][l]);
    }
    
    for (l = 0; l < count; l++) {    // shuffle array
        temp = nlink[2][l];
        randomIndex = rand() % count; // random number between 0 and N-1
        
        nlink[2][l]= nlink[2][randomIndex];  //shuffle the array to have a random ordering
        nlink[2][randomIndex] = temp;
    }
    
    for(i=0;i<count;i++) // print the sequence on a file
        fprintf(fp_link,"%d\t%d\t%d\n",nlink[0][l],nlink[1][l],nlink[2][l]);
    
    printf("\nTotal number of links = %d\nPrinted on the file list_links.dat\n\n",count);
    
    
    return;
}

/***************************************************/

void generate_rand_config(char *argv[]){
    
    int i,j,N,temp,randomIndex;
    FILE *fp_rand;
    char filename[101];
    
    N=atoi(argv[1]);
    int array[N];
    
    for(j=1;j<= atoi(argv[2]);j++){
        sprintf(filename,"rand_configuration_%d.dat",j);
        fp_rand=fopen(filename,"w");
    
        for (i = 0; i < N; i++) {     // fill array
            array[i] = i;       // create an ordered array of elements
        }
    
        for (i = 0; i < N; i++) {    // shuffle array
            temp = array[i];
            randomIndex = rand() % N; // random number between 0 and N-1
        
            array[i]= array[randomIndex];  //shuffle the array to have a random ordering
            array[randomIndex] = temp;
        }
    
        for(i = 0; i < N; i++)
            fprintf(fp_rand,"%d\n",array[i]);
        }
    
    return;
    
}









