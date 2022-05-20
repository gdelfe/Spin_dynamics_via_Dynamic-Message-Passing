//
//  generate_rand_fields.c
//  
//
//  Created by gino on 2/16/16.
//
//  Generate a list of N random field with bimodal distribution \pm h 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#define directory "DMP_data" //directory containing the data

double get_rand(void); // Rand numb generator
void generate_field(int N, double h); // print random binomial fields on a file


/* =========================================================================== */
/*                                  MAIN                                       */
/* =========================================================================== */

int main(int argc, char *argv[]){
    
    int N;
    double h;
    
    chdir(directory); // Move to the directory
    srand(time(NULL)); //initialize get_rand(), random numb generator
    N = atoi(argv[1]);
    h = atof(argv[2]);
    
    generate_field(N,h);
    
    
    return 0;
}

/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */

void generate_field(int N,double h){

    int i;
    FILE *fp_fields;
    
    fp_fields = fopen("rand_fields.dat","w");
    
    if(fp_fields == NULL){
        fprintf(stderr,"problem opening file\n");
        exit(errno);
    }
    
    for(i=0;i<N;i++){
        
        if(get_rand()<0.5)
            fprintf(fp_fields,"%.2lf\n",h);
        
        else
            fprintf(fp_fields,"%.2lf\n",-h);
        
        }
    
    fclose(fp_fields);
    
    return;

}


/***************************************************/


double get_rand(void){
    
    return(1.0 * rand()/(1.0 + RAND_MAX));
    
}








