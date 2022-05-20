/* ================================ */
/* == This code considers an Ising Chain and computes the value of the energy
 from the initial value to the equilibrium one. To achive this goal it first computes the value
 of the correlation R_ij at any time using the Glauber's formula, which is a differential
 equation formula for R_ij, and so the code uses a Runge-Kutta4  method to integrate this equation.
 Once R(t) is computed for anytime the code computes the Energy using E= - 1/N \sum_{i} R_{i,i+1}.
 The chain considered is closed, i.e. i=0 is equal to i=N */
/* ================================ */

/* difference with Ising_rk4.c : here The energy is computed inside the time cicle (k) and printed at any time */

#include<stdlib.h>
#include<stddef.h>
#include<stdio.h>
#include<math.h>


#define Temp 2.  // Temperature
#define dt 0.05 // Time step
#define m0 1. // Initial magnetization
#define NT 200 // Numb of time step iteration
#define N 100 // Matrix dimension, number of spins
#define OUT1 "GLAUBER_en_corr.dat" //file used to print the energy results

double ***create_matrix(void);
void initiate(double ***R);
void compute_kappa(double ***R, double ***kappa);
void compute_temp(double ***R, double ***kappa, double ***temp); // function of derivatives
void compute_energy(double ***R, double *E,double *corr, FILE *en_file, int k);
void free_matrix(double ***mat);

double g; //gamma constant, for Glauber = tanh(2*\beta*J)
int delta;  

int main(void)
{
    int i,j,t;
    double ***R, *E,*corr;        //correlation matrix - R=2xNxN where 2 is referred to vectors at time t, and t+1
    double ***kappa1,***kappa2,***kappa3,***kappa4;    // matrixes for kappa=2xNxN
    double ***temp1,***temp2,***temp3;     // matrixes for temp variables 2xNxN
    
    FILE *en_file;
    en_file = fopen(OUT1,"w");              // file for Energy data
    
    g=tanh(2./Temp);
    
    
    /* == Creation of Space for the Energy and Correlation Matrix - dynamic allocation == */
    /* first elements of the matrix is 0, last one is N-1 which coincide with 0 for the boundary conditions */
    R=create_matrix();
    kappa1=create_matrix();
    kappa2=create_matrix();
    kappa3=create_matrix();
    kappa4=create_matrix();
    temp1=create_matrix();
    temp2=create_matrix();
    temp3=create_matrix();
    E=(double *)malloc(N*sizeof(double));     /* Create space for the Energy vector */
    corr=(double *)malloc(N*sizeof(double));     /* Create space for the correlation vector */	
    /* ========================================================= */
    
    /* ===== INITIAL CONDITION ===== */
    initiate(R);
    E[0]=0.;
    corr[0]=0;
    compute_energy(R,E,corr,en_file,0); // Compute and print the initial value for the Energy
   
    /* ====== INTEGRATION FOR THE CORRELATION R_ij - RK fourth order ====== */
    for(t=1; t<NT; ++t){
        if(t % 10 == 0) printf("\n iteration numb: %d",t);
           E[t]=0.;
        
        delta=1;
        /* kappa1 --- compute_kappa(R,kappa,k) k1= h*f(R) */
        compute_kappa(R,kappa1);
        
        /* temp1 --- compute_temp(R,kappa,temp,k) t1= R+0.5*k1 */
        compute_temp(R,kappa1,temp1);
        
        /* kappa2 --- compute_kappa(R,kappa,k) k2=h*f(R+0.5*k1)=h*f(temp1) */
        compute_kappa(temp1,kappa2);
        
        /* temp2 --- compute_temp(R,kappa,temp,k) t2= R+0.5*k2 */
        compute_temp(R,kappa2,temp2);
        
        /* kappa3 --- compute_kappa(R,kappa,k) k3= h*f(R+0.5*k2) = h*f(temp2) */
        compute_kappa(temp2,kappa3);
        
        delta=2; /*with this command, in the temp computation I compute K3 and not 0.5*k3*/
        /* temp3 --- compute_temp(R,kappa,temp,k) t3= y1+k3 */
        compute_temp(R,kappa3,temp3);
        
        /* kappa4 --- compute_kappa(R,kappa,k) k4= h*f(temp3)*/
        compute_kappa(temp3,kappa4);
        
        /* ============= ITERATION y(t+1)= y(t)+1/6[k1(t)+ 2*k2(t) + 2*k3(t) + k4(t)] ============================ */
        for (i=0;i<N;++i){ // Iteration over the elements
            for (j=0;j<N;++j){
            
                R[1][i][j]=R[0][i][j]+1./6*(kappa1[0][i][j]+2*kappa2[0][i][j]+2*kappa3[0][i][j]+kappa4[0][i][j]);
               // R[1][j][i]=R[1][i][j]; 
            }
        }
        R[1][N-1][0]=R[1][0][N-1]=R[1][0][0]; // Fixed condition for the corners of the matrix: R_{n-1,0 }= R_{00} = R_{0,n-1}  }
        /* =========================================================*/
        
       for(i=0; i<N; ++i){
            for(j=0; j<N; ++j){
                R[0][i][j]=R[1][i][j]; //Exchange for the R time correlation values before to compute the Energy
            }
        }
        compute_energy(R,E,corr,en_file,t);
        }
    
    /* Print information about the code and the used value of the variables */
    printf("\n\n Ising chain with N=%d  spins. Number of iterations NT=%d \n" , N, NT);
    printf(" Temp=%g,  initial magnetization=%g, time step dt=%g\n\n",Temp,m0,dt);
    printf(" The Energy values are in the file: '%s', time vs E(t) \n\n",OUT1);
    
    /* Compute and print Energy on a file */
    printf("\n Energia stampata \n\n");
    
    
    /* ==== Return of the Space for the correlation matrix ====*/
    free_matrix(R);
    free_matrix(kappa1);
    free_matrix(kappa2);
    free_matrix(kappa3);
    free_matrix(kappa4);
    free_matrix(temp1);
    free_matrix(temp2);
    free_matrix(temp3);
    
    return 0;
}






/* ------------------------------------------------------------------------------------------------------  */
/* ===== FUNCTIONS DEFINITIONS * ====  */
/* ------------------------------------------------------------------------------------------------------  */

/* ===============================  */
/* Allocate space for the matrix    */
/* ===============================  */

double ***create_matrix(void)
{
    double ***mat;
    int i,j;
    
    mat=(double ***)malloc(2*sizeof(double**));            // dynamic allocation of memory for matrix R, at time t and t+1
    
    for(i=0; i<2; ++i){
        mat[i] = (double **)malloc(N*sizeof(double*)); // dynamic allocation for the N spins of the chain
        for(j=0; j<N; ++j){
            mat[i][j] = (double *)malloc(N*sizeof(double)); //dynamic allocation for the N spins of the chain
        }
    }
    
    return(mat);
}

/* ==========================  */
/* Initiate the matrix         */
/* ==========================  */

void initiate(double ***R){
    
    int i,j;
    
    for (i=0;i<N;i++){
        R[0][i][i]=1; // Autocorrelation for diagonal elements
        for (j=0;j<i;j++){
            R[0][i][j]=pow(m0,2); // Initial value=m0^2 for the other elements
            R[0][j][i]=R[0][i][j]; //Symmetric condition, the correlation matrix is symmetric
        }
    }
    
    for(j=0; j<N-1; ++j){
        R[0][N-1][j]=R[0][0][j];  // Boudary condition: first column = last column */
    }
    R[0][0][N-1]=R[0][N-1][0]; // Boundary condition for the corners of the matrix
    
    return;
}

/* ==========================  */
/* Kappa computation           */
/* ==========================  */
void compute_kappa(double ***R, double ***kappa){
    
    int i,j;
    int a,b,c,d,e,f;
    
    for (i=0;i<N;++i){ // Iteration over the elements
        for (j=0;j<N;++j){ // Iteration over the elements
            
            a=i; /* a,b,c,d,e,f: sets of parameters which takes account of the fact that the chain is closed - boundary conditions */
            b=j;
            c=j-1;
            d=j+1;
            e=i-1;
            f=i+1;
            
            if(i==0)
                e=N-2;
            if(i==N-1){
                a=0;
                f=1;}
            if(j==0)
                c=N-2;
            if(j==N-1){
                b=0;
                d=1;}
            
            // k1= h*f(x_n,y_n)
            kappa[0][i][j]=dt*(-2.*(R[0][a][b])+0.5*g*(R[0][a][c]+R[0][a][d]+R[0][e][b]+R[0][f][b])); /* Fromulas from Glauber's paper */
            //    kappa[k][j][i]=kappa[k][i][j]; //Symmetry condition for the correlation matrix
        }
        kappa[0][i][i]=0.0; //dyagonal elements are not changed: R[i][i]=1 at any time
    }
    kappa[0][N-1][0]=kappa[0][0][N-1]=kappa[0][0][0]; // Fixed condition for the corners of the matrix: R_{n-1,0 }= R_{00} = R_{0,n-1}  }
    
    return;
}

/* =========================================  */
/* Computation of the temporary variable TEMP */
/* =========================================  */
/* here delta is always = 1 except when one compute k3, in that case delta=2 */

void compute_temp(double ***R, double ***kappa, double ***temp){
    
    int i,j;
    
    for (i=0;i<N;++i){ 
        for (j=0;j<N;++j){
            
            temp[0][i][j]= R[0][i][j]+delta*0.5*(kappa[0][i][j]);              /* temp_n(i) = y_n(i)+1/2*k1 */
            //   temp[k][j][i]=temp[k][i][j];        //Symmetric statement
            }
    }
    
    return;
}

/* =========================================  */
/* Energy Computation and Printing            */
/* =========================================  */

void compute_energy(double ***R, double *E, double *corr, FILE *en_file, int t){
    
    int i,j;
    
    /* Compute the average energy at a given time k and the average correlation C(t,t+2) */
    // averaged over all the spins // 
    for (i=0;i<N-2;i++){
        E[t]=E[t]+R[0][i][i+1];         // nearest neighbour correlation 
        corr[t]=corr[t]+R[0][i][i+2];  // next-nearest neighbours correlation
    }
    E[t]=E[t]+R[0][N-2][N-1];
    
    E[t]= - E[t]/((double)(N-1)); //Normalization for the Energy
    corr[t]=corr[t]/((double)(N-1));
    fprintf(en_file,"%d\t%lf\t%lf\t%lf\n", t, E[t],corr[t],R[0][0][2]); // print t vs E(t), Corr_{i,i+2}(t) AVERAGED, corr_{0,2}(t)
    
    return;
}

/* ==========================  */
/* DEALLOCATE SPACE        ==== */
/* ==========================  */

void free_matrix(double ***mat){
    
    int i,j;
    
    for(i=0; i<2; ++i){
        for(j=0; j<N; ++j){
            free(mat[i][j]);
        }
        free(mat[i]);
    }
    free(mat);
    
    return;
}




