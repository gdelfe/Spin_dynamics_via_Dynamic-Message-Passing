//
// This code runs Monte Carlo simulations with glauber dynamics of a kinetic Ising model on a graph read from an imput file.
// Total and local magnetization m_tot and <s_i(t)>, correlations < s_i(t) s_j (t) >, < s_i(t) s_j (t+1) > and < s_i(t) s_j (t+2) >
// are computed from the initial time to the final time, therefore not only at equilibrium but also during the transient dynamics.

// ./a.out N M0 BETA Tfinal Teq Nsample WSAMPLE S NS NNS
// Gino Del Ferraro, October 2015, Stockholm.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#define directory "DMP_data" //directory containing the data

struct variable{
    int degree; // number of neighbours of a given site
    int spin_t0; // value of the spin at time s
    int spin_t1; // value of the spin at time s+1
    int spin_t2; //value of the spin at time s+2
    int *neigh; //neighbours of a given site
    double *J;     //couplings between i and all its neighbours
};

void convert_parameter(int argc, char **argv);

void allocate_memory(struct variable **site, double **m_eq, double ***m_loc, double ****C, double ****C01, double ****C02);
void read_ERRG(struct variable *site, char **argv);

void initialize_ERRG(struct variable *site);
int MC_Glauber(struct variable * site, double beta);
void swap_spins(struct variable *site);

double magn(struct variable *site);
void correlations_t0(struct variable *site, double **m_loc, double ***C, double ***C01, double ***C02);
void correlations_t1(struct variable *site, double **m_loc, double ***C, double ***C01, double ***C02);
void correlations_t2(struct variable *site, double *m_eq, double **m_loc, double ***C, double ***C01, double ***C02, int t);

void print_on_file(double *m_eq, double **m_loc, double ***C, double ***C01, double ***C02, char **argv);
double get_rand(void);

int N,c,Nsample,WSAMPLE,Tfinal,Teq,S,NS,NNS;
double BETA,M0;

/* ================================================ */
/*                    MAIN                         */
/* ================================================ */

int main(int argc, char *argv[]){
    
    chdir(directory); // Move to the directory
    
    int i,j,t;
    double m;
    struct variable *site;
    double *m_eq, **m_loc, ***C, ***C01, ***C02;
    FILE *fp_sample, *fp_av;

    srand(time(NULL)); //initialize get_rand(), random numb generator
    convert_parameter(argc,argv);
    allocate_memory(&site,&m_eq,&m_loc,&C,&C01,&C02);
    read_ERRG(site,argv);
    

    //--- Cicle over the samples ------ //
    for(i=0;i<Nsample;++i){
        initialize_ERRG(site); //initialize all the spins
        
        //--- compute correlations at t = 0 ---- //
        correlations_t0(site,m_loc,C,C01,C02);
        
        // MC dynamics for t = 1  //
        MC_Glauber(site,BETA);
        correlations_t1(site,m_loc,C,C01,C02);
        swap_spins(site);

       // ---- MC time dynamics for t>1 ----- //
        for(t=2;t<Tfinal;t++){
        
            MC_Glauber(site,BETA);
            correlations_t2(site,m_eq,m_loc,C,C01,C02,t);
            swap_spins(site); // spin(s+1) --> spin(s), spin(s+2) --> spin(s+1), where s is time.
           
        }
  //      if(i%WSAMPLE==0)   printf("sample = %d\n",i);
        
    }
    print_on_file(m_eq,m_loc,C,C01,C02,argv);

    printf("\nN = %d, m0 = %lf, beta = %lf, Tfinal = %d, Nsample = %d, WSAMPLE = %d, file = %s\n\n",N,M0,BETA,Tfinal,Nsample,WSAMPLE,argv[11]);
    
    free(m_loc);
    free(C);
    free(C01);
    free(C02);
    
    return 0;
}


    

/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */

void convert_parameter(int argc, char **argv){
    
    if(argc>1){
        
        N=atoi(argv[1]);
        M0=atof(argv[2]);
        BETA=atof(argv[3]);
        Tfinal=atoi(argv[4]);
        Teq=atoi(argv[5]);
        Nsample=atoi(argv[6]);
        WSAMPLE=atoi(argv[7]);
        S=atoi(argv[8]);
        NS=atoi(argv[9]);
        NNS=atoi(argv[10]);
    }
    
    return;
    
}

/***************************************************/

void allocate_memory(struct variable **site, double **m_eq, double ***m_loc, double ****C, double ****C01, double ****C02){
    
    int i,j;
    
 
    
    /* allocate memory for the array of structures */
    *site=(struct variable *)calloc(N,sizeof(struct variable));        //all.mem. for an array of structures
    
    /* allocate memory for the local magnetization for each spin */
    *m_loc =(double **)malloc(N*sizeof(double *)); // vector for the local magnetization 
    *m_eq =(double *)malloc(N*sizeof(double)); // vector for the equilibrium local magnetization
    
    /* allocate memory for the correlation at same time, different times */
    *C = (double ***)malloc(N*sizeof(double **));
    *C01 = (double ***)malloc(N*sizeof(double **));
    *C02 = (double ***)malloc(N*sizeof(double **));
    
    for(i=0;i<N;i++){
        
        (*site)[i].neigh=(int *)malloc(N*sizeof(int));           // For each element of the array, all.mem. for each spin's neighbours
        (*site)[i].J=(double *)malloc(N*sizeof(double));              //For each element of the array, all. mem. for its couplings
        (*site)[i].degree=0;                                    //set the initial degree equal to zero

        (*C)[i]= (double **)malloc(N*sizeof(double *));
        (*C01)[i]= (double **)malloc(N*sizeof(double *));
        (*C02)[i]= (double **)malloc(N*sizeof(double *));
        
        (*m_loc)[i]= (double *)calloc(Tfinal,sizeof(double));
        
        for(j=0;j<N;j++){
            
            (*C)[i][j]=(double *)calloc(Tfinal,sizeof(double)); // for each element of the matrix C_ij allocate a vector of dimension equal to the total time of dynamics
            (*C01)[i][j]=(double *)calloc(Tfinal,sizeof(double));
            (*C02)[i][j]=(double *)calloc(Tfinal,sizeof(double));
        }
    }
    
    return ;

}


/***************************************************/

void read_ERRG(struct variable *site, char **argv){
    
    int i,j,nn;
    double J1,J2;
    FILE *fp_file;
    char filename[101];
    
    sprintf(filename,"%s_MC_graph.dat",argv[11]);
    fp_file=fopen(filename,"r");
    
    if(fp_file==NULL){
        fprintf(stderr,"PROBLEM OPENING FILE %s\n\n" ,"MC_graph.dat");
        exit(errno);
    }
    
    for(i=0;i<N;i++){
        fscanf(fp_file,"%d",&site[i].degree);
       // printf("degree %d = %d\n",i,site[i].degree);
        for(j=0;j<site[i].degree;j++){
            fscanf(fp_file,"%d%d%lf%lf",&i,&nn,&J1,&J2);
            site[i].neigh[j]=nn;
            site[i].J[nn]=J1;
            site[nn].J[i]=J2;
      //        printf("%d ---> %d \t J[%d][%d] = %lf \t J[%d][%d] = %lf\n",i,nn,i,nn,site[i].J[nn],nn,i,site[nn].J[i]);
        }
        
    }
    
    return;
}

/***************************************************/

//This function initialize the spin to a given Magnetization M0 read from
//stdout. All the spins are set to 1 and then a portion of them is changed to -1
//in order to have the total M0. This is done by creating an array, shuffling it
// and taking the random index of the array: random number without repetition.

//Current time t=1

void initialize_ERRG(struct variable *site){
 
    int i,temp,randomIndex,rand_site;
    int *array;
    double M_in;
    
   // M_in= 0.5*(M0+1.);
   // printf("M_in = %lf (int) %d (1-M_in) = %lf (1-M_in)*N = %lf \n\n",M_in,(int)((1.-M_in)*N),(1.-M_in),round((1.-M_in)*N));
    
    FILE *fp_conf;
    fp_conf=fopen("initial_configuration.dat","r");
    
    
    //write initial magnetization on stdout
    double m0=0.;
    for(i=0;i<N;i++){
        fscanf(fp_conf,"%d",&site[i].spin_t1);
        m0+= site[i].spin_t1;
        site[i].spin_t0=0;
        site[i].spin_t2=0;
    }
        
    
    //-------- Random initialization that guarantees global m0 //
 /*
  
  array=(int *)malloc(N*sizeof(int));
  
  for (i = 0; i < N; i++) {     // fill array
  array[i] = i;       // create an ordered array of elements
  site[i].spin_t1=1; // initialize all the spin to 1
  site[i].spin_t0=0;
  site[i].spin_t2=0;
  }
  
  
    for (i = 0; i < N; i++) {    // shuffle array
        temp = array[i];
        randomIndex = rand() % N; // random number between 0 and N-1
        
        array[i]= array[randomIndex];  //shuffle the array to have a random ordering
        array[randomIndex] = temp;
    }
    
    //change to -1 only a fraction of the spins
    for(i=0;i<round((1-M_in)*N);i++){
        rand_site=array[i]; 
        site[rand_site].spin_t1=-1; //the site is picked from the shuffle array in order not to have repetition
    }
   // printf("\n\n");
    
 /*
  // ---------------------------------------------------------- //
    
    // write initial magnetization on stdout 
 /*   double m0=0.;
    for(i=0;i<N;i++){
        m0+= site[i].spin_t1;
    }
    printf(" m0 = %lf \n",m0/N); 
 */
    return;
}

/***************************************************/

// Current time t_1, flipped spins go in t_2

int MC_Glauber(struct variable *site, double beta){
    
    int i,j,nn;
    double pr,sum_neigh;
    
    // PARALLEL UPDATE RULE: all the spins updated at the same time //
    // One Metropolis sweep on all the N spins at the same time     //
    
    for(i=0;i<N;++i){ // cicle on all the i's
        sum_neigh=0.; //initialize the sum of the neighbours spins to zero
        
        //cicle on all the i-neighbours
        for(j=0;j<site[i].degree;j++){
            nn= site[i].neigh[j]; //nearest neighbourd site of site i
            sum_neigh += site[nn].J[i] * site[nn].spin_t1; // \sum_j J_{ij} s_j
        }
            
        // rate probability //
        pr= 0.5*(1-site[i].spin_t1*tanh(beta*sum_neigh));
        // Check wheather to accept or not //
        if(pr> get_rand()) //if GLAUBER condition satisfied
            site[i].spin_t2= -site[i].spin_t1; // spin flip
        
        else
            site[i].spin_t2= site[i].spin_t1; // NO spin flip
        
    
    }
    return 1;
}

/***************************************************/

void swap_spins(struct variable *site){
    
    int i;
    
    for(i=0;i<N;i++){
        site[i].spin_t0= site[i].spin_t1;  //First we swap 1 to 0
        site[i].spin_t1= site[i].spin_t2;  // and then 2 to 1 
        
    }
    return;
}


/***************************************************/

// Compute magnetization and correlation at the inizial time.
// OBS: Both need to be normalized by the Number of Sample when printed
// OBS: at this step the inizial time is time 1

void correlations_t0(struct variable *site, double **m_loc, double ***C, double ***C01, double ***C02){
    
    int i;
    
    for(i=0;i<N;i++){
        m_loc[i][0]+= site[i].spin_t1; // sum of local magnetization of spin i at time zero for different samples
    }
    
    //Correlation t,t
    C[S][NNS][0]+=site[S].spin_t1 * site[NNS].spin_t1; // Correlation C_{s,ns} at t-1 and t
    
    // Correlation t,t-1
    C01[S][NS][0]+=site[S].spin_t1 * site[NS].spin_t1; // Correlation C_{s,ns} at t-1 and t 
    C01[NS][S][0]+=site[NS].spin_t1 * site[S].spin_t1; // Correlation C_{s,ns} at t-1 and t
    C01[NS][NNS][0]+=site[NS].spin_t1 * site[NNS].spin_t1; // Correlation C_{s,ns} at t-1 and t
    C01[NNS][NS][0]+=site[NNS].spin_t1 * site[NS].spin_t1; // Correlation C_{s,ns} at t-1 and t
    
    C02[S][S][0]+= site[S].spin_t1 * site[S].spin_t1; //Correlation C_ij at t-2 and t
    C02[NS][NS][0]+= site[NS].spin_t1 * site[NS].spin_t1; //Correlation C_ij at t-2 and t
    C02[NNS][NNS][0]+= site[NNS].spin_t1 * site[NNS].spin_t1; //Correlation C_ij at t-2 and t
    
    
    return;
}

/***************************************************/

// Compute magnetization and correlations at any time 1 and between 0 and 1
// OBS: Both need to be normalized by the Number of Sample when printed
// and the current time is t=2

void correlations_t1(struct variable *site, double **m_loc, double ***C, double ***C01, double ***C02){

    int i;
    
    for(i=0;i<N;i++){
        m_loc[i][1]+= site[i].spin_t2; // sum of local magnetization of spin i at time zero for different samples
    }
    
    // Correlation t,t
    C[S][NNS][1]+=site[S].spin_t2 * site[NNS].spin_t2; // Correlation C_{s,ns} at t-1 and t
    
    // Correlation t,t-1
    C01[S][NS][1]+=site[S].spin_t2 * site[NS].spin_t1; // Correlation C_ij at t-1 and t
    C01[NS][S][1]+=site[NS].spin_t2 * site[S].spin_t1; // Correlation C_ij at t-1 and t
    C01[NS][NNS][1]+=site[NS].spin_t2 * site[NNS].spin_t1; // Correlation C_{s,ns} at t-1 and t
    C01[NNS][NS][1]+=site[NNS].spin_t2 * site[NS].spin_t1; // Correlation C_{s,ns} at t-1 and t
   
    // Auto-correlation t,t-2
    C02[S][S][1]+= site[S].spin_t2 * site[S].spin_t1; //Correlation C_ij at t-2 and t
    C02[NS][NS][1]+= site[NS].spin_t2 * site[NS].spin_t1; //Correlation C_ij at t-2 and t
    C02[NNS][NNS][1]+= site[NNS].spin_t2 * site[NNS].spin_t1; //Correlation C_ij at t-2 and t

    return;
}

/***************************************************/

// Compute magnetization and correlations at any time 1 and between 0 and 1
// OBS: Both need to be normalized by the Number of Sample when printed
// and the current time is t=2

void correlations_t2(struct variable *site, double *m_eq, double **m_loc, double ***C, double ***C01, double ***C02, int t){
    
    int i;
   
    for(i=0;i<N;i++){
        m_loc[i][t]+= site[i].spin_t2; // sum of local magnetization of spin i at time zero for different samples
       // if(t>Teq)
       //     m_eq[i] += site[i].spin_t2; // sum for the logal magnetization at equilibrium. Sum over the samples and over time
    }
    
    // Correlation t,t 
    C[S][NNS][t]+=site[S].spin_t2 * site[NNS].spin_t2; // Correlation C_{s,ns} at the same time t
    
    // Correlation t, t-1
    C01[S][NS][t]+=site[S].spin_t2 * site[NS].spin_t1; // Correlation C_ij at t-1 and t
    C01[NS][S][t]+=site[NS].spin_t2 * site[S].spin_t1; // Correlation C_ij at t-1 and t
    C01[NS][NNS][t]+=site[NS].spin_t2 * site[NNS].spin_t1; // Correlation C_ij at t-1 and t
    C01[NNS][NS][t]+=site[NNS].spin_t2 * site[NS].spin_t1; // Correlation C_ij at t-1 and t
    
    // Auto-correlation t,t-2
    C02[S][S][t]+= site[S].spin_t2 * site[S].spin_t0; //Correlation C_ij at t-2 and t
    C02[NS][NS][t]+= site[NS].spin_t2 * site[NS].spin_t0; //Correlation C_ij at t-2 and t  
    C02[NNS][NNS][t]+= site[NNS].spin_t2 * site[NNS].spin_t0; //Correlation C_ij at t-2 and t  
    
    
    return;
}

/***************************************************/

void print_on_file(double *m_eq, double **m_loc, double ***C, double ***C01, double ***C02, char **argv){

    int i,j,t,count;
    double m_tot;
    
    FILE *fp_meq, *fp_ml, *fp_mtot, *fp_C, *fp_C01, *fp_C02;
    char filename[101];
    
        sprintf(filename,"%s_global_magn.dat",argv[12]);
        fp_mtot=fopen(filename,"w");
    
        sprintf(filename,"%s_local_magn.dat",argv[12]);
        fp_ml=fopen(filename,"w");
    
        sprintf(filename,"%s_eq_magn.dat",argv[12]);
        fp_meq=fopen(filename,"w");

        sprintf(filename,"%s_corr.dat",argv[12]);
        fp_C=fopen(filename,"w");
    
        sprintf(filename,"%s_corr01.dat",argv[12]);
        fp_C01=fopen(filename,"w");
        
        sprintf(filename,"%s_corr02.dat",argv[12]);
        fp_C02=fopen(filename,"w");

    
        fprintf(fp_C,"# (1)time (2)C_{i,nnn}(t,t)disc (3)C_{i,nnn}(t,t)conn \n");
        fprintf(fp_C01,"# (1)time (2)C_{i,nn}(t,t-1)disc (3)C_{i,nn}(t,t-1)conn (4)C_{nn,i}(t,t-1)disc (5)C_{nn,i}(t,t-1)conn (6)C_{nn,nnn}(t,t-1)disc (7)C_{nn,nnn}(t,t-1)conn (8)C_{nnn,nn}(t,t-1)disc (9)C_{nnn,nn}(t,t-1)conn\n");
        fprintf(fp_C02,"# (1)time (2)C_{i,i}(t,t-2)disc (3)C_{i,i}(t,t-2)conn (4)C_{nn,nn}(t,t-2)disc (5)C_{nn,nn}(t,t-2)conn (6)C_{nnn,nnn}(t,t-2)disc (7)C_{nnn,nnn}(t,t-2)conn\n");
    
    
    
    // --------- t = 0 ------------ //
    m_tot=0. ;
    fprintf(fp_ml,"0 ");
    for(i=0;i<N;i++){
        
        m_tot+=m_loc[i][0];  // total magnetization at every give time t
        fprintf(fp_ml,"%lf ",m_loc[i][0]/((double)Nsample));
        fprintf(fp_meq,"%lf\t",m_eq[i]/(Nsample*(Tfinal-Teq)));
    }
    fprintf(fp_meq,"\n");
    fprintf(fp_ml,"\n");
    fprintf(fp_mtot,"0\t%lf\n",m_tot/((double) N*Nsample));
    
    fprintf(fp_C,"0\t%lf\t%lf",C[S][NNS][0]/Nsample, (C[S][NNS][0] - m_loc[S][0]*m_loc[NNS][0]/Nsample)/Nsample);
    
    fprintf(fp_C01,"0\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",C01[S][NS][0]/Nsample,(C01[S][NS][0]-m_loc[S][0]*m_loc[NS][0]/Nsample)/Nsample, C01[NS][S][0]/Nsample,(C01[NS][S][0]-m_loc[NS][0]*m_loc[S][0]/Nsample)/Nsample,C01[NS][NNS][0]/Nsample,(C01[NS][NNS][0]-m_loc[NS][0]*m_loc[NNS][0]/Nsample)/Nsample,C01[NNS][NS][0]/Nsample,(C01[NNS][NS][0]-m_loc[NNS][0]*m_loc[NS][0]/Nsample)/Nsample);
    
    fprintf(fp_C02,"0\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",C02[S][S][0]/Nsample,(C02[S][S][0] - m_loc[S][0]*m_loc[S][0]/Nsample)/Nsample,C02[NS][NS][0]/Nsample,(C02[NS][NS][0]-m_loc[NS][0]*m_loc[NS][0]/Nsample)/Nsample,C02[NNS][NNS][0]/Nsample,(C02[NNS][NNS][0]-m_loc[NNS][0]*m_loc[NNS][0]/Nsample)/Nsample);
    
    fprintf(fp_C,"\n");
    fprintf(fp_C01,"\n");
    fprintf(fp_C02,"\n");
    
     // --------- t = 1 ------------ //
    
    m_tot=0. ;
    fprintf(fp_ml,"1 ");
    for(i=0;i<N;i++){
        
        m_tot+=m_loc[i][1];  // total magnetization at every give time t
        fprintf(fp_ml,"%lf ",m_loc[i][1]/((double)Nsample));
    }
    fprintf(fp_ml,"\n");
    fprintf(fp_mtot,"1\t%lf\n",m_tot/((double) N*Nsample));
    
    fprintf(fp_C,"1\t%lf\t%lf",C[S][NNS][1]/Nsample, (C[S][NNS][1] - m_loc[S][1]*m_loc[NNS][1]/Nsample)/Nsample);
    
    fprintf(fp_C01,"1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",C01[S][NS][1]/Nsample,(C01[S][NS][1]-m_loc[S][1]*m_loc[NS][0]/Nsample)/Nsample, C01[NS][S][1]/Nsample,(C01[NS][S][1]-m_loc[NS][1]*m_loc[S][0]/Nsample)/Nsample,C01[NS][NNS][1]/Nsample,(C01[NS][NNS][1]-m_loc[NS][1]*m_loc[NNS][0]/Nsample)/Nsample,C01[NNS][NS][1]/Nsample,(C01[NNS][NS][1]-m_loc[NNS][1]*m_loc[NS][0]/Nsample)/Nsample);
    
    fprintf(fp_C02,"1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",C02[S][S][1]/Nsample,(C02[S][S][1] - m_loc[S][1]*m_loc[S][0]/Nsample)/Nsample,C02[NS][NS][1]/Nsample,(C02[NS][NS][1]-m_loc[NS][1]*m_loc[NS][0]/Nsample)/Nsample,C02[NNS][NNS][1]/Nsample,(C02[NNS][NNS][1]-m_loc[NNS][1]*m_loc[NNS][0]/Nsample)/Nsample);
    
    fprintf(fp_C,"\n");
    fprintf(fp_C01,"\n");
    fprintf(fp_C02,"\n");
    
    
    // --------- t >= 2 ------------ //
    for(t=2;t<Tfinal;t++){
        
        // magnetizations ---------------- //
        m_tot=0. ;
        fprintf(fp_ml,"%d ",t);
        for(i=0;i<N;i++){
            
            m_tot+=m_loc[i][t];  // total magnetization at every give time t
            fprintf(fp_ml,"%lf ",m_loc[i][t]/((double)Nsample));
        }
        
        fprintf(fp_ml,"\n");
        fprintf(fp_mtot,"%d\t%lf\n",t,m_tot/((double) N*Nsample));
        
        // correlations ------------------ //
        fprintf(fp_C,"%d\t%lf\t%lf",t,C[S][NNS][t]/Nsample, (C[S][NNS][t] - m_loc[S][t]*m_loc[NNS][t]/Nsample)/Nsample);
        
        fprintf(fp_C01,"%d \t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",t,C01[S][NS][t]/Nsample,(C01[S][NS][t]-m_loc[S][t]*m_loc[NS][t-1]/Nsample)/Nsample, C01[NS][S][t]/Nsample,(C01[NS][S][t]-m_loc[NS][t]*m_loc[S][t-1]/Nsample)/Nsample,C01[NS][NNS][t]/Nsample,(C01[NS][NNS][t]-m_loc[NS][t]*m_loc[NNS][t-1]/Nsample)/Nsample,C01[NNS][NS][t]/Nsample,(C01[NNS][NS][t]-m_loc[NNS][t]*m_loc[NS][t-1]/Nsample)/Nsample);
        
        
        fprintf(fp_C02,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",t,C02[S][S][t]/Nsample,(C02[S][S][t] - m_loc[S][t]*m_loc[S][t-2]/Nsample)/Nsample,C02[NS][NS][t]/Nsample,(C02[NS][NS][t]-m_loc[NS][t]*m_loc[NS][t-2]/Nsample)/Nsample,C02[NNS][NNS][t]/Nsample,(C02[NNS][NNS][t]-m_loc[NNS][t]*m_loc[NNS][t-2]/Nsample)/Nsample);
        
        fprintf(fp_C,"\n");
        fprintf(fp_C01,"\n");
        fprintf(fp_C02,"\n");
        
    }
 
    
    fclose(fp_C);  
    fclose(fp_C01);
    fclose(fp_C02);
    fclose(fp_ml);
    fclose(fp_mtot);
    
    return;
}

/***************************************************/

double get_rand(void){

    return(1.0 * rand()/(1.0 + RAND_MAX));
}







