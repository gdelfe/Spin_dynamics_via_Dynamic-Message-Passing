//
//  DMP_new.c
//  
//
//  Created by gino on 8/24/15.
//
//  ./a.out N M0 beta Tfinal S NS FileRead FilePrint

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#define directory "DMP_data" //directory containing the
#define OPEN 0
#define CLOSE 1

/* -------------------------------------------  */

struct variable{ //Structure for Network topology
        
    int degree; // number of neighbours of a given site
    int *neigh; //neighbours of a given site
    double *J;     // couplings J_{ij}
};
/* -------------------------------------------  */

struct dynamic_mp{ // Structure for dynamic message passing algorithm
    
    double **to; //value of the message from the site to another site
};

/* -------------------------------------------  */
void allocate_network(struct variable **site);
void allocate_DMP(struct variable *site, struct dynamic_mp **F, struct dynamic_mp **F0, struct dynamic_mp **T, double ***G, double ***G0, double ***w, double ***m_loc);
int read_ERRG(struct variable *site, char **argv);
/* ------------------------------------------------- */
void get_parameters(int argc, char **argv);
void open_close_files( FILE **fp_w, FILE **fp_T, FILE **fp_F, FILE **fp_G, FILE **fp_mi, FILE **fp_mtot, FILE **fp_corr, char **argv[], int var);
/* ------------------------------------------------- */
void initialize_F_G(struct variable *site, struct dynamic_mp *F0, double **G0);
void compute_F0_G0(int i, int degree, int cur_spin, int ind, double prod, struct dynamic_mp *F, double **G, struct variable *site);
void print_F_G(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F, double **G, FILE *fp_w, FILE *fp_T, FILE *fp_F, FILE *fp_G, int t);
void print_mathematica(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F, double **G, FILE *fp_w, FILE *fp_T, FILE *fp_F, FILE *fp_G, int t);
/* ------------------------------------------------- */
void get_w(struct variable *site, double **w);
/* ------------------------------------------------- */
void get_T(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F);
void normalize_T(struct variable *site, struct dynamic_mp *T);
/* ------------------------------------------------- */
void get_F(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F0, struct dynamic_mp *F);
/* ------------------------------------------------- */
void get_G(struct variable *site, double **w, struct dynamic_mp *T, double **G0 , double **G);
/* ------------------------------------------------- */
void magnetization(struct variable *site, double **G, FILE *fp_mi, FILE *fp_mtot, double **m_loc, int t);
/* ------------------------------------------------- */
void swap_F_G(struct variable *site, struct dynamic_mp *F0, struct dynamic_mp *F, double **G0, double **G);
/* ------------------------------------------------- */
void init_magn_corr(FILE *fp_corr,FILE *fp_mtot,FILE *fp_mi, double **m_loc);
void correlations(struct variable *site, double **w, double **G, FILE *fp_corr, double **m_loc, int t);
/* ------------------------------------------------- */
double get_rand(void);
/* ------------------------------------------------- */
void showbits(unsigned int x);

int N,Tfinal,S,NS;
double BETA,M0;


/* =========================================================================== */
/*                                  MAIN                                       */
/* =========================================================================== */

int main(int argc, char *argv[]){
    
    chdir(directory); // Move to the directory
    
    int i,l,t;
    double norm;
    struct variable *site;
    struct dynamic_mp *F, *F0, *T;
    double **w, **G, **G0, **m_loc;;
    
    printf("prova\n");
    
    FILE *fp_w,*fp_T,*fp_F, *fp_G, *fp_mi, *fp_mtot, *fp_corr;
    get_parameters(argc,argv); // get parameters from command line, stdin
    
    open_close_files(&fp_w,&fp_T,&fp_F,&fp_G,&fp_mi,&fp_mtot,&fp_corr,&argv,OPEN); // open files to write results
    
    // ------- ALLOCATE ----------- //
    allocate_network(&site); // Allocate memory for variable on the network
    read_ERRG(site,argv);   // read network topology from file
    allocate_DMP(site,&F,&F0,&T,&G,&G0,&w,&m_loc); //Allocate memory for Dynamic Message Passing algorithm
    
    

    // -------- INITIAL CONDITIONS AND FIXED VALUE ------- //
    get_w(site,w);
    initialize_F_G(site,F0,G0);
    
    printf(" ------- t = 0 ------- \n");
    init_magn_corr(fp_corr,fp_mtot,fp_mi,m_loc);  // compute and write initial correlation which will be Cij(0,0)
    
   // print_F_G(site,w,T,F0,G0,fp_w,fp_T,fp_F,fp_G,0);
  //  print_mathematica(site,w,T,F0,G0,fp_w,fp_T,fp_F,fp_G,0);
    
  
    // ----------- DYNAMICS -------------- //
    for(t=1;t<Tfinal;t++){
        
        printf(" ------- t = %d ------- \n",t);
 
        get_T(site,w,T,F0);
        normalize_T(site,T);

        get_F(site,w,T,F0,F);
     //   normalize_F(site,F);
    
        get_G(site,w,T,G0,G);
        magnetization(site,G,fp_mi,fp_mtot,m_loc,t);
        correlations(site,w,G,fp_corr,m_loc,t);
        
      //  print_F_G(site,w,T,F,G,fp_w,fp_T,fp_F,fp_G,t);
     //   print_mathematica(site,w,T,F,G,fp_w,fp_T,fp_F,fp_G,t);
        
        swap_F_G(site,F0,F,G0,G);
        
    }
    
    open_close_files(&fp_w,&fp_T,&fp_F,&fp_G,&fp_mi,&fp_mtot,&fp_corr,&argv,CLOSE); // close files to write results
    


    // printf("%d   %d     %d \n",1<<0,1<<1,1<<2);
 
    printf("\nN = %d, m0 = %lf, beta = %lf, Tfinal = %d, file = %s\n\n",N,M0,BETA,Tfinal,argv[7]);
    
    return 0;
    
}



/* =========================================================================== */
/*                          MODULE FUNCTIONS DEFINITIONS                       */
/* =========================================================================== */

// get parameter from stdout and convert them in the right format
void get_parameters(int argc, char **argv){
    
    if(argc>1){
        
        N=atoi(argv[1]);
        M0=atof(argv[2]);
        BETA=atof(argv[3]);
        Tfinal=atoi(argv[4]);
        S=atoi(argv[5]);
        NS=atoi(argv[6]);
    }
    
    return;
    
}

/*************************************************************************************/

void showbits(unsigned int x)
{
    int i;
    for(i=( sizeof(int)*2)-1; i>=0; i--)
        (x&(1<<i))?putchar('1'):putchar('0');
    
    printf("\n");
}

/*************************************************************************************/

// The function w depndends on the following variables w(si| {sj})

void get_w(struct variable *site, double **w){
    
    int i,si,cj,nj,nn,indw,degree;
    double sum;
    
    for(i=0;i<N;i++){
        degree = site[i].degree;
        for(si=0;si<2;si++){
        
            for(cj=0;cj<(int)pow(2,degree);cj++){ // for all the possible configuration
                sum=0.;
                for(nj=0; nj< degree; nj++){
                    nn = site[i].neigh[nj];
                    sum += site[nn].J[i]*(2*( (cj & (1<<nj)) != 0) - 1);
                    }
                indw = si*(1<<degree) | cj; //(si| 0 0 0) | (0| sj sj sj) = (si| sj sj sj)
                w[i][indw] = exp(BETA*(2*si-1)*sum)/(2*cosh(BETA*sum));
            }
        }
    }

    return;

}

/*************************************************************************************/

void initialize_F_G(struct variable *site, struct dynamic_mp *F0, double **G0){

    int i,j,l,nn;
    
    for(i=0;i<N;i++){
        compute_F0_G0(i,site[i].degree,0,0,1.,F0,G0,site);
    }
    
    return;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

void compute_F0_G0(int i, int degree, int cur_spin, int ind, double prod, struct dynamic_mp *F0, double **G0, struct variable *site){
    
    int k;
    if(cur_spin<degree){
        for(k=0;k<2;k++){
            compute_F0_G0(i,site[i].degree,cur_spin+1,2*ind+k,prod*0.5*(1+(2*k-1)*M0),F0,G0,site);
        }
    }
    else{
        int j,nn,s;
        for(j=0;j<site[i].degree;j++){
            nn=site[i].neigh[j];
            F0[i].to[nn][ind]=prod;
        }
        for(s=0;s<2;s++){
            G0[i][2*ind+s]=prod*0.5*(1+(2*s-1)*M0);
        }
    }
    return;
}

/*************************************************************************************/

// Initialize local magnetizations and correlations

void init_magn_corr(FILE *fp_corr,FILE *fp_mtot,FILE *fp_mi, double **m_loc){

    int si,sj,i;
    double corr=0;
    
    fprintf(fp_corr,"#(1)time (2)corr_S_NS disc (3) corr_S_NS conn \n");
    
    for(si=0;si<2;si++){
        for(sj=0;sj<2;sj++){
            corr += (2*si-1)*(2*sj-1)*0.5*(1+(2*si-1)*M0)*0.5*(1+(2*sj-1)*M0); // s_i s_j P(s_i)P(s_j)
        }
    }
    fprintf(fp_corr,"0\t%lf\t%lf\n",corr,corr-M0*M0); // corr(t,t-1) - m[S](t)m[NS](t-1)
    
    fprintf(fp_mi,"0\t");
    for(i=0;i<N;i++){
        fprintf(fp_mi,"%lf\t",M0); // print local magnetization
        
    }
    fprintf(fp_mi,"\n");
    fprintf(fp_mtot,"0\t%lf\n",M0); // print global magnetization
    
    m_loc[0][0] = M0;
    m_loc[1][0] = M0;

    return;

}


/*************************************************************************************/

/* Assume a structure (sk, sk, sk). Then w(si | sk sk sk sj) where sj can be placed anywhere among sk. 
 Whereas F(sip | sk sk sk). First of all we index: w(si | 0 0 0) with mask_si and F(sip | 0 0 0 ) with mask_ip. 
 Then with k, decimal for k = (sk sk sk). We create the rest for w the expression (sk sk sk) is modified inserting 
 sj in the position required (sk sk sj sk) for instance in pos 1. NOTE: the function works also for the leaves, 
 i.e. neighbours = 1. This is because we use a do-while instead that a while or a for loop.  */

void get_T(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F){
    
    int i, nj, si, ip, jp, nn, l;
    int k, indw, indw_j, indF;
    int  mask_jp, mask_ip, mask_si, degree;
    
    /* ---------------------- Cicles on the variable of the network i, nj = neigh of i ----------------------- */
    for(i=0;i<N;i++){ // cicle for all the sites
        degree = site[i].degree;
  //      printf("node i = %d\n",i);
        
        for(nj=0; nj<site[i].degree; nj++){  // nj is pos_jp, i.e neighbour of site i order by 0,1,2,3
            mask_jp = pow(2,nj)-1;    // 0000111 where #1 = nj
  //          printf("neigh_j numb[%d] = %d \n",nj,site[i].neigh[nj]);
            nn = site[i].neigh[nj];
            
            /* ------------------- Cicles on the variable dependencies of T(si|ip,jp) ------------------------ */
            l=0;
            for(si=0;si<2;si++){             
                mask_si = si*(1<<degree);   // either 001000 where 1 is in pos degree or it is 000000
                indw = mask_si; // (si | 0 0 0 0)
              
                for(ip=0;ip<2;ip++){
                    mask_ip = ip*(1<<(degree-1));  // either 00010 where 1 is in pos degree-1 or it is 00000
                    indF = mask_ip; // (sip | 0 0 0)
                
                    for(jp=0;jp<2;jp++){
                        T[i].to[nn][l]=0.;
    
                        /* -----------------------      CORE     ------------------------- */
                        /* ------    RHS of the equation = \sum_kp w(si| sk sk sj sk ) F(ip| sk sk sk)   ---------------- */
    
                        k=0;
      //                  printf("T[%d].to[%d][%d] = ",i,site[i].neigh[nj],l);
                        
                        do{         // k is a decimal representation for the binary k = (sk sk sk)
                                    // this is a cicle over all the k neighbours
        
                            indw_j = ((k & ~mask_jp)<<1) | (k & mask_jp); // create an new 0 position among the k spins like (sk sk 0 sk)
                            if(jp==1) indw_j |= (1<<nj); // changes the 0 bit to 1 only if jp=1 in the argument of the function, i.e. (sk sk 1 sk)
                                                        // Result of these two lines is (0 | sk sk sj sk)
        
                            indF = k | mask_ip;         // (0 | sk sk sk ) OR (sip | 0 0 0 ) = (sip | sk sk sk)
                            indw = indw_j | mask_si;    // (0 | sk sk sj sk) OR (si | 0 0 0 0) = (si | sk sk sj sk)
                            
                            T[i].to[nn][l] += F[i].to[nn][indF] * w[i][indw];
         //                   printf("F[%d]w[%d] + ",indF,indw);
        
                            k++;
                    
                        }while(k<(int)pow(2,degree-1));
                        l++;
        //                printf("\n");
                        /* -------------------------------------------------------------*/
                    }
                }
            }
        }
    }
    
    
    return;
    
}

/*************************************************************************************/

void get_F(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F0 , struct dynamic_mp *F){

    int i, k, si, ip, kp, nj, nk, pk, nn;
    int indw_j0, indw_j1, indw_T, indF0, indF1, indT;
    int mask_jp, mask_ip, mask_si, degree;
    double prod, F1;
    double wT;
    //double wT[100];
    
    /* ---------------------- Cicles on the variable of the network i, nj = neigh of i -----------------------*/
    for(i=0;i<N;i++){    // cicle for all the sites
        degree = site[i].degree;
      //  printf("\n\n\nnode i = %d  degree = %d \n\n",i,degree);
    
        for(nj=0; nj<degree; nj++){  // nj is pos_jp, i.e neighbour of site i order by 0,1,2,3
            mask_jp = pow(2,nj)-1; // 0000111 where #1 = nj
      //      printf("--------------- \n");
      //      printf("neigh_j[%d] = %d \n",nj,site[i].neigh[nj]);
            nn = site[i].neigh[nj];

            /* -----------------------------   LEFT HAND SIDE ---------------------------- */
            /* ---------- F_{i->nj}(si| sk sk sk) where variables are si and k = (sk sk sk) --------------*/
            for(si=0;si<2;si++){
                mask_si = si*(1<<degree); // either 00100 where 1 is in pos degree or it is 00000
            
                k=0;
                do{
                    indF1 = mask_si>>1 | k ; // (si | sk sk sk) -- shift to the right because F has one less neigh then w (mask_si is used for both F,w)
                    
                    F1=0.;
         //           printf("\nF[%d].to[%d][%d] = \n",i,nn,indF1);
                
                    /* -----------------------      RIGHT HAND SIDE     -------------------------*/
                    
                    for(ip=0;ip<2;ip++){   
                        mask_ip = ip*(1<<(degree-1)); // either 00010 where 1 is in pos degree or it is 00000
    
                        /* -----------------------      CORE     -------------------------*/
                        kp=0;
                        do{
                            indw_j0 = (((kp & ~mask_jp)<<1) | (kp & mask_jp)) | mask_si;  // (0 | sk sk 0 sk ) OR (si | 0 0 0 ) = (si | sk sk 0 sk)
                            indw_j1 = (indw_j0 | (1<<nj));      // (si | sk sk 0 sk ) OR (0 | 0 0 1 0 ) = (si | sk sk 1 sk)
                            indw_T = kp | mask_si;      // (0 | sk sk sk ) OR (si | 0 0 0 ) = (si | sk sk sk)
                
                       //      wT[indw_T] = w[i][indw_j0] + w[i][indw_j1];   // w(si| sk sk sk) = \sum_sj w(si | sk sk sj sk)
                             wT = 0.5*(w[i][indw_j0] + w[i][indw_j1]);   // w(si| sk sk sk) = \sum_sj w(si | sk sk sj sk)  --- 0.5 is to assure normalization of F --
                        
                            // ADD COMMENT ON THIS LOOP BECAUSE USE A SMART PROPERTY not easy to remember 
                            pk=0;  // position bit k-th in k = (sk sk sk)
                            prod=1.;
                            for(nk=0; nk<degree; nk++){
                                if(nk != nj){ // Esclude the neigh j from the computation of the \prod T[k \to i]
                                    indT = ip | ((kp &(1<<pk)) != 0)<<1 | ((k &(1<<pk)) != 0)<<2;
                                    prod *= T[site[i].neigh[nk]].to[i][indT]; // T[neigh_k].to[i][...]
                    //                printf("T[%d][%d]= %lf *",site[i].neigh[nk],indT,T[site[i].neigh[nk]].to[i][indT]);
                                    pk++;
                                }
                            }
                        
                            indF0 = kp | mask_ip;        // (0 | sk sk sk ) OR (ip | 0 0 0 ) = (ip | sk sk sk)
                        
                          //  F1 += prod * F0[i].to[nn][indF0] * wT[indw_T];
                            F1 += prod * F0[i].to[nn][indF0] * wT;
                //            printf("F0[%d]= %lf *0.5(w[%d]= %lf + w[%d]= %lf )  +  ",indF0,F0[i].to[nn][indF0],indw_j0,w[i][indw_j0],indw_j1,w[i][indw_j1]);
                            //   showbits(kp);
                
                            kp++;
                        
                        }while(kp<(int)pow(2,degree-1));
                        /* ------------------------------------------------------------- */
                    
                    } // cicle ip --- END RIGHT HAND SIDE
                    
                    F[i].to[nn][indF1] = F1;  // Assign value of F in the LHS
            //        printf("\n");
                    k++;
                
                }while(k<(int)pow(2,degree-1));
            } // cicle si
            
        } // cicle neigh_i (nj)
    } // cicle i


    return;

}

/*************************************************************************************/

void get_G(struct variable *site, double **w, struct dynamic_mp *T, double **G0 , double **G){
    
    int i, j, si, ip, jp, nj, nn, l;
    int indw, indG0, indG1, indT;
    int mask_ip, mask_si, degree;
    double prod, G_prod_T, G1, norm;

    
    /* ---------------------- Cicles on the variable of the network i, nj = neigh of i -----------------------*/
    for(i=0;i<N;i++){    // cicle for all the sites
        degree = site[i].degree;
     //   printf("\n\n\nnode i = %d  degree = %d \n\n",i,degree);
        
        
            /* -----------------------------   LEFT HAND SIDE ---------------------------- */
            /* ---------- G(si| sj sj sj) where variables are si and j = (sj sj sj) --------------*/
            for(si=0;si<2;si++){
                mask_si = si*(1<<degree); // either 00100 where 1 is in pos degree+1 or it is 00000
                
                j=0;
                do{     // cicle over the j variable on the LHS. One cicle, one value G1
                    indG1 = mask_si | j ; // (si | sj sj sj)
                    G1=0.;
                    
                //    printf("\nG1[%d][%d] = \n",i,indG1);
                    
                    /* -----------------------      RIGHT HAND SIDE     -------------------------*/
                    
                    jp=0;  // cicle on sj' (jp) on the RHS
                    do{
                        indw = mask_si | jp; // indw = (si| sj' sj' sj'), indeed jp = (sj' sj' sj') 
                            
                        // ------------ \sum_si' \prod_j T(sj | sj' si') G(si', {sj'}) ------------------- //
                        G_prod_T=0.;
                    //    printf("( ");
                        for(ip=0;ip<2;ip++){ // cicle on si'
                            mask_ip = ip*(1<<degree); // either 00010 where 1 is in pos degree+1 or it is 00000
                            indG0 = jp | mask_ip;        // (0 | sj sj sj ) OR (ip | 0 0 0 ) = (ip | sj sj sj)
                            
                            // Compute the product of T(sj | sj' si'). nj is the position of the j-th bit in j = (sj sj sj)
                            prod=1.;
                            for(nj=0; nj<degree; nj++){
                                indT = ip | ((jp &(1<<nj)) != 0)<<1 | ((j &(1<<nj)) != 0)<<2;
                                prod *= T[site[i].neigh[nj]].to[i][indT]; // T[neigh_j].to[i][...]
                          //      printf("T[%d].to[%d][%d]*",site[i].neigh[nj],i,indT);
                            }

                            G_prod_T += prod * G0[i][indG0];
                       //     printf(" G0[%d][%d] + ",i,indG0);
                        }
                    //    printf(")");
                            
                        G1 += G_prod_T * w[i][indw];
                     //   printf("*w[%d][%d] \n",i,indw);
                            
                        //   showbits(kp);
                            
                        jp++;
                            
                    }while(jp<(int)pow(2,degree));
                    /* ------------------------------------------------------------- */
                    
                    
                    G[i][indG1] = G1;  // Assign value of F in the LHS
                 //   printf("\n");
                    j++;
                    
                }while(j<(int)pow(2,degree));
            } // cicle si
            
    } // cicle i
    
    
    return;
    
}

/*************************************************************************************/

/* This function normalizes T: \sum_{s_i} T(s_i|s_i',s_j') = 1 */

void normalize_T(struct variable *site, struct dynamic_mp *T){
    
    int i,j,l,nn;
    double norm;
    
    for(i=0;i<N;i++){
        //       printf(" n_neigh function T = %d\n",site[i].degree);
        for(j=0;j<site[i].degree;j++){
            //            printf(" %d ---> %d\n",i,site[i].neigh[j]);
            nn=site[i].neigh[j];
            // norm=0.;
            /* Normalize T: \sum_{s_i} T(s_i|s_i',s_j') = 1 */
            for(l=0;l<4;l++){
                norm=T[i].to[nn][l]+T[i].to[nn][l+4];  // T(1|s_i',s_j')+T(-1|s_i',s_j')
                if(norm!=0){
                    T[i].to[nn][l]=T[i].to[nn][l]/norm;    // T(1|s_i',s_j')/norm
                    T[i].to[nn][l+4]=T[i].to[nn][l+4]/norm;  // T(-1|s_i',s_j')/norm
                }
            }
            
            
            //for(l=0;l<8;l++){
            //     printf("T[%d].to[%d][%d]= %lf \n",i,nn,l,T[i].to[nn][l]);
            //}
            
        }
        
    }
    
    return;
    
    
}

/*************************************************************************************/

void swap_F_G(struct variable *site, struct dynamic_mp *F0, struct dynamic_mp *F, double **G0, double **G){


    int i,j,l,nn,degree;
    
    for(i=0;i<N;i++){
        degree=site[i].degree;
        
        for(l=0;l<(int)pow(2,degree+1);l++){
            G0[i][l]=G[i][l];
        }
        
        for(j=0;j<degree;j++){
            nn = site[i].neigh[j];
            for(l=0;l<(int)pow(2,degree);l++){
                F0[i].to[nn][l]=F[i].to[nn][l];
            }
        }
    }


    return;

}

/*************************************************************************************/

void magnetization(struct variable *site, double **G, FILE *fp_mi, FILE *fp_mtot, double **m_loc, int t){
    
    int i,l,size;
    double P_up, P_down, m_i, m_tot;
    
    fprintf(fp_mi,"%d\t",t); 
    
    m_tot=0;
    for(i=0;i<N;i++){
        size=(int)pow(2,site[i].degree+1)/2;
        P_up=P_down=0.;
        
        for(l=0;l<size;l++){
            
            P_down+=G[i][l]; // \sum_{s_j} G(1,{s_j})
            P_up+=G[i][2*size-1-l]; // \sum_{s_j} G(-1,{s_j})
        }
        
        m_i = P_up - P_down;
        
        if(i==S) m_loc[0][1] = m_i;
        if(i==NS) m_loc[1][1] = m_i;
        
        fprintf(fp_mi,"%lf\t    ",m_i); // print local magnetization
        m_tot += m_i;
        
    } //cicle i
    fprintf(fp_mi,"\n");
    fprintf(fp_mtot,"%d\t%lf\n",t,(double)m_tot/N); // print global magnetization
    
}

/*************************************************************************************/

void correlations(struct variable *site, double **w, double **G, FILE *fp_corr, double **m_loc, int t){

    int si, sj, k, l, degree, nj, indw, check;
    int size, dim, mask_1, mask_0, mask_jp, mask_si;
    double P_ij[4], corr;
    
    degree = site[S].degree;
    size = 1<<degree; // 2^degree or (1 0 0 0) with one in position degree
    dim =  1<<(degree+1); // 2^(degree+1)
    
    double P_neigh[dim],P_i_n[dim];
        
    // --------- find out to which nj-neighbour NS correspond to. ----- //
    check=0;
    for(l=0;l<degree;l++){ 
        if(NS==site[S].neigh[l]){
            nj=l;
            check=1;
        }
    }
    if(check==0){
        printf("ERROR - Chosen neighbours %d doesn't correspond to one of the %d-neighbours\n",NS,S);
        exit(errno);
    }
    // ------------------------------------------------------ //
    
    // ------ marginalize G(si',{sj'}) to get P({sj'}) and Get P(si,{sj'}) = w(si,{sj'}) * P({sj'}) ----------- //
    mask_jp = pow(2,nj)-1; // 0000111 where #1 = nj
    
    for(l=0;l<size;l++){
        
        // marginalize G(si',{sj'}) to get P({sj'})
        P_neigh[l] = G[S][l] + G[S][size+l];
       // printf("P_neigh[%d] = G[%d][%d] + G[%d][%d]\n",l,S,l,S,size+l);
        

        P_i_n[l] = w[S][l] * P_neigh[l];  // P(si,{sj'}) with si = -1
      //  printf("P_i_n[%d] = w[%d][%d] * P_neigh[%d]\n",l,S,l,l);
        indw = size | l ;
        P_i_n[indw] = w[S][indw] * P_neigh[l]; // P(si,{sj'}) with si = 1
      //  printf("P_i_n[%d] = w[%d][%d] * P_neigh[%d]\n\n",indw,S,indw,l);
    }
    

    //--------------  marginalize P(si, {sj'}) to get P(si,sj') ----------------- //
    P_ij[0]=P_ij[1]=P_ij[2]=P_ij[3]=0.;
    for(si=0;si<2;si++){
        mask_si = si*(1<<degree);
        
        k=0;
        
     //   printf("\nP_ij[%d] = \n",2*si);
      //  printf("\nP_ij[%d] = \n",2*si+1);
        
        
        do{
        
            mask_0 = (((k & ~mask_jp)<<1) | (k & mask_jp)) | mask_si;  // (0 | sk sk 0 sk ) OR (si | 0 0 0 ) = (si | sk sk 0 sk)
            mask_1 = (mask_0 | (1<<nj));      // (si | sk sk 0 sk ) OR (0 | 0 0 1 0 ) = (si | sk sk 1 sk)
            
            P_ij[(2*si)] += P_i_n[mask_0];
      //      printf("P_in[%d] +  ",mask_0);
            P_ij[(2*si+1)] += P_i_n[mask_1];
         //   printf("P_in[%d] + ",mask_1);
            k++;
            
        }while(k<(int)pow(2,degree-1));
    
    }
   // printf("\n");
    
    //-- compute correlations a time t and t-1 between S and NS ------- //
    corr=0.;
    for(si=0;si<2;si++){
        for(sj=0;sj<2;sj++){
            
            corr += (2*si-1)*(2*sj-1)*P_ij[(2*si+sj)];
        }
    }
    fprintf(fp_corr,"%d\t%lf\t%lf\n",t,corr,corr-m_loc[0][1]*m_loc[1][0]); // corr(t,t-1) - m[S](t)m[NS](t-1)
  

    
    m_loc[0][0] = m_loc[0][1];  // m_loc[S](t-1) = m_loc[S](t);
    m_loc[1][0] = m_loc[1][1];  // m_loc[NS](t-1) = m_loc[NS](t);
    
    return;

}


    
/*************************************************************************************/

void open_close_files( FILE **fp_w, FILE **fp_T, FILE **fp_F, FILE **fp_G, FILE **fp_mi, FILE **fp_mtot, FILE **fp_corr, char **argv[], int var){
    
    char filename[101];
    
    // Open files with the name given from std input: useful if you want
    // to run more simulations at the same time --> write on different files
    if(var==0){
   
/*
        sprintf(filename,"%s_w_values.dat",(*argv)[8]);
        *fp_w=fopen(filename,"w");
        
        sprintf(filename,"%s_T_values.dat",(*argv)[8]);
        *fp_T=fopen(filename,"w");
         
        sprintf(filename,"%s_F_values.dat",(*argv)[8]);
        *fp_F=fopen(filename,"w");
         
        sprintf(filename,"%s_G_values.dat",(*argv)[8]);
        *fp_G=fopen(filename,"w");
         
  */       
         
        sprintf(filename,"%s_local_magn.dat",(*argv)[8]);
        *fp_mi=fopen(filename,"w");
        
        sprintf(filename,"%s_global_magn.dat",(*argv)[8]);
        *fp_mtot=fopen(filename,"w");
        
        sprintf(filename,"%s_corr01.dat",(*argv)[8]);
        *fp_corr=fopen(filename,"w");
        
        if(fp_corr==NULL){
            fprintf(stderr,"PROBLEM OPENING FILE %s\n\n" ,"XXX_MC_graph.dat");
            exit(errno);
        }
        
        
        
    }
    // Close files
    if(var==1){
    /*
        fclose(*fp_w);
        fclose(*fp_T);
        fclose(*fp_F);
        fclose(*fp_G);
      */ 
        
        fclose(*fp_mi);
        fclose(*fp_mtot);
        fclose(*fp_corr);

    }
    return;
    
}

/*************************************************************************************/

int read_ERRG(struct variable *site, char **argv){
    
    int i,j,nn, d_max;
    double J1,J2;
    FILE *fp_file;
    
    char filename[101];
    
    sprintf(filename,"%s_MC_graph.dat",argv[7]);
    fp_file=fopen(filename,"r");
    
    if(fp_file==NULL){
        fprintf(stderr,"PROBLEM OPENING FILE %s\n\n" ,"XXX_MC_graph.dat");
        exit(errno);
    }
    
    d_max=1;
    for(i=0;i<N;i++){
        fscanf(fp_file,"%d",&site[i].degree);
      //  printf("degree %d = %d\n",i,site[i].degree);
        for(j=0;j<site[i].degree;j++){
            fscanf(fp_file,"%d%d%lf%lf",&i,&nn,&J1,&J2);
            site[i].neigh[j]=nn;
            site[i].J[nn]=J1;
            site[nn].J[i]=J2;
    //        printf("%d ---> %d \t J[%d][%d] = %lf \t J[%d][%d] = %lf\n",i,nn,i,nn,site[i].J[nn],nn,i,site[nn].J[i]);
        }
        if(site[i].degree>d_max) // count the max degree of the network
            d_max=site[i].degree;
        
    }
    
    return (d_max);
}

/*************************************************************************************/

void allocate_DMP(struct variable *site, struct dynamic_mp **F, struct dynamic_mp **F0, struct dynamic_mp **T, double ***G, double ***G0, double ***w, double ***m_loc){
    
    int i,j;
    
    /* Matrix of dimension N*2^(d_i + 1) where d_i is degree of site i */
    /* to store the value of the function w(s_i | s_j, {s_k}) for any value of s   */
    *w=(double **)malloc(N*sizeof(double *));
    *G=(double **)malloc(N*sizeof(double *));
    *G0=(double **)malloc(N*sizeof(double *));
    
    *m_loc=(double **)malloc(2*sizeof(double *)); // local magn for spin S and NS
    (*m_loc)[0]=(double *)malloc(2*sizeof(double)); // local magn for spin S at time t and t-1
    (*m_loc)[1]=(double*)malloc(2*sizeof(double)); // local magn for spin S at time t and t-1
    
    
    *F=(struct dynamic_mp *)malloc(N*sizeof(struct dynamic_mp));
    *F0=(struct dynamic_mp *)malloc(N*sizeof(struct dynamic_mp)); // F at the previous time
    *T=(struct dynamic_mp *)malloc(N*sizeof(struct dynamic_mp));
        
    for(i=0;i<N;i++){
        
        (*w)[i]=(double *)malloc( ((int)pow(2,site[i].degree+1)) * sizeof(double)); //w[i][j], w_i(s_i|s_j,{s_k}_{\})
        (*G)[i]=(double *)malloc( ((int)pow(2,site[i].degree+1)) * sizeof(double)); // G[i][j], G_i(s_i,{s_j})
        (*G0)[i]=(double *)malloc( ((int)pow(2,site[i].degree+1)) * sizeof(double)); // G[i][j], G_i(s_i,{s_j}) at the previous time
        
        
        /* F[i].to[j] is F_{i \to ij }, F[i].to[j][l] cointains the value of */
        /* F_{i \to ij}(s_i,{s_k}_{\}) for a given realization of the spins  */
        (*F)[i].to=(double **)malloc( N * sizeof(double *)); // allego N spin perché l'indicizzazione che voglio non è ordinata 0,1,2 etc... ma può essere F[i].to[12], F[i].to[32], etc...
        (*F0)[i].to=(double **)malloc( N * sizeof(double *));
        
        (*T)[i].to=(double **)malloc( N * sizeof(double *));
        
        for(j=0;j<site[i].degree;j++){
            // nearest neighbour nn=site[i].neigh[j]
            /* All possible realization of F[i].to[j], which is F[i].to[j][l] */
            (*F)[i].to[site[i].neigh[j]]=(double *)malloc( ((int)pow(2,site[i].degree)) * sizeof(double)); // F_{i \to ij}(s_i,{s_k}_{\})
            (*F0)[i].to[site[i].neigh[j]]=(double *)malloc( ((int)pow(2,site[i].degree)) * sizeof(double));
            
            (*T)[i].to[site[i].neigh[j]]=(double *)malloc( (int)8 * sizeof(double)); // T[i].to[j][l] where l=8 values, 2^3 spin configurations
        }
    }
    
    return;
}


/*************************************************************************************/


void allocate_network(struct variable **site){
    
    int i;
    
    /* allocate memory for the array of structures */
    *site=(struct variable *)malloc(N*sizeof(struct variable));        //all.mem. for an array of structures
    for(i=0;i<N;i++){
        
        (*site)[i].neigh=(int *)malloc(N*sizeof(int));           // For each element of the array, all.mem. for each spin's neighbours
        (*site)[i].J=(double *)malloc(N*sizeof(double));              //For each element of the array, all. mem. for its couplings
        (*site)[i].degree=0;                                    //set the initial degree equal to zero
    }
    
    
    return;
    
}

/*************************************************************************************/

void print_F_G(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F, double **G, FILE *fp_w, FILE *fp_T, FILE *fp_F, FILE *fp_G, int t){
    
    int i,l,k,nn;
    
    fprintf(fp_w,"---------------------- t = %d -------------\n",t);
    fprintf(fp_T,"---------------------- t = %d -------------\n",t);
    fprintf(fp_G,"---------------------- t = %d -------------\n",t);
    fprintf(fp_F,"---------------------- t = %d -------------\n",t);
    
    for(i=0;i<N;i++){
        for(l=0;l<(int)pow(2,site[i].degree+1);l++){
            fprintf(fp_G,"G[%d][%d] = %lf\n",i,l,G[i][l]);
            fprintf(fp_w,"w[%d][%d] = %lf\n",i,l,w[i][l]);
            
        }
        fprintf(fp_G,"----------------------------\n");
        fprintf(fp_w,"----------------------------\n");
        
        for(k=0;k<site[i].degree;k++){
            for(l=0;l<(int)pow(2,site[i].degree);l++){
                nn=site[i].neigh[k];
                fprintf(fp_F,"F[%d]to[%d][%d] = %lf\n",i,nn,l,F[i].to[nn][l]);
            }
            fprintf(fp_F," - - - - - - - - - - - \n");
            
            for(l=0;l<8;l++){
                nn=site[i].neigh[k];
                fprintf(fp_T,"T[%d]to[%d][%d] = %lf\n",i,nn,l,T[i].to[nn][l]);
            }
            fprintf(fp_T," - - - - - - - - - - - \n");
            
            
        }
        fprintf(fp_F,"----------------------- \n");
        fprintf(fp_T,"----------------------- \n");
    }
    
    return;
}

/*************************************************************************************/

void print_mathematica(struct variable *site, double **w, struct dynamic_mp *T, struct dynamic_mp *F, double **G, FILE *fp_w, FILE *fp_T, FILE *fp_F, FILE *fp_G, int t){
    
    int i,j,l,nn,d_i;
    
    fprintf(fp_T,"*********************************** \n* \t\t t = %d \t\t\t*\n*********************************** \n",t);
    fprintf(fp_F,"*********************************** \n* \t\t t = %d \t\t\t*\n*********************************** \n",t);
    fprintf(fp_G,"*********************************** \n* \t\t t = %d \t\t\t*\n*********************************** \n",t);
    fprintf(fp_w,"*********************************** \n* \t\t t = %d \t\t\t*\n*********************************** \n",t);
    
    
    for(i=0;i<N;i++){
        
        d_i= site[i].degree;
        
        //    printf("i = %d N NEIGH = %d \n",i,site[i].degree);
        
        fprintf(fp_T,"i = %d \t N NEIGH = %d \n",i,site[i].degree);
        fprintf(fp_F,"i = %d \t N NEIGH = %d \n",i,site[i].degree);
        fprintf(fp_G,"i = %d \t N NEIGH = %d \n",i,site[i].degree);
        fprintf(fp_w,"i = %d \t N NEIGH = %d \n",i,site[i].degree);
        
        for(j=0;j<site[i].degree;j++){
            
            ///   printf(" %d --- > %d \n",i,site[i].neigh[j]);
            
            nn=site[i].neigh[j];
            fprintf(fp_T," %d ----> %d\n",i,nn);
            fprintf(fp_F," %d ----> %d\n",i,nn);
            
            
            for(l=0;l<8;l++){
                fprintf(fp_T,"T1[%d] = %lf;\n",l,T[i].to[nn][l]);
            }
            
            for(l=0;l<(int)pow(2,d_i);l++){
                fprintf(fp_F,"F0[%d] = %lf;\n",l,F[i].to[nn][l]);
            }
        }
        //     printf("------- \n");
        
        for(l=0;l<(int)pow(2,d_i+1);l++){
            fprintf(fp_w,"w[%d] = %lf;\n",l,w[i][l]);
            fprintf(fp_G,"G0[%d] = %lf;\n",l,G[i][l]);
        }
        
        fprintf(fp_G,"----------------\n");
        fprintf(fp_w,"----------------\n");
    }
    
    //  printf("----------------------\n----------------------\n----------------------\n");
    
    return;
    
}

/*************************************************************************************/

double get_rand(void){
    
    return(1.0 * rand()/(1.0 + RAND_MAX));
}








