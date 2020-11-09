//
//  main.cpp
//  incubation_continuous
//
//  Created by Piet Groeneboom on 07/11/2020.
//

//
//  Simulations of 1000 samples of n=1000 for
//  the distribution of the incubation time
//
//  Created by Piet Groeneboom on Novemeber 5, 2020.
//  Copyright (c) 2020 Piet Groeneboom. All rights reserved.


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <random>
#include <string.h>

using namespace std;

#define SQR(x) ((x)*(x))

typedef struct
{
    double x;
    int index;
}
data0_object;


typedef struct
{
    double x;
    double y;
}
data_object;

typedef struct
{
    int i;
    int j;
}
index_object;

int     n,n1,ngrid,**N;
double  h,*loc,*tt,*tt_bootstrap,max_obs1,min_obs2;
double  alpha,*F,*F_bootstrap,*MLE,*MLE_bootstrap,*SMLE_bootstrap,*dens_bootstrap;
int     *index1,*index2,*index0,*ind,*ind1,*sec;
double  *data0,*data1,*data2,*data3,*data1_bootstrap,*data2_bootstrap;
double  a,b,*pp,*pp_bootstrap,*grid,*SMLE,*SMLE1,*dens1,*dens,*qq;
double  *cumw,*cs,*yy,*yy_new,*grad,*w;
double  inprod,partsum,M,M1;

int     compare (const void * a, const void * b);
int     compare2(const void *a, const void *b);

double  KK(double x);
double  K(double x);
double  dens_estimate(double A, double B,  int njumps, double jumploc[], double p[], double h, double u);
double  bdf(double A, double B,  int m, double tt[], double p[], double h, double u);
double  dens_estimate2(int n, double A, double B1, double B2, double data1[], double data2[], double u, double v, double h1, double h2);

double  criterion(int n1, int **M, double yy[]);
double  f_alpha(int n1, int **N, double alpha);
double  f_alpha_prime(int n1, int **N, double alpha);
int     fenchelviol(int n1, double b[], double nabla[], double tol, double *inprod, double *partsum);

void    isoreg(int n1, int **N, double yy[], double F[], double cumw[], double grad[], double w[]);
void    cumsum(int n1, double yy[], double cumw[], double grad[], double w[]);
void    weights(int n1, int **N, double yy[], double w[]);
void    gradient(int n1, int **N, double yy[], double grad[]);
void    transfer(int first, int last, double a[], double b[]);
double  golden(int n1, int **N, double (*f)(int,int**,double));

int     compute_mle(int n, double data1[], double data2[], double F[], double tt[], double pp[]);
double  MSE_dens(int n, int B, double h1, double h2);
double  MSE_dens2(int ngrid, double B, double h);
int     compare (const void * a, const void * b);
void    data_incub(int n, double a, double b, double data_Weibull[], int seed);

double Weibull_df(double a, double b, double x);
double Weibull_inv(double a, double b, double u);
double Weibull_dens(double a, double b, double x);
void   data_exp(int n, double data1[], double data2[], double data3[], int seed);
void   sort_data(int n, double data1[], double data2[]);
void   sort_data0(int m, double data0[], int index0[]);
void   sort_index(int n, int index1[]);


int main(int argc, const char * argv[])
{
    int i,iter,NumIt,seed;
    double  a1;
    FILE    *outputfile = NULL;
    
    n=1000;
    M=30;
    M1=20;
    
    NumIt= 1000;
    
    cout << endl;
    
    a=3.035140901;
    b=0.002619475;
    
    ngrid=10;

    grid = new double[ngrid+1];
    dens = new double[ngrid+1];
    
    N = new int *[2*n+2];
    for (i=0;i<2*n+2;i++)
        N[i] = new int[3];
    
    for (i=0;i<=ngrid;i++)
        grid[i] = 2.0+i*1.0;

    data0 = new double[2*n];
    data1 = new double[n];
    data2 = new double[n];
    data3 = new double[n];
    
    index0 = new int[2*n];
    index1 = new int[2*n+1];
    index2 = new int[2*n+1];
    ind    = new int[2*n+1];
    ind1    = new int[2*n+1];
    sec  = new int[2*n+1];
    
    for (i=0;i<2*n;i++)
        index1[i]=i;
    
    F =  new double[n+2];
    
    MLE =  new double[ngrid+1];

    pp = new double[2*n+2];
    tt = new double[2*n+2];
     
    cumw = new double[2*n+1];
    cs = new double[2*n+1];
    yy = new double[2*n+2];
    yy_new = new double[2*n+2];
    grad = new double[2*n+1];
    w = new double[2*n+1];
    
    h=3.4;
    
    cout << "bandwidth = " << 3.4 << endl << endl;
    
    outputfile=fopen("dens_values.txt","w");
    rewind(outputfile);
    
    for (iter=0;iter<NumIt;iter++)
    {
        cout << iter+1 << endl;
 
        seed = rand();
        data_exp(n,data1,data2,data3,seed);
        n1 = compute_mle(n,data1,data2,F,tt,pp);
        for (i=0;i<=ngrid;i++)
        {
            a1=dens_estimate(0.0,50.0,n1+1,tt,pp,grid[i],h);
            fprintf(outputfile,"%15.10f",a1);
        }
        
        fprintf(outputfile,"\n");
    }
    
    
    fclose(outputfile);
    
    printf("\n\n");
    
    outputfile=fopen("MLE.txt","w");
    rewind(outputfile);
       
    for (i=0;i<=n1+1;i++)
        fprintf(outputfile,"%15.10f     %15.10f\n",tt[i],F[i]);
     
    fclose(outputfile);
    
    outputfile=fopen("data_exp.txt","w");
    rewind(outputfile);
    
    for (i=0;i<n;i++)
        fprintf(outputfile,"%15.10f   %15.10f\n",data1[i],data2[i]);
        
    fclose(outputfile);
    
    
    outputfile=fopen("ICM_masses.txt","w");
    rewind(outputfile);
    
    for (i=1;i<=n1+1;i++)
        fprintf(outputfile,"%15.10f     %15.10f\n",tt[i],pp[i]);
        
    fclose(outputfile);
    
    // free memory
    
    for (i=0;i<2*n+2;i++)
        delete[] N[i];
    delete[] N;
    
    delete[] dens; delete[] grid; delete[] F;  delete[] MLE;
    delete[] cumw; delete[] cs; delete[] yy; delete[] yy_new; delete[] grad; delete[] w;
    delete[] tt; delete[] pp; delete[] data0; delete[] data1; delete[] data2; delete[] data3;
    delete[] index0; delete[] index1; delete[] index2; delete[] ind; delete[] ind1;
    delete[] sec;
    
    return 0;
}


double golden(int n1, int **N, double (*f)(int,int**,double))
{
    double a,b,eps=1.0e-10;
    
    a=0;
    b=1;
    
    double k = (sqrt(5.0) - 1.0) / 2;
    double xL = b - k*(b - a);
    double xR = a + k*(b - a);
    
    while (b-a>eps)
    {
        if ((*f)(n1,N,xL)<(*f)(n1,N,xR))
        {
            b = xR;
            xR = xL;
            xL = b - k*(b - a);
        }
        else
        {
            a = xL;
            xL = xR;
            xR = a + k * (b - a);
        }
    }
    return (a+b)/2;
    
}

double f_alpha(int n1, int **N, double alpha)
{
    int i;
    
    for (i=1;i<=n1;i++)
        yy_new[i]=(1-alpha)*F[i]+alpha*yy[i];

    return criterion(n1,N,yy_new);
}

double f_alpha_prime(int n1, int **N, double alpha)
{
    int        i;
    double    sum;
    
    sum=0;
            
    for (i=1;i<=n1;i++)
            sum -= (yy[i]-F[i])*grad[i]/((1-alpha)*F[i]+alpha*yy[i]);
    
    return sum;
}


int fenchelviol(int n1, double yy[], double grad[], double tol, double *inprod, double *partsum)
{
    double    sum,sum2;
    int    i;
    int    fenchelvioltemp;
    
    fenchelvioltemp = 0;
    
    sum=0;
    sum2=0;
    
    for (i=1;i<=n1;i++)
    {
        sum += grad[i];
        if (sum < sum2)
            sum2 = sum;
    }
    
    sum=0;
    for (i=1;i<=n1;i++)
        sum += grad[i]*yy[i];
    
    *inprod = sum;
    *partsum = sum2;
    
    if ((fabs(sum) > tol) || (sum2 < -tol) ) fenchelvioltemp = 1;
    
    return fenchelvioltemp;
}

void transfer(int first, int last, double a[], double b[])
{
    int    i;
    for (i = first; i<= last;i++)    b[i] = a[i];
}

double criterion(int n1, int **N, double yy[])
{
    int i;
    double sum=0;
    
    for (i=1;i<=n1;i++)
    {
        if (N[i][0]>0)
            sum -= N[i][0]*log(yy[i]);
        
        if (N[i][2]>0)
            sum -= N[i][2]*log(1-yy[i]);
        if (N[i][1]>0)
        {
            if (sec[i]>i)
                sum -= N[i][1]*log(yy[sec[i]]-yy[i])/2;
            else
                sum -= N[i][1]*log(yy[i]-yy[sec[i]])/2;
        }
    }
    
    return sum;
}

void gradient(int n1, int **N, double yy[], double grad[])
{
    int i;
    
    for (i=1;i<=n1;i++)
        grad[i]=0;
    
    
    for (i=1;i<=n1;i++)
    {
        if (N[i][0]>0)
            grad[i] += N[i][0]/yy[i];
        if (N[i][2]>0)
            grad[i] -= N[i][2]/(1-yy[i]);
        if (N[i][1]>0)
        {
            if (sec[i]>i)
                grad[i] -= N[i][1]/(yy[sec[i]]-yy[i]);
            else
                grad[i] += N[i][1]/(yy[i]-yy[sec[i]]);
        }
    }
}


void weights(int n1, int **N, double yy[], double w[])
{
    int i,j;
    
    for (j=1;j<=n1;j++)
        w[j]=0;
    
    for (i=1;i<=n1;i++)
    {
        if (N[i][0]>0)
            w[i] += N[i][0]/SQR(yy[i]);
        if (N[i][2]>0)
            w[i] += N[i][2]/SQR((1-yy[i]));
        if (N[i][1]>0)
            w[i] += N[i][1]/SQR(yy[sec[i]]-yy[i]);
    }
}


void cumsum(int n1, double yy[], double cumw[], double grad[], double w[])
{
    int    j;
    
    cumw[0]=0;
    cs[0]=0;
    
    for (j=1;j<=n1;j++)
    {
        cumw[j] = cumw[j-1]+w[j];
        cs[j]   = cs[j-1]+yy[j]*w[j]+grad[j];
    }
    
}

void convexminorant(int n1, double cumw[], double cs[], double yy[])
{
    int    i,j,m;
    
    yy[1] = cs[1]/cumw[1];
    for (i=2;i<=n1;i++)
    {
        yy[i] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
        if (yy[i-1]>yy[i])
        {
            j = i;
            while ((yy[j-1] > yy[i]) && (j>1))
            {
                j--;
                yy[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
                for (m=j;m<i;m++)
                    yy[m] = yy[i];
            }
        }
    }
    
    for (i=1;i<=n1;i++)
    {
        if (yy[i]<=0)
            yy[i]=0;
        if (yy[i]>=1)
            yy[i]=1;
    }
}

void isoreg(int n1, int **N, double yy[], double F[], double cumw[], double grad[], double w[])
{
    int i,iter;
    double tol=1.0e-10;
    
    gradient(n1,N,F,grad);
    
    iter=0;
    while (fenchelviol(n1,F,grad,tol,&inprod,&partsum) && iter<=1000)
    {
        iter++;
        transfer(0,n1,F,yy);
        gradient(n1,N,yy,grad);
        weights(n1,N,yy,w);
        cumsum(n1,yy,cumw,grad,w);
        convexminorant(n1,cumw,cs,yy);

        alpha=golden(n1,N,f_alpha);
        
        for (i=1;i<=n1+1;i++)
            F[i] = alpha*yy[i]+(1-alpha)*F[i];
    
    }
}

void   data_exp(int n, double data1[], double data2[], double data3[], int seed)
{
    int    i;
    double x,u,v;
    
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> dis_unif(0,1);
    
    for (i=0;i<n;i++)
    {
        x = M*dis_unif(gen);
        data3[i] = x;
        
        u = dis_unif(gen);
        v = x*dis_unif(gen);
        data2[i] = v+Weibull_inv(a,b,u);
        
        if (data2[i]<=x)
            data1[i]=0;
        else
            data1[i]=data2[i]-x;
    }
}


int compute_mle(int n, double data1[], double data2[], double F[], double tt[], double pp[])
{
    int i,j,n1;
    //FILE    *outputfile = NULL;
    
    for (i=0;i<n;i++)
    {
        data0[i]=data1[i];
        data0[n+i]=data2[i];
    }
    
    sort_data0(2*n,data0,index0);
    
    sort_index(2*n,index0);
    
    // index0 maps indices of data1 and data2 to indices of data0
    // The indices of data2 are shifted to n+i
    // We have: data0[index0[i]]=data1[i] and data0[index0[n+i]]=data2[i],
    // for (i=0,..,n.
    
    min_obs2= 1.0e10;
    for (i=0;i<n;i++)
    {
        if (data2[i]<min_obs2)
            min_obs2 = data2[i];
    }
    
    max_obs1= 0;
    for (i=0;i<n;i++)
    {
        if (data1[i]>max_obs1)
            max_obs1 = data1[i];
    }
    
    //printf("max_obs1 = %15.10f\n\n",max_obs1);
    
    //printf("min_obs2 = %15.10f\n\n",min_obs2);
    
    for (i=0;i<2*n;i++)
    {
        if (data0[i]<min_obs2)
            ind[i]=0;
    }
    
    tt[0]=0;
    ind[0]=0;
    
    j=0;
    
    for (i=0;i<2*n;i++)
    {
        if (data0[i]<min_obs2)
        {
            ind[i]=0;
        }
        else
        {
            if (data0[i]>=min_obs2 && data0[i]<=max_obs1)
            {
                if (i>0 && data0[i]>data0[i-1])
                    j++;
                ind[i]=j;
                tt[j]=data0[i];
            }
        }
    }
    
    n1=j;
    
    for (i=0;i<2*n;i++)
    {
        if (data0[i]>max_obs1)
            ind[i]=n1+1;
    }
    
    tt[n1+1]=max_obs1+1;
    
    for (i=0;i<2*n;i++)
        ind1[i]=ind[index0[i]];

    
    for (i=0;i<=2*n+1;i++)
    {
        for (j=0;j<3;j++)
            N[i][j]=0;
    }
    
    for (i=0;i<n;i++)
    {
        if (ind1[n+i]<=n1)
        {
            if (data1[i]<min_obs2)
            {
                N[ind1[n+i]][0]++;
                sec[ind1[n+i]]=0;
            }
            else
            {
                if (data1[i]>=min_obs2 && data2[i]<=max_obs1)
                {
                    N[ind1[i]][1]++;
                    N[ind1[n+i]][1]++;
                    sec[ind1[i]]=ind1[n+i];
                    sec[ind1[n+i]]=ind1[i];
                }
            }
        }
        else
        {
            if (data1[i]>=min_obs2)
            {
                N[ind1[i]][2]++;
                sec[ind1[i]]=n1+1;
            }
        }
    }
    
    F[0]=0.0;
    for (i=1;i<=n1;i++)
         F[i]=i*1.0/(n1+1);
     
    F[n1+1]=1;
     
    yy[0]=yy_new[0]=0;
    yy[n1+1]=yy_new[n1+1]=1;
     
    isoreg(n1,N,yy,F,cumw,grad,w);
    
    pp[0]=0;
    for (i=1;i<=n1+1;i++)
        pp[i] = F[i]-F[i-1];
    
    /*for (i=0;i<n1+2;i++)
        delete[] N[i];
    delete[] N;*/
    
    return n1;
}


int compare(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int compare2(const void *a, const void *b)
{
    double x = *(double*)a;
    double y = *(double*)b;
    
    if (x < y)
        return -1;
    if (x > y)
        return 1;
    return 0;
}
void sort_data0(int m, double data0[], int index0[])
{
    int i;
    data0_object *obs;
    
    obs = new data0_object[m];
    
    for (i=0;i<m;i++)
    {
        obs[i].x = data0[i];
        obs[i].index = i;
    }
    
    qsort(obs,m,sizeof(data0_object),compare2);
    
    for (i=0;i<m;i++)
    {
        data0[i]=obs[i].x;
        index0[i]=obs[i].index;
    }
    
    delete[] obs;
}

void sort_data(int m, double data1[], double data2[])
{
    int i;
    data_object *obs;
    
    obs = new data_object[m];
    
    for (i=0;i<m;i++)
    {
        obs[i].x = data1[i];
        obs[i].y = data2[i];
    }
    
    qsort(obs,m,sizeof(data_object),compare2);
    
    for (i=0;i<m;i++)
    {
        data1[i]=obs[i].x;
        data2[i]=obs[i].y;
    }
    
    delete[] obs;
}

void sort_index(int m, int index1[])
{
    int i;
    index_object *obs;
    
    obs = new index_object[m];
    
    for (i=0;i<m;i++)
    {
        obs[i].i = index1[i];
        obs[i].j = i;
    }
    
    qsort(obs,m,sizeof(index_object),compare);
    
    for (i=0;i<m;i++)
        index1[i]=obs[i].j;
    
    delete[] obs;
}

double KK(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y = (16.0 + 35*x - 35*pow(x,3) + 21*pow(x,5) - 5*pow(x,7))/32.0;
    else
    {
        if (x>1)
            y=1;
        else
            y=0;
        
    }
    
    return y;
}

double K(double x)
{
    double u,y;
    
    u=x*x;
    
    if (u<=1)
        y=(35.0/32)*pow(1-u,3);
    else
        y=0.0;
    
    return y;
}

double dens_estimate(double A, double B,  int m, double t[], double p[], double u, double h)
{
    int k;
    double      t1,sum;
    
    sum=0;
    
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        sum += K(t1)*p[k]/h;
    }
    
    return fmax(0,sum);
}


double Weibull_df(double a, double b, double x)
{
    if (x>0 && x<= M1)
        return 1-exp(-b*pow(x,a))/(1-exp(-b*pow(M1,a)));
    else
    {
        if (x>M1)
            return 1;
        else
            return 0;
    }
}

double Weibull_inv(double a, double b, double u)
{

    if (u<1 && u>0)
        return pow(-log((1 - exp(-b*pow(M1,a)))*(1/(1-exp(-b*pow(M1,a)))-u))/b,1.0/a);
    else
        return 0;

    //return pow(-log(1-u)/b,1.0/a);;
}


double Weibull_dens(double a, double b, double x)
{
    if (x>0 || x<=M1)
        return a*b*exp(-b*pow(x,a))*pow(x,a-1)/(1-exp(-b*pow(M1,a)));
    else
        return 0;
}
