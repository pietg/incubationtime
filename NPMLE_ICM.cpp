//
//  The iterative convex minaorant algorithm for computing the MLE for estimating
//  the distribution of the incubation time
//
//  Created by Piet Groeneboom on 23/09/2020.
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
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#define SQR(x) ((x)*(x))

int     n,n1,ngrid,nloc,*loc,*tt,n1_bootstrap,*tt_bootstrap,max_obs1,min_obs2;
double  alpha,*F,*F_bootstrap,*SMLE_bootstrap,*dens_bootstrap;
int     *index1,*data0,*data1,*data2,*data3,*data1_bootstrap,*data2_bootstrap;
double  a,b,*pp,*pp_bootstrap,*grid,*SMLE,*dens,*qq;
double  *cumw,*cs,*yy,*yy_new,*grad,*w,*data_infection;
double  inprod,partsum;

int     compare (const void * a, const void * b);
double  KK(double x);
double  K(double x);
double  dens_estimate(double A, double B,  int njumps, int jumploc[], double p[], double h, double u);
double  bdf(double A, double B,  int m, int tt[], double p[], double h, double u);
int     compare(const void *a, const void *b);
void    sort_data(int n, int data0[]);

double  criterion(int n1, int **M, double yy[]);
double  f_alpha(int n1, int **N, double alpha);
double  f_alpha_prime(int n1, int **N, double alpha);
int     fenchelviol(int n1, double yy[], double grad[], double tol, double *inprod, double *partsum);
void    isoreg(int n1, int **N, double yy[], double F[], double cumw[], double grad[], double w[]);
void    cumsum(int n1, double yy[], double cumw[], double grad[], double w[]);
void    weights(int n1, int **N, double yy[], double w[]);
void    gradient(int n1, int **N, double yy[], double grad[]);
void    transfer(int first, int last, double a[], double b[]);
double  golden(int n1, int **N, double (*f)(int,int**,double));

int     compute_mle(int n, int data1[], int data2[], double F[], int tt[], double pp[]);
void    data_bootstrap(int n, int data1[], int data2[], int data1_bootstrap[], int data2_bootstrap[], int seed);
double  MSE_SMLE(int n, int B, double h1, double h2);
double  MSE_dens(int n, int B, double h1, double h2);

void    data_infect(int n, double data_infection[], int seed);
double  data_smooth();


// [[Rcpp::export]]

List NPMLE(DataFrame input)
{
    int i,B,seed;
    double  a1,sum,min,max,locmax,locmin,MSE;
    
    // number of bootstrap samples
    
    B=1000;
    
    DataFrame DF = Rcpp::DataFrame(input);
    IntegerVector data01 = DF["V1"];
    IntegerVector data02 = DF["V2"];
    
    // determine the sample size
    
    n = (int)data01.size();
    
    seed = 2;
    
    data1 = new int[n];
    data2 = new int[n];
    data3 = new int[n];
    
    for (i=0;i<n;i++)
    {
        data1[i]=data01[i];
        data2[i]=data02[i];
        data3[i]=data2[i]-data1[i];
    }
    
    data1_bootstrap = new int[n];
    data2_bootstrap = new int[n];
    data_infection  = new double[n];
    data0 = new int[2*n];
    index1 = new int[n];
    
    F =  new double[n+2];
    F_bootstrap =  new double[n+2];
    
    pp = new double[n+2];
    tt = new int[n+2];
    pp_bootstrap = new double[n+2];
    tt_bootstrap = new int[n+2];
      
    cumw = new double[n+1];
    cs = new double[n+1];
    yy = new double[n+2];
    yy_new = new double[n+2];
    grad = new double[n+1];
    w = new double[n+1];
    
    ngrid=140;

    grid = new double[ngrid+1];
    SMLE = new double[ngrid+1];
    dens = new double[ngrid+1];
    SMLE_bootstrap = new double[ngrid+1];
    dens_bootstrap = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i] = i*0.1;
    
    n1 = compute_mle(n,data1,data2,F,tt,pp);
    
    sum=max=locmax=0;
    for (i=0;i<=ngrid;i++)
        dens[i]= dens_estimate(0.0,14.0,n1+1,tt,pp,grid[i],4);
    
    min=1.0e10;
    locmin=0;
    
    Rcout << endl;
    
    Rcout << "MSE of density estimates for different bandwidths:" << endl << endl;
    
    for (i=0;i<=20;i++)
    {
        MSE = MSE_dens(n,B,3.0,3+i*0.2);
        
        Rcout << setprecision(10) << setw(15) << 3+i*0.2 << setprecision(10) << setw(15) << MSE << endl;
        
        if (MSE<min)
        {
            min = MSE;
            locmin=3+i*0.2;
        }
    }
    
    Rcout << endl;
    
    for (i=0;i<=ngrid;i++)
    {
        SMLE[i] = bdf(0.0,14.0,n1+1,tt,pp,grid[i],3.6);
        dens[i]= dens_estimate(0.0,14.0,n1+1,tt,pp,grid[i],locmin);
        
        a1 = dens[i];
        if (max<a1)
        {
            max=a1;
            locmax=grid[i];
        }
    }
    
    
    NumericMatrix out1 = NumericMatrix(n1+2,2);
    
    for (i=0;i<=n1+1;i++)
    {
        out1(i,0)=tt[i];
        out1(i,1)=F[i];
    }
    
    sum=0;
    for (i=1;i<=n1+1;i++)
        sum += tt[i]*pp[i];
    
    double out2 = sum;
    
    double out3 = locmin;
    
    double out4 = locmax;
    
    NumericMatrix out5 = NumericMatrix(ngrid+1,2);
    
    for (i=0;i<=ngrid;i++)
    {
        out5(i,0)=grid[i];
        out5(i,1)=SMLE[i];
    }
    
    NumericMatrix out6 = NumericMatrix(ngrid+1,2);
     
     for (i=0;i<=ngrid;i++)
     {
         out6(i,0)=grid[i];
         out6(i,1)=dens[i];
     }

    
    // make the list for the output
    
    List out = List::create(Named("MLE")=out1,Named("mean")=out2,Named("bandwidth")=out3,Named("locmax")=out4,Named("SMLE")=out5,Named("dens")=out6);

    // free memory
    
    delete[] dens; delete[] pp; delete[] grid; delete[] F; delete[] SMLE;
    delete[] cumw; delete[] cs; delete[] yy; delete[] yy_new; delete[] grad; delete[] w,
    delete[] tt; delete[] data1; delete[] data2; delete[] data3; delete[] data_infection;
    
    return out;
}


void data_infect(int n, double data_infection[], int seed)
{
    int    i,j;
    
    std::mt19937_64 gen(seed);
    
    for (i=0;i<n;i++)
    {
        j=data3[i];
        std::uniform_real_distribution<double> dis_unif(0.0,j*1.0);
        data_infection[i] = dis_unif(gen);
    }
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
    double	sum,sum2;
    int	i;
    int	fenchelvioltemp;
    
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
    int	i;
    for (i = first; i<= last;i++)	b[i] = a[i];
}

double criterion(int n1, int **N, double yy[])
{
    int i,j;
    double sum=0;
    
    for (i=0;i<=n1;i++)
    {
        for (j=i+1;j<=n1+1;j++)
        {
            if (N[i][j]>0)
                sum -= N[i][j]*log(yy[j]-yy[i]);
        }
    }
    
    return sum;
}

void gradient(int n1, int **N, double yy[], double grad[])
{
    int i,j;
    
    for (i=1;i<=n1;i++)
        grad[i]=0;
    
    
    for (i=1;i<=n1;i++)
    {
        for (j=i+1;j<=n1+1;j++)
        {
            if (N[i][j]>0)
                grad[i] -= N[i][j]/(yy[j]-yy[i]);
        }
    }
        
    for (i=1;i<=n1;i++)
    {
        for (j=0;j<i;j++)
        {
            if (N[j][i]>0)
                grad[i] += N[j][i]/(yy[i]-yy[j]);
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
        for (j=i+1;j<=n1+1;j++)
        {
            if (N[i][j]>0)
                w[i] += N[i][j]/SQR(yy[j]-yy[i]);
        }
    }
        
    for (i=1;i<=n1;i++)
    {
        for (j=0;j<i;j++)
        {
            if (N[j][i]>0)
                w[i] += N[j][i]/SQR(yy[i]-yy[j]);
        }
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

int compute_mle(int n, int data1[], int data2[], double F[], int tt[], double pp[])
{
    int i,j,n1,m,**N;
    
    for (i=0;i<n;i++)
    {
        data0[i]=data1[i];
        data0[i+n]=data2[i];
    }
    
    qsort(data0,2*n,sizeof(int),compare);
    
    min_obs2=(int)1.0e5;
    for (i=0;i<n;i++)
    {
        if (data2[i]<min_obs2)
            min_obs2 = data2[i];
    }
    
    min_obs2 -= 1;
    
    max_obs1=0;
    for (i=0;i<n;i++)
    {
        if (data1[i]>max_obs1)
            max_obs1 = data1[i];
    }

    
    for (i=0;i<2*n;i++)
    {
        if (data0[i]<=min_obs2)
            index1[data0[i]]=0;
    }
    
    tt[0]=0;
    
    j=0;
    for (i=1;i<2*n;i++)
    {
        if (data0[i]>data0[i-1] && data0[i]>min_obs2 && data0[i]<=max_obs1)
        {
            j++;
            index1[data0[i]]=j;
            tt[j]=data0[i];
        }
    }
    
    n1=j;
    
    for (i=0;i<2*n;i++)
    {
        if (data0[i]>max_obs1)
            index1[data0[i]]=n1+1;
    }
    
    tt[n1+1]=max_obs1+1;
    
    // n1 is number of characterizing parameters
    
    m = n1+1;
    
    N = new int *[m+1];
     
    for (i=0;i<m+1;i++)
        N[i] = new int[m+1];
    
    for (i=0;i<=m;i++)
    {
        for (j=0;j<=m;j++)
            N[i][j]=0;
    }
    
    for (i=0;i<n;i++)
    {
        if (index1[data2[i]]<=n1)
        {
            if (data1[i]<=min_obs2)
                N[0][index1[data2[i]]]++;
            else
            {
                if (data2[i]>min_obs2 && index1[data2[i]]<=m)
                    N[index1[data1[i]]][index1[data2[i]]]++;
            }
        }
    }
    
    for (i=0;i<n;i++)
    {
        if (index1[data2[i]]>n1 && data1[i]>min_obs2)
                 N[index1[data1[i]]][n1+1]++;
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
    
    for (i=0;i<m+1;i++)
        delete[] N[i];
    delete[] N;
    
    return n1;
}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

void data_bootstrap(int m, int data1[], int data2[], int data1_bootstrap[], int data2_bootstrap[], int seed)
{
    int    i,v;
        
    data_infect(n,data_infection,seed);
    
    for (i=0;i<m;i++)
    {
        v = (round)(data_infection[i]+data_smooth());
        if (v > data3[i])
        {
            data1_bootstrap[i]= v-data3[i];
            data2_bootstrap[i]=v;
        }
        else
        {
            data1_bootstrap[i]=0;
            data2_bootstrap[i]=v;
        }
    }
}

double MSE_SMLE(int n, int B, double h1, double h2)
{
    int i,iter,seed;
    double MSE;
    
    MSE=0;
    
    for (iter=1;iter<=B;iter++)
    {
        seed = rand();
        data_bootstrap(n,data1,data2,data1_bootstrap,data2_bootstrap,seed);
        n1_bootstrap = compute_mle(n,data1_bootstrap,data2_bootstrap,F_bootstrap,
                                   tt_bootstrap,pp_bootstrap);
        
        for (i=1;i<=ngrid;i++)
        {
            MSE += SQR(bdf(0.0,14.0,n1_bootstrap+1,tt_bootstrap,pp_bootstrap,grid[i],h2)
                     -bdf(0.0,14.0,n1+1,tt,pp,grid[i],h1))*0.1;
        }
    }
    
    return MSE/B;
}

double MSE_dens(int n, int B, double h1, double h2)
{
    int i,iter,seed;
    double MSE;
        
    MSE=0;
    
    for (iter=1;iter<=B;iter++)
    {
        seed = rand();
        data_bootstrap(n,data1,data2,data1_bootstrap,data2_bootstrap,seed);
        n1_bootstrap = compute_mle(n,data1_bootstrap,data2_bootstrap,F_bootstrap,
                                   tt_bootstrap,pp_bootstrap);
        
        for (i=1;i<=ngrid;i++)
            MSE += SQR(dens_estimate(0.0,14.0,n1_bootstrap+1,tt_bootstrap,pp_bootstrap,grid[i],h2)
                       -dens_estimate(0.0,14.0,n1+1,tt,pp,grid[i],h1))*0.1;
    }
    
    return MSE/B;
}

double dens_estimate(double A, double B,  int m, int t[], double p[], double u, double h)
{
    int k;
    double      t1,t2,t3,sum;
    
    sum=0;
    
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k]-2*A)/h;
        t3=(2*B-u-t[k])/h;
        sum += (K(t1)+K(t2)+K(t3))*p[k]/h;
        //sum += K(t1)*p[k]/h;
    }
    
    return fmax(0,sum);
}


double bdf(double A, double B, int m, int t[], double p[], double u, double h)
{
    int       k;
    double    t1,t2,t3,sum;
    
    sum=0;
    
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k]-2*A)/h;
        t3=(2*B-u-t[k])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
        //sum+= KK(t1)*p[k];
    }

    return fmax(sum,0);
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




double data_smooth()
{
    int i,seed;
    double v,w,c;
    
    seed = rand();
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> dis(0,14);
    std::uniform_real_distribution<double> dis_unif(0.0,1.0);
    
    i=0;
    w=0;
    
    while (i<1)
    {
        w=dis(gen);
        v=dis_unif(gen);
        c = dens_estimate(0.0,14.0,n1+1,tt,pp,w,4);
        if (v<14*c/100)
            i++;
    }
    return w;
}

