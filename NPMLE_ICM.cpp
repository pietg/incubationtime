//
//  The iterative convex minaorant algorithm for computing the MLE for estimating
//  the distribution of the incubation time
//
//  Created by Piet Groeneboom on 16/07/2020.
//  Copyright (c) 2020 Piet Groeneboom. All rights reserved.

//

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#define SQR(x) ((x)*(x))

int n,n2;
double alpha;
int *tt;
double  *dens,*pp,*grid,*F,*MLE,*SMLE,h;
double  *cumw,*cs,*yy,*yy_new,*grad,*w;
double  inprod,partsum;

double  f_alpha(double alpha);
double  f_alpha_prime(double alpha);
int     fenchelviol(double b[], double nabla[], double tol, double *inprod, double *partsum);

double  criterion(double yy[]);
void    convexminorant(double cumw[], double cs[], double yy[]);
void    isoreg(double yy[], double F[], double cumw[], double grad[], double w[]);
void    cumsum(double yy[], double cumw[], double grad[], double w[]);
void    weights(double yy[], double w[]);
void    gradient(double yy[], double grad[]);
void    transfer(int first, int last, double a[], double b[]);
double  golden(double (*f)(double));

double bdf(double A, double B, int m, int t[], double p[], double u, double h);
double dens_estimate(double A, double B,  int m, int t[], double p[], double u, double h);
double KK(double x);
double K(double x);


// [[Rcpp::export]]

List NPMLE()
{
    int     i,ngrid;
    double  a,sum,max,locmax,endpoint=14;

    // number of parameters to be estimated
    
    n=6;
    
    // number of points needed for MLE, extended to the left (by two points)
    // and to the right (by one point, where it becomes equal to 1)
    
    n2=9;
    
     F =  new double[n+1];
     MLE = new double[n2+1];
     pp = new double[n+2];
     tt = new int[n+2];
     
     ngrid=140;

     grid = new double[ngrid+1];
     SMLE = new double[ngrid+1];
     dens = new double[ngrid+1];
     
     for (i=0;i<=ngrid;i++)
         grid[i] = i*0.1;
     
     cumw = new double[n+1];
     cs = new double[n+1];
     yy = new double[n+1];
     yy_new = new double[n+1];
     grad = new double[n+1];
     w = new double[n+1];
     
     for (i=1;i<=n;i++)
         tt[i]=i+2;
     
     for (i=0;i<=n;i++)
         F[i]=i*1.0/(n+1);
     
     isoreg(yy,F,cumw,grad,w);
     
     MLE[0]=MLE[1]=MLE[2]=0;
    
     for (i=1;i<=n;i++)
         MLE[i+2] = F[i];
    
    MLE[n2]=1;
    
    NumericMatrix out1 = NumericMatrix(n2+2,2);
    
    for (i=0;i<=n2;i++)
    {
        out1(i,0)=i;
        out1(i,1)=MLE[i];
    }
    
    out1(n2+1,0)=endpoint;
    out1(n2+1,1)=1;
    
    pp[0]=0;
    for (i=1;i<=n;i++)
        pp[i] = F[i]-F[i-1];
    for (i=1;i<=n;i++)
        tt[i] = i+2;
    
    pp[n+1]=1-F[n];
    tt[n+1]=n+3;
    
    sum=0;
    for (i=1;i<=n+1;i++)
        sum += tt[i]*pp[i];
    
    double out2 = sum;
    
    max=locmax=0;
    
    for (i=0;i<=ngrid;i++)
    {
        SMLE[i] = bdf(0.0,endpoint,n+1,tt,pp,grid[i],3.0);
        dens[i]=dens_estimate(0.0,14.0,n+1,tt,pp,grid[i],4.0);

        a = dens[i];
        if (max<a)
        {
            max=a;
            locmax=grid[i];
        }
    }

    double out3 = locmax;
    
    NumericMatrix out4 = NumericMatrix(ngrid+1,2);
    
    for (i=0;i<=ngrid;i++)
    {
        out4(i,0)=grid[i];
        out4(i,1)=SMLE[i];
    }
    
    NumericMatrix out5 = NumericMatrix(ngrid+1,2);
     
     for (i=0;i<=ngrid;i++)
     {
         out5(i,0)=grid[i];
         out5(i,1)=dens[i];
     }

    
    // make the list for the output
    
    List out = List::create(Named("MLE")=out1,Named("mean")=out2,Named("locmax")=out3,Named("SMLE")=out4,Named("dens")=out5);

    // free memory
    
    delete[] dens; delete[] pp; delete[] grid; delete[] F; delete[] MLE; delete[] SMLE;
    delete[] cumw; delete[] cs; delete[] yy; delete[] grad; delete[] w, delete[] tt;
    
    return out;
}

void isoreg(double yy[], double F[], double cumw[], double grad[], double w[])
{
    int i;
    double tol=1.0e-10;
    
    gradient(F,grad);
    
    while (fenchelviol(F,grad,tol,&inprod,&partsum))
    {
        transfer(1,n,F,yy);
        gradient(yy,grad);
        weights(yy,w);
        cumsum(yy,cumw,grad,w);
        convexminorant(cumw,cs,yy);
        
        if (f_alpha_prime(alpha)<=0)
            alpha=1;
        else
            alpha=golden(f_alpha);
        
        for (i=1;i<=n;i++)
            F[i] = alpha*yy[i]+(1-alpha)*F[i];
    }
}

double f_alpha(double alpha)
{
    int i;
    
    for (i=1;i<=n;i++)
        yy_new[i]=(1-alpha)*F[i]+alpha*yy[i];

    return criterion(yy_new);
}

double f_alpha_prime(double alpha)
{
    int        i;
    double    sum;
    
    sum=0;
            
    for (i=1;i<=n;i++)
            sum -= (yy[i]-F[i])/((1-alpha)*F[i]+alpha*yy[i]);
    
    return sum;
}


int fenchelviol(double yy[], double grad[], double tol, double *inprod, double *partsum)
{
    double    sum,sum2;
    int    i;
    int    fenchelvioltemp;
    
    fenchelvioltemp = 0;
    
    sum=0;
    sum2=0;
    
    for (i=1;i<=n;i++)
    {
        sum += grad[i];
        if (sum < sum2)
            sum2 = sum;
    }
    
    sum=0;
    for (i=1;i<=n;i++)
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

double criterion(double yy[])
{
    double sum=0;
    
    sum -= log(yy[1])+3*log(yy[2])+4*log(yy[3])+2*log(yy[6])+2*log(yy[2]-yy[1])
            +log(yy[3]-yy[1])+log(yy[4]-yy[3])+log(yy[4]-yy[2])+log(yy[5]-yy[4])+
            +log(yy[5]-yy[2])+2*log(yy[6]-yy[3])+log(yy[6]-yy[5])+9*log(1-yy[1])
            +4*log(1-yy[2])+3*log(1-yy[3])+6*log(1-yy[4])+3*log(1-yy[5])+3*log(1-yy[6]);
    
    return sum;
}

void gradient(double yy[], double grad[])
{
    grad[1] = 1/yy[1]-9/(1-yy[1])-2/(yy[2]-yy[1])-1/(yy[3]-yy[1]);
    grad[2] = 3/yy[2]-4/(1-yy[2])+2/(yy[2]-yy[1])-1/(yy[4]-yy[2])-1/(yy[5]-yy[2]);
    grad[3] = 4/yy[3]-3/(1-yy[3])+1/(yy[3]-yy[1])-1/(yy[4]-yy[3])-1/(yy[6]-yy[3]);
    grad[4] = 1/(yy[4]-yy[2])-6/(1-yy[4])+1/(yy[4]-yy[3])-1/(yy[5]-yy[4]);
    grad[5] = 1/(yy[5]-yy[2])-3/(1-yy[5])+1/(yy[5]-yy[4])-1/(yy[6]-yy[5]);
    grad[6] = 2/(yy[5]-yy[2])-3/(1-yy[6])+2/(yy[6]-yy[3])+1/(yy[6]-yy[5]);
}


void weights(double yy[], double w[])
{
   w[1] = 1/SQR(yy[1])+9/SQR(1-yy[1])+2/SQR(yy[2]-yy[1])+1/SQR(yy[3]-yy[1]);
   w[2] = 3/SQR(yy[2])+4/SQR(1-yy[2])+2/SQR(yy[2]-yy[1])+1/SQR(yy[4]-yy[2])+1/SQR(yy[5]-yy[2]);
   w[3] = 4/SQR(yy[3])+3/SQR(1-yy[3])+1/SQR(yy[3]-yy[1])+1/SQR(yy[4]-yy[3])+1/SQR(yy[6]-yy[3]);
   w[4] = 1/SQR(yy[4]-yy[2])+6/SQR(1-yy[4])+1/SQR(yy[4]-yy[3])+1/SQR(yy[5]-yy[4]);
   w[5] = 1/SQR(yy[5]-yy[2])+3/SQR(1-yy[5])+1/SQR(yy[5]-yy[4])+1/SQR(yy[6]-yy[5]);
   w[6] = 2/SQR(yy[5]-yy[2])+3/SQR(1-yy[6])+2/SQR(yy[6]-yy[3])+1/SQR(yy[6]-yy[5]);
}


void cumsum(double yy[], double cumw[], double grad[], double w[])
{
    int    j;
    
    cumw[0]=0;
    cs[0]=0;
    
    for (j=1;j<=n;j++)
    {
        cumw[j] = cumw[j-1]+w[j];
        cs[j]   = cs[j-1]+yy[j]*w[j]+grad[j];
    }
    
}

void convexminorant(double cumw[], double cs[], double yy[])
{
    int    i,j,m;
    
    yy[1] = cs[1]/cumw[1];
    for (i=1+1;i<=n;i++)
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
}

double golden(double (*f)(double))
{
    double a,b,eps=1.0e-10;
    
    a=0;
    b=1;
    
    double k = (sqrt(5.0) - 1.0) / 2;
    double xL = b - k*(b - a);
    double xR = a + k*(b - a);
    
    while (b-a>eps)
    {
        if ((*f)(xL)<(*f)(xR))
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


double bdf(double A, double B, int m, int t[], double p[], double u, double h)
{
    int           k;
    double        t1,t2,t3,sum;
    
    
    sum=0;
    for (k=1;k<=m;k++)
    {
        t1=(u-t[k])/h;
        t2=(u+t[k]-2*A)/h;
        t3=(2*B-u-t[k])/h;
        sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
    }
    return sum;
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
    }
    
    return fmax(0,sum);
}


