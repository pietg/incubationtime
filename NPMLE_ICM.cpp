//
//  The iterative convex minaorant algorithm for computing the MLE for estimating
//  the distribution of the incubation time
//
//  Created by Piet Groeneboom on 16/07/2020.
//  Copyright (c) 2020 Piet Groeneboom. All rights reserved.


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

int     n,n1,n2,ngrid,nloc,*tt,**N;
double  alpha,*F;
int     *data1,*data2;
double  *pp,*grid,*MLE,*SMLE,*dens;
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

List NPMLE(DataFrame input)
{
    int     i,j;
    double  a,sum,max,locmax,endpoint=14;

    // number of parameters to be estimated
    
    n1=6;
    
    // number of points needed for MLE, extended to the left (by two points)
    // and to the right (by one point, where it becomes equal to 1)
    
    n2=9;
    
    DataFrame DF = Rcpp::DataFrame(input);
    IntegerVector data01 = DF["V1"];
    IntegerVector data02 = DF["V2"];
    
    // determine the sample size
    
    n = (int)data01.size();
    
    data1 = new int[n];
    data2 = new int[n];
    
    for (i=0;i<n;i++)
    {
        data1[i]=data01[i];
        data2[i]=data02[i];
    }
    
    F =  new double[n1+2];
    MLE = new double[n1+4];
    pp = new double[n1+2];
    tt = new int[n1+2];
      
    cumw = new double[n1+1];
    cs = new double[n1+1];
    yy = new double[n1+2];
    yy_new = new double[n1+2];
    grad = new double[n1+1];
    w = new double[n1+1];
    
    ngrid=140;

    grid = new double[ngrid+1];
    SMLE = new double[ngrid+1];
    dens = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i] = i*0.1;
     
    N = new int *[n1+2];
    for (i=0;i<n1+2;i++)
        N[i] = new int[n1+2];
    
    for (i=0;i<n1+2;i++)
    {
        for (j=0;j<n1+2;j++)
            N[i][j]=0;
    }
    
    for (i=0;i<n;i++)
    {
        if (data2[i]<=8)
        {
            if (data1[i]<=2)
                N[0][data2[i]-2]++;
            else
            {
                if (data2[i]>2 && data2[i]<=8)
                    N[data1[i]-2][data2[i]-2]++;
            }
        }
    }
    
    for (i=0;i<n;i++)
    {
        if (data2[i]>8 && data1[i]>2)
                N[data1[i]-2][7]++;
    }
     
    for (i=1;i<=n1;i++)
         tt[i]=i+2;
     
    for (i=0;i<=n1;i++)
         F[i]=i*1.0/(n1+1);
    
    yy[0]=yy_new[0]=0;
    yy[n1+1]=yy_new[n1+1]=1;
     
    isoreg(yy,F,cumw,grad,w);
     
     MLE[0]=MLE[1]=MLE[2]=0;
    
     for (i=1;i<=n1;i++)
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
    for (i=1;i<=n1;i++)
        pp[i] = F[i]-F[i-1];
    for (i=1;i<=n1;i++)
        tt[i] = i+2;
    
    pp[n1+1]=1-F[n1];
    tt[n1+1]=n1+3;
    
    sum=0;
    for (i=1;i<=n1+1;i++)
        sum += tt[i]*pp[i];
    
    double out2 = sum;
    
    max=locmax=0;
    
    for (i=0;i<=ngrid;i++)
    {
        SMLE[i] = bdf(0.0,endpoint,n1+1,tt,pp,grid[i],3.0);
        dens[i]=dens_estimate(0.0,14.0,n1+1,tt,pp,grid[i],4.0);

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
    
    for (i=0;i<n1+2;i++)
        delete[] N[i];
    delete[] N;
        
    delete[] dens; delete[] pp; delete[] grid; delete[] F; delete[] MLE; delete[] SMLE;
    delete[] cumw; delete[] cs; delete[] yy; delete[] yy_new; delete[] grad; delete[] w, delete[] tt;
    delete[] data1; delete[] data2;
    
    return out;
}

void isoreg(double yy[], double F[], double cumw[], double grad[], double w[])
{
    int i,iter;
    double tol=1.0e-10;
    
    gradient(F,grad);
    
    iter=0;
    while (fenchelviol(F,grad,tol,&inprod,&partsum) && iter<=1000)
    {
        iter++;
        transfer(1,n1,F,yy);
        gradient(yy,grad);
        weights(yy,w);
        cumsum(yy,cumw,grad,w);
        convexminorant(cumw,cs,yy);
        
        if (f_alpha_prime(alpha)<=0)
            alpha=1;
        else
            alpha=golden(f_alpha);
        
        for (i=1;i<=n1;i++)
            F[i] = alpha*yy[i]+(1-alpha)*F[i];
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

double f_alpha(double alpha)
{
    int i;
    
    for (i=1;i<=n1;i++)
        yy_new[i]=(1-alpha)*F[i]+alpha*yy[i];

    return criterion(yy_new);
}

double f_alpha_prime(double alpha)
{
    int        i;
    double    sum;
    
    sum=0;
            
    for (i=1;i<=n1;i++)
            sum -= (yy[i]-F[i])*grad[i]/((1-alpha)*F[i]+alpha*yy[i]);
    
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

double criterion(double yy[])
{
    int i,j;
    double sum=0;
    
    for (i=0;i<=6;i++)
    {
        for (j=i+1;j<=7;j++)
        {
            if (N[i][j]>0)
                sum -= N[i][j]*log(yy[j]-yy[i]);
        }
    }
    
    return sum;
}

void gradient(double yy[], double grad[])
{
    int i,j;
    
    for (i=1;i<=6;i++)
        grad[i]=0;
    
    
    for (i=1;i<=6;i++)
    {
        for (j=i+1;j<=7;j++)
        {
            if (N[i][j]>0)
                grad[i] -= N[i][j]/(yy[j]-yy[i]);
        }
    }
        
    for (i=1;i<=6;i++)
    {
        for (j=0;j<i;j++)
        {
            if (N[j][i]>0)
                grad[i] += N[j][i]/(yy[i]-yy[j]);
        }
    }
}


void weights(double yy[], double w[])
{
    int i,j;
    
    for (j=1;j<=6;j++)
        w[j]=0;
    
    
    for (i=1;i<=6;i++)
    {
        for (j=i+1;j<=7;j++)
        {
            if (N[i][j]>0)
                w[i] += N[i][j]/SQR(yy[j]-yy[i]);
        }
    }
        
    for (i=1;i<=6;i++)
    {
        for (j=0;j<i;j++)
        {
            if (N[j][i]>0)
                w[i] += N[j][i]/SQR(yy[i]-yy[j]);
        }
    }
}


void cumsum(double yy[], double cumw[], double grad[], double w[])
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

void convexminorant(double cumw[], double cs[], double yy[])
{
    int    i,j,m;
    
    yy[1] = cs[1]/cumw[1];
    for (i=1+1;i<=n1;i++)
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


double bdf(double A, double B, int m, int t[], double p[], double u, double h)
{
    int            k;
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


