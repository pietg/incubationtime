//
//  The Weibull maximum likelihood method for estimating the distribution of the incubation time
//
//  Created by Piet Groeneboom on 16/07/2020.
//  Copyright (c) 2020 Piet Groeneboom. All rights reserved.
//

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

double  K1(double t, double w, double h);
double  Weibull_df(double x, double a, double b, double M1);
int     LUPDecompose(double **A, int N, double Tol, int *P);
void    LUPSolve(double **A, int *P, double *b, int N, double *x);


// [[Rcpp::export]]

List inteq(double point)
{
    int     i,j,ngrid,*P;
    double  *phi,*psi,M,M1;
    double  alpha,beta,**A;
    double  *grid,tol=1.0e-9;
    double  step,t,e,w,bandwidth,variance;
    FILE         *outputfile = NULL;
    char         filename[128];
    
    // point of evaluation
    
    t = (double)point;
    
    if (t<2 || t>11)
        Rcpp::stop("Please take a value between 2 and 11!\n");
    
    alpha=3.035140901;
    beta=0.002619475;
    M=30.0;
    M1=20;
    
    ngrid=1999;
    step=0.01;
    
    grid = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        grid[i]=i*step;
    
    bandwidth=3.4;
    
    A = new double *[ngrid+1];
    for (i=0;i<ngrid+1;i++)
        A[i] = new double[ngrid+1];
    
    P = new int[ngrid+2];
    
    phi = new double[ngrid+1];
    psi = new double[ngrid+1];
    
    for (i=0;i<=ngrid;i++)
        phi[i]=psi[i]=0;
    
    for (i=0;i<=ngrid;i++)
        psi[i] = K1(t,grid[i],bandwidth);
    
    for (i=0;i<=ngrid;i++)
    {
        for (j=0;j<=ngrid;j++)
            A[i][j]=0;
    }
    
    A[0][0] = -log(M/grid[1])/(M*Weibull_df(grid[1],alpha,beta,M1));
    
    for (i=1;i<=ngrid;i++)
        A[i][i] = -log(M/grid[i])/(M*Weibull_df(grid[i],alpha,beta,M1));
    
    for (i=1;i<=ngrid;i++)
    {
        w=grid[i];
        for (j=1;j<=i;j++)
        {
            e=grid[j];
            if (i+j<=ngrid)
                A[i][i+j] += step/(M*e*(Weibull_df(w+e,alpha,beta,M1)-Weibull_df(w,alpha,beta,M1)));
            A[i][i] -= step/(M*e*(Weibull_df(w+e,alpha,beta,M1)-Weibull_df(w,alpha,beta,M1)));
        }
        for (j=1;j<=i;j++)
        {
            e=grid[j];
            A[i][i-j]   += step/(M*e*(Weibull_df(w,alpha,beta,M1)-Weibull_df(w-e,alpha,beta,M1)));
            A[i][i]   -= step/(M*e*(Weibull_df(w,alpha,beta,M1)-Weibull_df(w-e,alpha,beta,M1)));
        }
        
        for (j=i+1;j<=ngrid;j++)
        {
            e=grid[j];
            if (i+j<=ngrid)
                A[i][i+j]   += step/(M*e*(Weibull_df(w+e,alpha,beta,M1)-Weibull_df(w,alpha,beta,M1)));
            A[i][i]   -= step/(M*e*(Weibull_df(w+e,alpha,beta,M1)-Weibull_df(w,alpha,beta,M1)));
        }
    }
    
        
    LUPDecompose(A,ngrid+1,tol,P);
    LUPSolve(A,P,psi,ngrid+1,phi);
    
    sprintf(filename, "phi.txt");
    outputfile = fopen(filename, "w");
    rewind(outputfile);
    
    for(i=0;i<=ngrid;i++)
        fprintf(outputfile,"%12.8f   %15.10f     %15.10f\n", grid[i],phi[i],K1(t,grid[i],bandwidth));
    
    fclose(outputfile);
    
    variance=0;
    
    for (i=1;i<=ngrid;i++)
        variance += SQR(phi[i])*log(M/grid[i])*step/(M*Weibull_df(grid[i],alpha,beta,M1));
    
    for (i=1;i<=ngrid;i++)
    {
        for (j=0;j<i;j++)
            variance += SQR(phi[i]-phi[j])*SQR(step)/((grid[i]-grid[j])*M*(Weibull_df(grid[i],alpha,beta,M1)-Weibull_df(grid[j],alpha,beta,M1)));
    }
    
    
    NumericMatrix out1 = NumericMatrix(ngrid+1,2);
    
    for (i=0;i<=ngrid;i++)
    {
        out1(i,0)=grid[i];
        out1(i,1)=phi[i];
    }
    
    double out2 = pow(1000,-3.0/7)*variance;
    
    // make the list for the output, containing the two estimates and the log likelihood
    
    List out = List::create(Rcpp::Named("phi")=out1,Rcpp::Named("variance")=out2);

    // free memory
    
    for (i=0;i<ngrid+1;i++)
        delete[] A[i];
    delete[] A;
    
    delete[] phi; delete[] psi; delete[] P; delete[] grid;
    return out;
}

double Weibull_df(double x, double a, double b, double M1)
{
    if (x>0 && x<=M1)
        return (1-exp(-b*pow(x,a)))/(1-exp(-b*pow(M1,a)));
    else
    {
        if (x>M1)
            return 1;
        else return 0;
    }
}

double K1(double t, double w, double h)
{
    if (fabs((t-w)/h)<=1)
        return 105*SQR(SQR(h)-SQR(t-w))*(w-t)/(16*pow(h,7));
    else
        return 0;
}

// matrix routines from https://en.wikipedia.org/wiki/LU_decomposition

int LUPDecompose(double **A, int N, double Tol, int *P)
{
    int i, j, k, imax;
    double maxA, *ptr, absA;
    
    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    
    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;
        
        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }
        
        if (maxA < Tol) return 0; //failure, matrix is degenerate
        
        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            
            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;
            
            //counting pivots starting from N (for determinant)
            P[N]++;
        }
        
        for (j = i + 1; j < N; j++)
        {
            A[j][i] /= A[i][i];
            
            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }
    
    return 1;  //decomposition done
}

// INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
// OUTPUT: x - solution vector of A*x=b

void LUPSolve(double **A, int *P, double *b, int N, double *x)
{
    for (int i = 0; i < N; i++)
    {
        x[i] = b[P[i]];
        
        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }
    
    for (int i = N - 1; i >= 0; i--)
    {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];
        
        x[i] = x[i] / A[i][i];
    }
}



