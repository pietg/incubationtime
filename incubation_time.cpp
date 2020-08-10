//
//  Estimation of incubation time as inspired by the 2018 paper of Britton and Scalia Tomba
//
//
//  Created by Piet Groeneboom on July-20-2020.
//  Copyright (c) 2020 Piet. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <random>
#include <fstream>
#include <string.h>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

int n;
double *alpha_init,*alpha;
unsigned long seed;

int *time_index1,*time_index2,NumIt;
double  **times2;
double  kappa,lambda,crate,p;
s
double  criterion(int m, double alpha[]);
int     hooke(int m, double startpt[], double endpt[], double rho, double eps,
                int itermax, double f(int m, double alpha[]));
double best_nearby (int m, double delta[], double point[], double prevbest,
                    double f(int m, double alpha[]), int *funevals);

// [[Rcpp::export]]

List incubation(int N, int numit)
{
    int m,i,j,iter,nIterations;
    double rho,eps;
    
    n = (int)N;
    NumIt=numit;
    seed=100;
    
    m=3;
    
    crate=0.0725;
    p=0.5;
    kappa=1.9808;
    lambda=0.1738;
    
    rho=0.5;
    eps=1.0e-10;
    nIterations=100;
     
    iter=0;
    
    alpha       = new double[m];
    alpha_init  = new double[m];
    
    alpha_init[0]=0.75;
    alpha_init[1]=2.0;
    alpha_init[2]=0.2;
    
    time_index1 = new int[n];
    time_index2 = new int[n];
    
    times2 = new double *[n];
    
    NumericMatrix out1 = NumericMatrix(m,NumIt);
        
    for (j=0;j<NumIt;j++)
    {
        seed++;
        data(n);
        
        iter=hooke(m,alpha_init,alpha,rho,eps,nIterations,criterion);
        
        for (i=0;i<m;i++)
            out1(i,j)=alpha[i];
        
        Rcout <<  "\t" <<  j+1;
        
        for (i=0;i<m;i++)
            Rcout << setprecision(10) << setw(15) << alpha[i];
        
        Rcout << endl;
        
        if (j<NumIt-1)
        {
            for (i=0;i<n;i++)
                delete[] times2[i];
        }
    }
    
    List out2(n);
    
    for (i=0;i<n;i++)
    {
        NumericVector time(time_index2[i]+1);
        
        time[0]=0;
        for (j=1;j<time_index2[i]+1;j++)
            time[j]=times2[i][j];
        out2[i] = time;
    }
    
    IntegerVector out3 = IntegerVector(n);
    for (i=0;i<n;i++)
        out3[i]=time_index2[i];
      
    List out = List::create(Rcpp::Named("parameters")=out1,Rcpp::Named("sample")=out2,Rcpp::Named("number_of_times")=out3);
    
    delete[] times2;
       
    delete[] time_index1; delete[] time_index2;
    delete[] alpha_init; delete[] alpha;
    
    return out;
}

void data(int n)
{
    int    i,j,k,f,m;
    double  s,t;

    vector<vector<double>> times(n);
    
    mt19937_64 gen(seed);
    geometric_distribution<int> geodis(p);
    gamma_distribution<double> gammadis(kappa,1.0/lambda);
    exponential_distribution<> expdis(crate);
    
    for (i=0;i<n;i++)
    {
        k=m=geodis(gen)+1;
        times[i] = vector<double>(k);
        times[i][0]=0;
        if (k>1)
        {
            for (j=1;j<k;j++)
                times[i][j]=times[i][j-1]+expdis(gen);
        }
            
        s=times[i][k-1]+gammadis(gen);
        t=times[i][k-1]+expdis(gen);
        
        while (t<s)
        {
            times[i].push_back(t);
            k++;
            t=times[i][k-1]+expdis(gen);
        }
        
        f=k;
        
        time_index1[i]=m;
        time_index2[i]=k;
        
        times2[i]= new double[k+1];
        for (j=0;j<k;j++)
            times2[i][j]=times[i][j];
        
        times[i].push_back(s);
        times2[i][k]=times[i][k];
    }
}

double criterion(int m, double alpha[])
{
    int i,j,k;
    double sum1,sum=0.0;
    
    for (i=0;i<n;i++)
    {
        sum1=0;
        k=time_index2[i];
        for (j=0;j<k;j++)
            sum1 += alpha[0]*pow(1-alpha[0],j)*gammadens(times2[i][k]-times2[i][j],alpha[1],1.0/alpha[2]);
        sum -= log(sum1);
    }
    return sum;
}

double   gammadens(double x, double gamma1, double gamma2)
{
    if (x>0)
        return pow(x,gamma1-1)*exp(-x/gamma2)/(tgamma(gamma1)*pow(gamma2,gamma1));
    else
        return 0;
}

int hooke(int m, double startpt[], double endpt[], double rho, double eps,
          int itermax, double f(int m, double alpha[]))
{
    double *delta,fbefore;
    int i,iters,keep,funevals,count;
    double newf,*newx,steplength,tmp;
    bool verbose = false;
    double *xbefore;
    
    delta = new double[m];
    newx = new double[m];
    xbefore = new double[m];
    
    for ( i = 0; i < m; i++ )
        xbefore[i] = newx[i] = startpt[i];
    
    for ( i = 0; i < m; i++ )
    {
        if ( startpt[i] == 0.0 )
            delta[i] = rho;
        else
            delta[i] = rho*fabs(startpt[i]);
    }
    
    funevals = 0;
    steplength = rho;
    iters = 0;
    
    
    fbefore = f(m,newx);
    funevals = funevals + 1;
    newf = fbefore;
    
    while ( iters < itermax && eps < steplength )
    {
        iters = iters + 1;
        
        if (verbose)
        {
            cout << "\n";
            cout << "  FUNEVALS, = " << funevals
            << "  F(X) = " << fbefore << "\n";
            
            for ( i = 0; i < m; i++ )
            {
                cout << "  " << i + 1
                << "  " << xbefore[i] << "\n";
            }
        }
        //
        //  Find best new alpha, one coordinate at a time.
        //
        for ( i = 0; i < m; i++ )
            newx[i] = xbefore[i];
        
        
        
        newf = best_nearby(m,delta,newx,fbefore,f,&funevals);
        //
        //  If we made some improvements, pursue that direction.
        //
        keep = 1;
        count=0;
        
        while (newf<fbefore && keep == 1 && count<=100)
        {
            count++;
            for ( i = 0; i < m; i++ )
            {
                //
                //  Arrange the sign of DELTA.
                //
                if ( newx[i] <= xbefore[i] )
                    delta[i] = - fabs(delta[i]);
                else
                    delta[i] = fabs(delta[i]);
                //
                //  Now, move further in this direction.
                //
                tmp = xbefore[i];
                xbefore[i] = newx[i];
                newx[i] = newx[i] + newx[i] - tmp;
            }
            
            fbefore = newf;
            
            newf = best_nearby(m,delta,newx,fbefore,f,&funevals);
            //
            //  If the further (optimistic) move was bad...
            //
            if (fbefore <= newf)
                break;
            //
            //  Make sure that the differences between the new and the old points
            //  are due to actual displacements; beware of roundoff errors that
            //  might cause NEWF < FBEFORE.
            //
            keep = 0;
            
            for ( i = 0; i < m; i++ )
            {
                if ( 0.5 * fabs(delta[i]) < fabs(newx[i]-xbefore[i]))
                {
                    keep = 1;
                    break;
                }
            }
        }
        
        if (eps <= steplength && fbefore <= newf)
        {
            steplength = steplength * rho;
            for ( i = 0; i < m; i++ )
                delta[i] = delta[i] * rho;
        }
        
    }
    
    for ( i = 0; i < m; i++ )
        endpt[i] = xbefore[i];
    
    delete [] delta;
    delete [] newx;
    delete [] xbefore;
    
    return iters;
}

double best_nearby(int m, double delta[], double point[], double prevbest,
                    double f(int m, double alpha[]), int *funevals)
{
    double ftmp,minf,*z;
    int i;
    
    z = new double[m];
    
    minf = prevbest;
    
    for ( i = 0; i < m; i++ )
        z[i] = point[i];
    
    for ( i = 0; i < m; i++ )
    {
        z[i] = point[i] + delta[i];
        
        ftmp = f(m,z);
        *funevals = *funevals + 1;
        
        if ( ftmp < minf )
            minf = ftmp;
        else
        {
            delta[i] = - delta[i];
            z[i] = point[i] + delta[i];
            ftmp = f(m,z);
            *funevals = *funevals + 1;
            
            if ( ftmp < minf )
                minf = ftmp;
            else
                z[i] = point[i];
        }
    }
    
    for ( i = 0; i < m; i++ )
        point[i] = z[i];
    
    delete [] z;
    
    return minf;
}

