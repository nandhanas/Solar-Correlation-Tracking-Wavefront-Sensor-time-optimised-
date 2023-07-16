/*
NAME:           surfacefit
PURPOSE:        TO return an fftw_complex array of input multiplied with surfacefit constant
SAMPLE CALLING SEQUENCE:
                result =surfacefit(in,kk) 
INPUTS:         in:fftw_compex array of input values
                h :double precision array of surfacefit constant
OUTPUTS:        res:fftw_complex array of manipulated input


RETURNS:        fftw_complex* in= fftw_array of input multiplied with surfacefit constant


HISTORY:       Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/


#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include <stdlib.h>
#include<complex.h>
#define nz 1
#define degree 1

extern int nx,ny,n2,m;      //accessing global variables declared in fccorr.c
extern double *ut;          //accessing global variable declared in preprocess.h

fftw_complex* surfacefit(fftw_complex *in,double *kk)
{ 
    int i,j,k,p,q,e;

//kx = float(kk # reform( z, m, 1))
//finding kx: reform the data to m,1

   fftw_complex *reformin;
   reformin=fftw_malloc( sizeof (fftw_complex)*m*1);
   e=0;
   for(p=0;p<nx;p++)
   {
      for(q=0;q<ny;q++)
      {
          reformin[e][0]=in[p*ny+q][0];
          reformin[e][1]=in[p*ny+q][1];
          e++;
      }
   }

// # operation on kk and reform input

   double *kx;
   kx=(double*)fftw_malloc( sizeof (double)*n2*1);
   for(i=0;i<n2;i++) {
        for(j=0;j<1;j++) {
            kx[i*1+j]=0;
            for(k=0;k<m;k++) {
              kx[i*1+j]+=kk[i*m+k]*reformin[k*1+j][0];

            }
        }
    }

//here max_deg=0 so perform reform( kx,degree+1,degree+1)

   double *reformkx;
   reformkx=(double*)fftw_malloc( sizeof (double)*(degree+1)*(degree+1));
   e=0;
   for(p=0;p<(degree+1);p++)
   {
      for(q=0;q<(degree+1);q++)
      {
           reformkx[p*(degree+1)+q]=kx[e];
           e++;
      }
   }

//this gives the coefficients of the bilinear fitting

 /*printf("Coefficients\n");
     for ( i = 0; i < (degree+1); i++ )
     {
           for ( j = 0; j <(degree+1); j++ )
           {
                 printf ( "  %lf  ",reformkx[i*(degree+1)+j] );
           }
           printf("\n");
    }
 */

//reform(reform(kx,n2) # ut, nx, ny)
//changing reform kx to kx which is already there
//# operation on kx and ut

  double *kxut;
  kxut=(double*)fftw_malloc( sizeof (double)*m*1);
  for(i=0;i<1;i++) {
         for(j=0;j<m;j++) {
              kxut[i*m+j]=0;
              for(k=0;k<n2;k++) {
                kxut[i*m+j]+=kx[i*n2+k]*ut[k*m+j];

              }
         }
   }

//reform(kxut,nx,ny)

  double *reformkxut;
  reformkxut=(double*)fftw_malloc( sizeof (double)*nx*ny);
  e=0;
  for(p=0;p<nx;p++)
  {
    for(q=0;q<ny;q++)
    {
      reformkxut[p*ny+q]=kxut[e];
      e++;
    }
  }

// subtracting sfit from input

  for(i=0;i<nx;i++)
  {
     for(j=0;j<ny;j++)
     {
         in[i*ny+j][0]=in[i*ny+j][0]-reformkxut[i*ny+j];
     }
  }

  return in; 
}
