/*
NAME:           hamming_function
PURPOSE:        To return the hamming function
SAMPLE CALLING SEQUENCE:
                result = hamming_function() 
INPUTS:         NULL
OUTPUTS:       res: peak correlation coefficient


RETURNS:       float *h1: a floating point array which consists of hamming function


HISTORY:       Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/

#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include <stdlib.h>
#include<complex.h>
#define PI 3.14159265358979323846
#define nz 1
#define degree 1

double *ut;          //global array declaration to use it in surfacefit function
extern int m,n2;    //Using global variable declared in fccorr.c
extern int nx,ny;   //Using global variable declared in fccorr.c

float* hamming_function()
{
  int i,j;
  float *a,*b;
  float *h1;
  a = fftw_malloc ( sizeof ( float ) * nx );
  b = fftw_malloc ( sizeof ( float ) * ny  );
  h1= fftw_malloc ( sizeof ( double ) * nx * ny );
  for(i=0;i<nx;i++)
  {
    a[i]=0.54+(0.54-1)*cos(2*PI*(i/(nx-1)));
    b[i]=0.54+(0.54-1)*cos(2*PI*(i/(ny-1)));
   
  }
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
      h1[i*ny+j]=a[i]*b[j];
    }
  }
  return h1;
}

/*
NAME:           imgprocess_darkimage
PURPOSE:        To return an int array of value 0 for dark image correction
SAMPLE CALLING SEQUENCE:
                result = imgprocess_darkimage() 
INPUTS:         NULL
OUTPUTS:        res:an int array of value 0


RETURNS:        int *d1= an int array of value 0


HISTORY:        Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/

int* imgprocess_darkimage()
{
   int *d1,i,j;
   d1 = fftw_malloc ( sizeof ( int ) * nx * ny );
   for ( i = 0; i < nx; i++ )
   {
      for ( j = 0; j < ny; j++ )
      {
         d1[i*ny+j] = 0;
      }
   }
   return d1;
}

/*
NAME:           imgprocess_gaintable
PURPOSE:        To return an int array of value 1 for gain table correction
SAMPLE CALLING SEQUENCE:
                result = imgprocess_gaintable() 
INPUTS:         NULL
OUTPUTS:        res:an int array of value 1


RETURNS:        int *f1= an int array of value 1


HISTORY:        Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/

int* imgprocess_gaintable()
{
   int *f1,i,j;
   f1 = fftw_malloc ( sizeof ( int ) * nx * ny );
   for ( i = 0; i < nx; i++ )
   {
      for ( j = 0; j < ny; j++ )
      {
         f1[i*ny+j] = 1;
      }
   }
   return f1;
}

/*
NAME:           sfit_kk
PURPOSE:        TO return a double array of surface fit constant
SAMPLE CALLING SEQUENCE:
                result = sfit_kk()() 
INPUTS:         NULL
OUTPUTS:        res:a double array of surface fit constant


RETURNS:        double *kk= a double array of surface fit constant


HISTORY:        Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/

double* sfit_kk()
{
   int i,j,k;
   n2=pow((degree+1),2);
   m=nx*ny;
   float *x,*y;
   x = fftw_malloc ( sizeof ( float ) * nx * ny );
   y = fftw_malloc ( sizeof ( float ) * nx * ny );
   float *findgenx,*repx,*findgeny,*repy;
   findgenx =fftw_malloc (sizeof (float) *nx);
   repx =fftw_malloc (sizeof (float) *ny);
   findgeny =fftw_malloc (sizeof (float) *ny);
   repy =fftw_malloc (sizeof (float) *nx);
   for(i=0;i<nx;i++)
   {
      repy[i]=1;
      findgenx[i]=i;
   }
   for(i=0;i<ny;i++)
   {
      repx[i]=1;
      findgeny[i]=i;
   }

//# function (columnwise multiplication in IDL)

   for(i=0;i<nx;i++) {
        for(j=0;j<ny;j++) {
            x[j*nx+i]=0;
            y[j*nx+i]=0;
            for(k=0;k<nz;k++) {
              x[j*nx+i]+=findgenx[k*nx+i]*repx[j*nz+k];
               y[j*nx+i]+=repy[k*nx+i]*findgeny[j*nz+k];

            }
        }
    }


//

    double *xy;
    ut=(double*)fftw_malloc( sizeof (double)*n2*m);
    xy=fftw_malloc( sizeof (double) * nx * ny);
    int kl=0,e=0;
    int p,q;
    double *addr;
    
//finding ut[kl,0]=reform(x^i,y^j,1,m)  
  
    for(i=0;i<=degree;i++)
    {
      for(j=0;j<=degree;j++)
      {
        for(p=0;p<nx;p++)
        { 
          for(q=0;q<ny;q++)
          {
             xy[p*ny+q]=pow(x[p*ny+q],i) * pow(y[p*ny+q],j);
          }
        }

 //ut reforming

        e=0;
        for(p=0;p<nx;p++)
        {
         for(q=0;q<ny;q++)
         {
           ut[kl*m+e]=xy[p*ny+q];
           e++;
         }
        }
        kl++;     
      }
    }


//kk = invert(ut # transpose(ut)) # ut
//transpose of ut

   double *transut;
   transut=(double*)fftw_malloc( sizeof (double)*n2*m);
   for(i=0; i<n2; ++i)
        for(j=0; j<m; ++j)
        	{
           	 transut[j*n2+i] = ut[i*m+j];
        	}


//column order multiply of ut and transut

   double *mulut;
   float matrix[10][10];
   mulut=(double*)fftw_malloc( sizeof (double)*n2*n2);
   for(i=0;i<n2;i++) {
        for(j=0;j<n2;j++) {
            mulut[i*n2+j]=0;
            for(k=0;k<m;k++) {
              mulut[i*n2+j]+=ut[i*m+k]*transut[k*n2+j];
              matrix[i][j]=mulut[i*n2+j];
            }
        }
    }

//inverse of mulut

   int n=n2;
   double a,ratio;
  for(i = 0; i < n; i++){
     for(j = n; j < 2*n; j++){
        if(i==(j-n))
             matrix[i][j] = 1.0;
               else
             matrix[i][j] = 0.0;
         }
      }
  for(i = 0; i < n; i++){
     for(j = 0; j < n; j++){
        if(i!=j){
              ratio = matrix[j][i]/matrix[i][i];
          for(k = 0; k < 2*n; k++){
              matrix[j][k] -= ratio * matrix[i][k];
              }
          }
        }
    }
  for(i = 0; i < n; i++){
      a = matrix[i][i];
      for(j = 0; j < 2*n; j++){
             matrix[i][j] /= a;
      } 
  }

  for(i = 0; i < n; i++){
     for(j = n,k=0; j < 2*n,k<n; j++,k++){
          mulut[i*n+k]=matrix[i][j];
        }
  }


//# of inverse mulut and ut to get kk

  double *kk;
  kk=(double*)fftw_malloc( sizeof (double)*n2*m);

  for(i=0;i<n2;i++) {
        for(j=0;j<m;j++) {
            kk[i*m+j]=0;
            for(k=0;k<n2;k++) {
              kk[i*m+j]+=mulut[i*n2+k]*ut[k*m+j];

            }
        }
    }
  return kk;
}
