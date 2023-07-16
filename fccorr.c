/*
NAME:           fccorr
PURPOSE:	Perform a cross correlation match on two images using Fourier transform
SAMPLE CALLING SEQUENCE:
		result = fccorr(img1,img2,flag) 
INPUTS:        img1: 2-D image
               img2: 2-D image
               flag:
                    if(flag==1)performs surfacefit
                    else       doesn't perform surfacefit
OUTPUTS:       res: peak correlation coefficient


RETURNS:      float* res[3]:  res(0) = x-shift in pixels FROM img1 to img2. (converted into int)
                              res(1) = y-shift in pixels FROM img1 to img2. (converted into int)
                              res(2) = peak cross correlation co-efficient.


HISTORY:       Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/
#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include<stdlib.h>
#include<complex.h>
#include "preprocess.h"   //file which contains functions that returns constants which are same for a given dimension
#include "hamming.h"      //file which contains hamming function which manipulates input array according to the hamming function
#include "correlation.h" // file which contains functions which performs cross correlation process and retunrs the row,column and max(op array)
#include "surfacefit.h" //file which contain function that performs surface fit

int  nx=32,ny=32;   	//global variable for dimensions of the input array
int *d,*f;  	        //global array for dark image and gain table correction
double *kk; 	       //global array for surface fit constant
float *h;             //global array for hamming function
int n2,m;             //global variable for surfacefit function

float* fccorr(fftw_complex *in1,fftw_complex *in2,int flag)
{

  int i,j,k; 
  int nip;
  nip=nx*ny;
  double sum1=0,sum2=0;
  float *res;
  for ( i = 0; i < nx; i++ )
   {
     for ( j = 0; j < ny; j++ )
     {
     //dark image subtraction and gain table correction
       in1[i*ny+j][0] =(in1[i*ny+j][0]-d[i*ny+j])/f[i*ny+j]; 
       in2[i*ny+j][0] =(in2[i*ny+j][0]-d[i*ny+j])/f[i*ny+j];
   

     //sum of all points in an image
       sum1+=in1[i*ny+j][0];   
       sum2+=in2[i*ny+j][0];

    } 
  }

  //mean
  sum1=sum1/nip;
  sum2=sum2/nip;
  
 //subtracting mean from the image
 for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      in1[i*ny+j][0]=in1[i*ny+j][0]-sum1;
      in2[i*ny+j][0]=in2[i*ny+j][0]-sum2;
    }
  }
 

 //bilinear fitting
 if(flag==1)
 {
 in1=surfacefit(in1,kk);
 in2=surfacefit(in2,kk);
 }
 
 //hamming window
  in1=hamming(in1,h);
  in2=hamming(in2,h);
  

 for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
    //adding mean to the image
      in1[i*ny+j][0]=in1[i*ny+j][0]+sum1;
      in2[i*ny+j][0]=in2[i*ny+j][0]+sum2;

    //dividing the image by mean
       in1[i*ny+j][0]/=sum1;
       in2[i*ny+j][0]/=sum2;
    }
  }
  res=cor(in1,in2);

  return res;
}
int main()
{
  int i,j,k,n=1,p;
  fftw_complex *in1,*in2,*out;
  int num1,num2;
  float *res;
  d =imgprocess_darkimage(); //dark image constant
  f =imgprocess_gaintable(); //gain table correction constant
  h= hamming_function();     //hamming function constant
  kk=sfit_kk();             //surfacefit constant
  unsigned int seed = 123456789;
  FILE *fp,*fptr,*fi1,*fi2;
  fptr = fopen("exetimeadaptive8.dat", "w");
  fp=fopen("datapeakinter.dat","a");
  in1 = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
  in2 = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
  out = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
  srand ( seed );
  for(k=0;k<n;k++)
  {
    fi1=fopen("img8.dat","r");          //ip1 file reading
    fi2=fopen("img4.dat","r");          //ip2 file reading

    for ( i = 0; i < nx; i++ )
    {
      for ( j = 0; j < ny; j++ )
      {
        p= fscanf(fi1,"%d",&num1);
        p= fscanf(fi2,"%d",&num2);
        in1[i*ny+j][0] =num1;
        in1[i*ny+j][1] = 0;
        in2[i*ny+j][0] =num2;
        in2[i*ny+j][1] = 0;
     }
   }
 //  for(k=0;k<n;k++)
 //  {
   clock_t t;
   t = clock();
   res=fccorr(in1,in2,0);
   t = clock() - t;
   double time_taken = ((double)t)/CLOCKS_PER_SEC; 
   fprintf(fptr,"%f\n", time_taken);
  // fprintf(fp,"Testcase=%d\n",k);
   fprintf(fp,"img8&img4\n\n");
   fprintf(fp,"column=%f\nrow=%f\nmax=%f\n\n",res[0],res[1],res[2]);
 //  }
   fclose(fi2);
   fclose(fi1);
  } 
  fftw_free ( in1 );
  fftw_free (in2);
  fclose(fptr);
  fclose(fp);
  return 0;
}
