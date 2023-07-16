/*
NAME:           complexfftw
PURPOSE:        To return an fftw_complex array of input on which Fourier transform is performed
SAMPLE CALLING SEQUENCE:
                result =complexfftw(in,sign) 
INPUTS:         in:fftw_compex array of input values
                sign :
                   if(sign==1)
                      performes forward fourier transform
                   if(sign==-1)
                      performes reverse fourier transform
OUTPUTS:        res:fftw_complex array of manipulated input


RETURNS:        fftw_complex* in= fftw_array of input on which Fourier transform is performed


HISTORY:       Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/


#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include <stdlib.h>
#include<complex.h>
#define Pi 3.14159
extern int nx,ny;      //accessing global variable declared in fccorr.c
fftw_complex* complexfftw(fftw_complex *in,int sign)
{
  int i;
  int j;
  fftw_complex *out;
  fftw_plan plan;
  out = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
  if(sign==1)
  plan= fftw_plan_dft_2d ( nx, ny, in, out, FFTW_FORWARD, 
    FFTW_ESTIMATE );
  if(sign==-1)
    plan= fftw_plan_dft_2d ( nx, ny, in, out, FFTW_BACKWARD,
    FFTW_ESTIMATE );
  fftw_execute ( plan);
  fftw_destroy_plan ( plan);
  return out;
}

/*
NAME:           cor
PURPOSE:        Perform a cross correlation match on two images using Fourier transform
SAMPLE CALLING SEQUENCE:
                result = fccorr(img1,img2) 
INPUTS:        img1: 2-D image
               img2: 2-D image
               
OUTPUTS:       res: peak correlation coefficient


RETURNS:       res[3]:  res(0) = x-shift in pixels FROM img1 to img2.
                        res(1) = y-shift in pixels FROM img1 to img2.
                        res(2) = peak cross correlation co-efficient.


HISTORY:       Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/

float *cor(fftw_complex *in1,fftw_complex *in2)
{
  int i,j,k,x,y;
  float *res;
  res=fftw_malloc(sizeof(double*)*3);
  float max;
  fftw_complex *out1,*out2,*out3;
  out1 = complexfftw(in1,1);
  out2 = complexfftw(in2,1);
  out3 = fftw_malloc(sizeof(fftw_complex)*nx*ny);
 
//multipling FFT[in1]*conj(FFT[in2])

  for(i=0;i<nx;i++)
   {
     for(j=0;j<ny;j++)
      {
	
      		 out3[i*ny+j][0]=(out1[i*ny+j][0]*out2[i*ny+j][0])-(out1[i*ny+j][1]*(-out2[i*ny+j][1]));
      		 out3[i*ny+j][1]=(out1[i*ny+j][0]*(-out2[i*ny+j][1]))+(out1[i*ny+j][1]*out2[i*ny+j][0]);
           
      }
  }

  int row=0,column=0;
//FFT-1[FFT[in1]*conj(FFT[in2])]

  out3=complexfftw(out3,-1);
//Finding the maximum present in the given Output array
  max=out3[0][0];
  for (i=0; i<nx; i++)
  {
     for (j=0; j<ny; j++)
     {
         if(out3[i*ny+j][0]>=max)
         {
             max=out3[i*nx+j][0];
             row=i;column=j;
         }

      }
  }
  printf("row=%d,column=%d\n",row,column);
//Finding the row shift and column shift

/*  if(row>(nx/2))
     row=row-nx;
    if(column>(ny/2))
     column=column-ny;
   
    res[0]=column;
    res[1]=row;
    res[2]=max;

*/

//peak interpolation

  float dr,cm,rm;
  if(max!=1)
  {
    if(((row*column)!=0)&&(row<(nx-1))&&(column<(ny-1)))
    {
       dr=max*2-out3[row*ny+(column-1)][0]-out3[row*ny+(column+1)][0];
       cm=(column-0.5)+(max-out3[row*ny+(column-1)][0])/dr;
       dr=max*2-out3[(row-1)*ny+column][0]-out3[(row+1)*ny+column][0];
       rm=(row-0.5)+(max-out3[(row-1)*ny+column][0])/dr;
    }
    else
    {
        cm=column;
        rm=row;
    }
    res[0]=cm-ny/2;
    res[1]=rm-nx/2;
   }
   else
   {
      cm=column;
      rm=row;
      res[0]=cm-ny/2;
      res[1]=rm-nx/2;
   }
  res[2]=max;

  return res;

}

