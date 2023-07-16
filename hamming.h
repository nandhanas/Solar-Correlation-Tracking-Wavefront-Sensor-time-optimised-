/*
NAME:           hamming
PURPOSE:        To return an fftw_complex array of input multiplied with hamming function
SAMPLE CALLING SEQUENCE:
                result =hamming(in,h) 
INPUTS:         in:fftw_compex array of input values
                h :floating point array of hamming function
OUTPUTS:        res:fftw_complex array of manipulated input


RETURNS:        fftw_complex* in= fftw_array of input multiplied with hamming function


HISTORY:       Written by Nandhana in C language to improve the execution time of a function already present in IDL in May 2017
*/


#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include <stdlib.h>
#include<complex.h>
#define PI 3.14159265358979323846

extern int nx,ny;          //accessing global varible declared in fccorr.c

fftw_complex* hamming(fftw_complex* in,float *h)
{
   int i,j;
   for(i=0;i<nx;i++)
   {
     for(j=0;j<ny;j++)
     {
        in[i*ny+j][0]*=h[i*ny+j];
        in[i*ny+j][1]*=h[i*ny+j];
     }
  }
  return in;
}
