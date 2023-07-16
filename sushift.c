#include<stdio.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include <stdlib.h>
#include<complex.h>
#define Pi 3.14159
int nx=4,ny=4;
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

fftw_complex *shift(fftw_complex *img,int siftx,int sifty)
{
   int i,j,k;
   float temp;
   //xshift
   for(k=0;k<siftx;k++)
   {
      for(i=0;i<nx;i++)
      {
          temp=img[i*ny+(ny-1)][0];
          for(j=ny-1;j>0;j--)
                img[i*ny+j][0]=img[i*ny+(j-1)][0];
          img[i*ny+0][0]=temp;
      }
   }
   //column shift
   for(k=0;k<sifty;k++)
   {
      for(i=0;i<ny;i++)
      {
            temp=img[(nx-1)*ny+i][0];
            for(j=nx-1;j>0;j--)
               img[j*ny+i][0]=img[(j-1)*ny+i][0];
             img[j*ny+i][0]=temp;
      }
   }
  return img;
}
fftw_complex* sushift(fftw_complex *img1,int siftx,int sifty)
{
  int i,j,k;
  img1=shift(img1,nx/2,ny/2);
  img1= complexfftw(img1,-1);
  img1=shift(img1,nx/2,ny/2);
  fftw_complex *fimg,*out,*expo;
  double img,real;
  fimg=fftw_malloc(sizeof(fftw_complex)*nx*ny);
  out=fftw_malloc(sizeof(fftw_complex)*nx*ny);
  expo=fftw_malloc(sizeof(fftw_complex)*nx*ny);
  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
       real=0*2*Pi*((i-(nx/2))*(-siftx))/nx+((j-(ny/2))*(-sifty))/ny;;
       img=1*2*Pi*((i-(nx/2))*(-siftx))/nx+((j-(ny/2))*(-sifty))/ny;
//       printf("%f\n",img);
       expo[i*ny+j][0]=exp(real)*cos(img);
       expo[i*ny+j][1]=exp(real)*sin(img);  
       fimg[i*ny+j][0]=img1[i*ny+j][0]*expo[i*ny+j][0]+img1[i*ny+j][1]*expo[i*ny+j][1];
       fimg[i*ny+j][1]=img1[i*ny+j][0]*expo[i*ny+j][1]+img1[i*ny+j][1]*expo[i*ny+j][0];
 //      printf("%f %f\n",fimg[i*ny+j][0],fimg[i*ny+j][1]);
    }
  }
  out=shift(fimg,nx/2,ny/2);
  out=complexfftw(out,1);
  out=shift(out,nx/2,ny/2);
  return out;
  return img1;
}
int main()
{
  int i,j,k,n=1,p;
  fftw_complex *in1,*in2,*out;
  int num1;
  FILE *fp,*fptr,*fi1,*fi2;
  in1 = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
  out = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );

    fi1=fopen("ip2.dat","r");          //ip1 file reading

    for ( i = 0; i < nx; i++ )
    {
      for ( j = 0; j < ny; j++ )
      {
        p= fscanf(fi1,"%d",&num1);
        in1[i*ny+j][0] =num1;
        in1[i*ny+j][1] =0;
     }
   }
   out=sushift(in1,1,1);
   fclose(fi1);
 for(i=0;i<nx;i++)
 {
   for(j=0;j<ny;j++)
   {
     printf("%f %f\n",out[i*ny+j][0],out[i*ny+j][1]);
   }
 } 
  fftw_free ( in1 );
  return 0;
}

