#include<stdio.h>
int main()
{
   int arr[4][3]={2,3,4,
                  5,6,7,
                  8,9,0,
                 1,6,3};

  int i,j,nx=4,ny=3,k;
  int sx,sy;
  int temp;
  printf("enter x shift and y shift");
  scanf("%d %d",&sx,&sy);
 for(i=0;i<nx;i++)
 {
   for(j=0;j<ny;j++)
      printf("%d\t",arr[i][j]);
   printf("\n");
}
//row shift

for(k=0;k<sx;k++)
{
for(i=0;i<nx;i++)
{
  temp=arr[i][ny-1];
  for(j=ny-1;j>0;j--)
     arr[i][j]=arr[i][j-1];
  arr[i][0]=temp;
}
}

//column shift
for(k=0;k<sy;k++)
{
for(i=0;i<ny;i++)
{
  temp=arr[nx-1][i];
  for(j=nx-1;j>0;j--)
     arr[j][i]=arr[j-1][i];
  arr[j][i]=temp;
}
}

/*   int *flat=&arr;
   for(k=0;k<sx;k++)
   {
    temp=flat[0];
    for(i=1;i<nx*ny;i++)
      flat[i-1]=flat[i];

    flat[nx*ny-1]=temp;
   }
 */
 for(i=0;i<nx;i++)
 {
   for(j=0;j<ny;j++)
     printf("%d\t",arr[i][j]);
   printf("\n");
}
}
