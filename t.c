#include<stdio.h>
#include<math.h>
int main()
{
 printf("%lf\n",0.54-0.46*cos(2*3.14*2/4));
 printf("%lf\n",0.54-0.46*cos(2*3.14*(2.0/4.0)));
 printf("%lf\n",2*3.14*2/4);
 printf("%lf\n",2*3.14*(2/4));
 printf("%lf\n",2.0/4.0);
 printf("%lf\n",0.54-0.46*cos(2*3.14));
 return 0;
}
