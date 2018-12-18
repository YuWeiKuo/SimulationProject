#include<stdio.h>
#include<math.h>
#include<stdlib.h>


double sf(double, double, double);
float bessj0(float);
int main()
{
	
	FILE *f;
	long i,j,t;
	double g[6000],r[6000],L,S[30000],K[30000];
	int a,conf;
	double V,N;
	double density;
	double p,q,l,m,s,k,k_min,k_max,dK;


  	k_min=0.01;
	k_max=10*M_PI,dK;

	dK = 0.001;

	double R,Lmax,Lmin,H;

	N = 242.131117;      /*******************/
	R = 46.859999999999999;
	L = 44.271219000000002;
	Lmax = 3.5;
	Lmin = 0.0;
	H = Lmax - Lmin;
	density = N/(R*L*(Lmax-Lmin));
	double final_term;
	double first_term;

	f = fopen("gxyOO_7.dat","r");
	
	conf = 1;

	for(i = 0;i < 30000;i++) S[i] = 0.0;
	for(i = 0;i < 30000;i++) K[i] = 0.0;

	for(int t =0;t < conf;t++){
			for(i = 0;i < 6000;i++){

				a=fscanf(f,"%lf %lf",&r[i],&g[i]);

			}
			

	
			k=0.0001;
			for(j = 0;j < 30000;j++){
				s = 0.0;
				K[j] = k;
	
				first_term=0.0;
				final_term=0.0;
				for(i = 0;i < 5999;i++){
					if(r[i]<0.5*L){
						if(i==0){
							first_term=first_term+sf(g[0],k,r[0]);
						}
						s = s + sf(g[i],k,r[i]);
						if(abs(r[i]-0.5*L)<=0.001){
							final_term=final_term+sf(g[i],k,r[i]);
						}
					}
					/*if(abs(r[i]-0.5*L)<=0.00001){
						final_term=final_term+sf(g[i],k,r[i]);
					}*/
					else continue;
					
		
		
				}
				
				S[j] = S[j] + 1.0 + 2.0*density*M_PI*0.5*(r[1]-r[0])*(first_term+final_term+2.0*s);
				k = k + dK;
			}	


	}



	for(i = 0;i < 30000;i++)
			printf("%lf %lf\n",K[i],S[i]/conf);


}

double sf(double g, double k, double r)
{
	double x;
	x = k*r;
	return (g-1.00)*bessj0(x)*r;

}


float bessj0(float x)
{
float ax,z;
double xx,y,ans,ans1,ans2;
if ((ax=fabs(x)) < 8.0) {
y=x*x;
ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
+y*(59272.64853+y*(267.8532712+y*1.0))));
ans=ans1/ans2;
} else {
z=8.0/ax;
y=z*z;
xx=ax-0.785398164;
ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
+y*(-0.2073370639e-5+y*0.2093887211e-6)));
ans2 = -0.1562499995e-1+y*(0.1430488765e-3
+y*(-0.6911147651e-5+y*(0.7621095161e-6
-y*0.934945152e-7)));
ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
}
return ans;
}
