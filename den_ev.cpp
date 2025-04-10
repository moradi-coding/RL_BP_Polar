
//density evolution for regular/irregular codes based on Gaussian approximation

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <cmath>

using namespace std;

#define pi 3.14159265358979323846

int main() { 

	int i,i2,j,len_vndeg,len_cndeg,iter,iter_max=100,flg,flg2,cnt,opt=0,num; //opt 1 does VN degree optimization
	float R,a,b,mu_CN=0, mu_0,*mu_vn,*mu_cn,*lambda_prime,sum,phi_mu_vni,tmp,Perr,x,Qx,inc,sigma,sigma_old,tR;

	//int vndeg[]={2,3,6,7,20},cndeg[]={8,9};
	//float lambda[]={0.234,0.212,0.147,0.103,0.304}, rho[]={0.719,0.281}; //fraction of VNs and CNs with those degrees
	
	int vndeg[]={1,2,3}, cndeg[]={4}; tR=0.571; //target rate [7,4] Hamming
	float lambda[]={0.429,0.429,0.142}, rho[]={1}; 
	
	//int vndeg[]={3}, cndeg[]={6}; tR=0.5; //target rate (3,6) LDPC
	//float lambda[]={1}, rho[]={1}; 

	//int vndeg[]={4,5}; float lambda[]={0,0};
	//int cndeg[]={60,61,62,63}; float rho[]={0.143,0.476,0.33,0.0477};

	float **lambda_sv; //let the program find the optimal VN deg. dist. (if opt=1)
	

	len_vndeg=sizeof(vndeg)/sizeof(vndeg[0]);
	len_cndeg=sizeof(cndeg)/sizeof(cndeg[0]);

	mu_vn=new float[len_vndeg];
	mu_cn=new float[len_cndeg];
	lambda_prime=new float[len_vndeg];
	
	a=0;
	for(i=0;i<len_cndeg;i++) 
		a+=rho[i]*cndeg[i];

	//cnt=0;
	inc=0.05;
	//num=1/inc-1;
	if(!opt) 
		num=1; //number of different VN degree dist.
	else
		num=20;
	lambda_sv=new float*[3]; 
	for(i=0;i<3;i++) 
		lambda_sv[i]=new float[num];
	cout<<"num: "<<num<<'\n';
	//lambda[0]=0.05; lambda[1]=0.95;
	
	for(i2=0;i2<num;i2++) {
		flg=1;
		if(opt) {
			lambda[0]+=inc; 
			if(len_vndeg==2) 
				lambda[1]=1-lambda[0];
			else {
				lambda[1]+=2*inc; 
				lambda[2]=1-(lambda[0]+lambda[1]);
			}		
			//cout<<"lambdas: "<<lambda[0]<<" "<<lambda[1]<<'\n';
			cout<<"lambda: "<<lambda[0]<<" "<<lambda[1]<<" "<<lambda[2]<<'\n';
		}
		
		sigma_old=0;
		//sigma=1.67;
		sigma=1.05; 
		while(flg /*&& ((lambda[2]>0 && len_vndeg==3) || len_vndeg==2)*/) {
			mu_0=2/(sigma*sigma);
			mu_CN=0;
			for(iter=0; iter<iter_max; iter++) { //no. of decoder iters.
				for(i=0;i<len_vndeg;i++) {
					mu_vn[i]=mu_0+(vndeg[i]-1)*mu_CN;	
					//cout<<"mu_vn: "<<mu_vn[i]<<'\n';
				}
			
				sum=0;
				for(i=0;i<len_vndeg;i++) {
					//if(mu_vn[i]>=0 && mu_vn[i]<=10)
						phi_mu_vni=exp(-0.4527*pow(mu_vn[i],0.86)+0.0218);
					//else
						//phi_mu_vni = sqrt(pi/mu_vn[i])*exp(-mu_vn[i]/4)*(1-10/7/mu_vn[i]);
					//cout<<"phi_mu_vni: "<<phi_mu_vni<<'\n';
					sum+=lambda[i]*phi_mu_vni;
				}
				//cout<<"sum: "<<sum<<'\n';
			
				for(j=0;j<len_cndeg;j++) {
					//mu_cn[j]=phi_inv(1-pow(1-sum,cndeg[j]-1));
					tmp= 1-pow(1-sum,cndeg[j]-1);
					if(tmp<0.00001)
						tmp=0.00001;
					mu_cn[j]=pow((0.0218-log(tmp))/0.4527,1/0.86); 
				}
				
				mu_CN=0;
				for(j=0;j<len_cndeg;j++)
					mu_CN+=rho[j]*mu_cn[j];
					//cout<<"mu_CN: "<<mu_CN<<'\n';			
			
				sum=0;
				for(i=0;i<len_vndeg;i++) 
					sum+=lambda[i]/vndeg[i];
				for(i=0;i<len_vndeg;i++) {
					lambda_prime[i]=(lambda[i]/vndeg[i])/sum;
					//cout<<"lambda_prime: "<<lambda_prime[i]<<'\n';
				}	
		
				Perr=0;
				for(i=0;i<len_vndeg;i++) {
					x=sqrt((mu_0+vndeg[i]*mu_CN)/2);
					//Perr+= lambda_prime[i]*Q(x);
					Qx=((1-exp(-1.98*x))*exp(-0.5*x*x))/(1.135*sqrt(2*pi)*x); //according to wiki
					Perr+= lambda_prime[i]*Qx;
				}
				cout<<"Perr: "<<Perr<<endl;
				if(Perr<=pow(10,-9)) {
					flg=1;
					break;
				}
				else
					flg=0;
			}
			b=0;
			for(i=0;i<len_vndeg;i++) 
				b+=lambda[i]*vndeg[i];
			flg2=0;
			R=1-b/a;
			cout<<"a, b, R: "<<a<<", "<<b<<", "<<R<<'\n';
			//cout<<"rate mismatch: "<<abs(R-tR)<<'\n';
				
			if(flg && sigma>sigma_old && abs(R-tR)<0.001) {
				sigma_old=sigma;
				sigma+=0.01; //keep incrementing until threshold is reached
				if(opt) {
					lambda_sv[0][i2]=lambda[0];
					lambda_sv[1][i2]=lambda[1];
					if(len_vndeg==3) 
						lambda_sv[2][i2]=lambda[2];
				}
				/*if(len_vndeg==2) 
					cout<<"lambda[0], lambda[1]: "<<lambda_sv[0][i2]<<" "<<lambda_sv[1][i2]<<'\n';
				else
					cout<<"lambda 0, 1, 2: "<<lambda_sv[0][i2]<<" "<<lambda_sv[1][i2]<<" "<<lambda_sv[2][i2]<<'\n';*/
				cout<<"sigma: "<<sigma<<'\n';
			}
			else {
				//flg2=1;
				break;
			}
		}
		
		//if(flg2) break;	
	}
}

