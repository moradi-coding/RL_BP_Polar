
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

string func(long n) {
	stringstream result;
	result << n;
	return result.str();
}


int main() { 

	int i,i2,j,len_vndeg,len_cndeg,iter,iter_max=100,flg,flg2,cnt,opt=0,num,m=134; //opt=1 does VN degree optimization
	float R,a,b,mu_CN=0, mu_0,*mu_vn,*mu_cn,*lambda_prime,sum,sum2=0,phi_mu_vni,tmp,Perr,x,Qx,inc,sigma,sigma_old,tR,
	sum_j1_to_i,epsilon=0.6,beta=0.9,beta_i, Q[iter_max], EbNo_BC[6], rwrd[iter_max];

	ofstream outf1;
	
	//int vndeg[]={1,2,3}, cndeg[]={4}; tR=0.571; //target rate [7,4] Hamming
	//float lambda[]={0.429,0.429,0.142}, rho[]={1}, Q[iter_max]; 
	
	//int vndeg[]={2}, cndeg[]={7}; tR=0.714; //(2,7) LDPC code
	//float lambda[]={1}, rho[]={1}; 
	
	int vndeg[]={2,3,4,5,6}, cndeg[]={4}; tR=0.143;
	float lambda[]={0.257996,0.285714,0.264392,0.153518,0.0383795}, rho[]={1}; 

	len_vndeg=sizeof(vndeg)/sizeof(vndeg[0]);
	len_cndeg=sizeof(cndeg)/sizeof(cndeg[0]);

	mu_vn=new float[len_vndeg];
	mu_cn=new float[len_cndeg];
	lambda_prime=new float[len_vndeg];
	
	a=0;
	for(i=0;i<len_cndeg;i++) 
		a+=rho[i]*cndeg[i];
	//lambda[0]=0.05; lambda[1]=0.95;
	EbNo_BC[0]=1; EbNo_BC[1]=2; EbNo_BC[2]=3; EbNo_BC[3]=4; EbNo_BC[4]=4.5; EbNo_BC[5]=5;
	num=sizeof(EbNo_BC)/sizeof(EbNo_BC[0]);
	cout<<"num: "<<num<<'\n';
	
	for(i2=0;i2<num;i2++) {
		flg=1;
		sigma=sqrt(1/(2*tR*pow(10,0.1*EbNo_BC[i2])));
		cout<<'\n'<<"sigma: "<<sigma<<endl;
		sum2=0;
		//while(flg /*&& ((lambda[2]>0 && len_vndeg==3) || len_vndeg==2)*/) {
		mu_0=2/(sigma*sigma);
		mu_CN=0;
		for(iter=0; iter<iter_max; iter++) { //no. of decoder iters.
			sum=0;
			for(i=0;i<len_vndeg;i++) {
				mu_vn[i]=mu_0+(vndeg[i]-1)*mu_CN; //mean of VN to CN message	
				//cout<<"mu_vn: "<<mu_vn[i]<<'\n';
			
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
				tmp=1-pow(1-sum,cndeg[j]-1);
				if(tmp<0.00001)
					tmp=0.00001;
				mu_cn[j]=pow((0.0218-log(tmp))/0.4527,1/0.86); 
			}
				
			mu_CN=0;
			for(j=0;j<len_cndeg;j++)
				mu_CN+=rho[j]*mu_cn[j];
			//cout<<"mu_CN: "<<mu_CN<<'\n';			
			
			/*sum=0;
			for(i=0;i<len_vndeg;i++) 
				sum+=lambda[i]/vndeg[i];
			for(i=0;i<len_vndeg;i++) {
				lambda_prime[i]=(lambda[i]/vndeg[i])/sum;
				cout<<"lambda_prime: "<<lambda_prime[i]<<'\n';
			}*/	
		
			Perr=0;
			for(i=0;i<len_vndeg;i++) {
				x=sqrt((mu_0+vndeg[i]*mu_CN)/2);
				//Perr+= lambda_prime[i]*Q(x);
				Qx=((1-exp(-1.98*x))*exp(-0.5*x*x))/(1.135*sqrt(2*pi)*x); //according to wiki
				//Perr+= lambda_prime[i]*Qx;
				Perr+= lambda[i]*Qx;
			}
			//cout<<"Perr: "<<Perr<<endl;
				
			//Q-value calculation
			rwrd[iter]=1-Perr;
			/*sum_j1_to_i=0;
			for(i=0;i<iter;i++)
				sum_j1_to_i+=(epsilon/m)*(m-1)*rwrd+(1-epsilon+epsilon/m)*rwrd;
			beta_i=pow(beta,iter);
			if(iter<iter_max-1)
				sum2+=beta_i*sum_j1_to_i;
				
			Q[iter]=rwrd+(1-beta)*sum2;*/
				
			if(Perr<=pow(10,-12)) {
				flg=1;
				break;
			}
			else
				flg=0;
		}
		/*b=0;
		for(i=0;i<len_vndeg;i++) 
			b+=lambda[i]*vndeg[i];
		flg2=0;
		R=1-b/a;
		cout<<"a, b, R: "<<a<<", "<<b<<", "<<R<<'\n';*/
			//cout<<"rate mismatch: "<<abs(R-tR)<<'\n';	
		//}
		//if(flg2) break;	
		
		//Q-value calculation
		for(int lmax=0;lmax<iter_max;lmax++) {
			sum2=0; 
			for(i=1;i<lmax;i++) {
				sum_j1_to_i=0;
				for(j=1;j<=i;j++)
					sum_j1_to_i+=(epsilon/m)*(m-1)*rwrd[j]+(1-epsilon+epsilon/m)*rwrd[j];	
					sum_j1_to_i*=pow(beta,i);
				sum2+=sum_j1_to_i;
			}	
			Q[lmax]=rwrd[0]+(1-beta)*sum2;
		}
		cout<<'\n'<<"rwrd: "; for(i=0;i<iter_max;i++) cout<<rwrd[i]<<" ";
		//cout<<'\n'<<"Q: "; for(i=0;i<iter_max;i++) cout<<Q[i]<<" ";
		outf1.open("QDE"+func(i2)+".txt"); //, std::ios_base::app
		for(i=0;i<iter_max;i++) 
			outf1<<Q[i]<< " "; 
		outf1.close();
	}
}

