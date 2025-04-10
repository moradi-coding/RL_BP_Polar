
//density evolution based on Gaussian approximation for sequential decoding

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
	srand(time(0));	
	int i,a,k,i2,j,len_vndeg,len_cndeg,iter,iter_max=100,flg,flg2,cnt,opt=0,num,m=134,**zeta,num_zeta=50,fn; //opt=1 does VN degree optimization
	float R,b,*mu_CN,mu_0,***mu_vn,**mu_cn,sum,sum1,sum2=0,**phi_mu_vni,tmp,Perr,x,Qx,inc,sigma,sigma_old,
	**rwrd,sum_j1_to_i,epsilon=0.6,beta=0.9,beta_i,*sum_rwd,*V, **Q, EbNo_BC[8];

	ofstream outf1,outf2;
	//int vndeg[]={1,2,3}, cndeg[]={4}; //target rate [7,4] Hamming
	//float lambda[]={0.429,0.429,0.142}, rho[]={1}; 
	
	//int vndeg[]={2}, cndeg[]={7}; R=0.714;  //(2,7) LDPC code
	//float lambda[]={1}, rho[]={1}; fn=1;
		
	//int vndeg[]={2,3,4,5,6}, cndeg[]={4,7}; R=0.501066;
	//float lambda[]={0.637527,0.226013,0.10661,0.0255864,0.00426439}, rho[]={0.641026,0.358974}; fn=2;
		
	//int vndeg[]={2,3,4,5,6}, cndeg[]={4,7}; R=0.287847; //for fn=5
	//float lambda[]={0.351812,0.34968,0.200426,0.0767591,0.021322}, rho[]={0.898204,0.101796}; fn=5;
	
	//int vndeg[]={2,3,4,5,6}, cndeg[]={4}; R=0.151386;
	//float lambda[]={0.264392,0.283582,0.266525,0.151386,0.0341151}, rho[]={0.994975,0.00502513}; fn=7;
	
	int vndeg[]={2,3,4,5,6}, cndeg[]={4}; R=0.143;
	float lambda[]={0.257996,0.285714,0.264392,0.153518,0.0383795}, rho[]={1}; fn=8;

	len_vndeg=sizeof(vndeg)/sizeof(vndeg[0]);
	len_cndeg=sizeof(cndeg)/sizeof(cndeg[0]);
	
	/*a=b=0;
	for(i=0;i<len_cndeg;i++) 
		a+=rho[i]*cndeg[i];
	for(i=0;i<len_vndeg;i++) 
		b+=lambda[i]*vndeg[i];
	R=1-b/a;*/
	cout<<"R: "<<R<<'\n';
	
	//mu_vn=new float*[len_vndeg]; 
	//for(i=0;i<len_vndeg;i++) 
		//mu_vn[i]=new float[m];
		
	mu_vn = new float**[len_vndeg];
	for(i = 0; i < len_vndeg; i++) {
    	mu_vn[i] = new float*[m];
   		for(j = 0; j< m; j++) {
    	    mu_vn[i][j] = new float[num_zeta];
    	    for(int k = 0; k < num_zeta; k++) { 
    	        mu_vn[i][j][k] = 0;
    	    }
    	}
	}
		
	mu_cn=new float*[len_cndeg];
	for(i=0;i<len_cndeg;i++) 
		mu_cn[i]=new float[m];
		
	phi_mu_vni=new float*[len_vndeg];
	for(i=0;i<len_vndeg;i++) 
		phi_mu_vni[i]=new float[num_zeta];
		
	Q=new float*[iter_max]; 
	for(i=0;i<iter_max;i++) 
		Q[i]=new float[m];
		
	rwrd=new float*[iter_max]; 
	for(i=0;i<iter_max;i++) 
		rwrd[i]=new float[m];
		
	mu_CN=new float[m];
	sum_rwd=new float[iter_max];
	V=new float[iter_max];
	
	zeta=new int*[len_vndeg]; 
	for(i=0;i<len_vndeg;i++) 
		zeta[i]=new int[m];
	
	
	
	//lambda[0]=0.05; lambda[1]=0.95;
	EbNo_BC[0]=-1; EbNo_BC[1]=0; EbNo_BC[2]=1; EbNo_BC[3]=2; EbNo_BC[4]=3; EbNo_BC[5]=4; EbNo_BC[6]=5; EbNo_BC[7]=6;
	num=sizeof(EbNo_BC)/sizeof(EbNo_BC[0]);
	cout<<"num: "<<num<<'\n';
	
	for(i2=0;i2<num;i2++) { //number of different SNRs
		flg=1;
		sigma=sqrt(1/(2*R*pow(10,0.1*EbNo_BC[i2])));
		cout<<'\n'<<"sigma: "<<sigma<<endl;
		mu_0=2/(sigma*sigma);
		for(iter=0; iter<iter_max; iter++) { //no. of decoder iters.
			for(a=0;a<m;a++) { //a is index of scheduled CN
				cnt=0;				
				
				while(cnt<num_zeta) { //randomly generate num_zeta vectors for each VN degree 
					for(i=0;i<len_vndeg;i++) 
						for(j=0;j<m;j++) zeta[i][j]=0; //refresh
							
					for(i=0;i<len_vndeg;i++) {
						sum2=0;
						while(sum2<vndeg[i]) {
							j=rand()%m;
							//cout<<"j: "<<j<<'\n';
							if(!zeta[i][j]) {
								flg=0;
								while(!flg) {
									k=rand()%vndeg[i];
									if(k+sum2<=vndeg[i]) {
										flg=1;
										zeta[i][j]=k;
										sum2+=zeta[i][j];
									}
									else
										flg=0;
								}
							}
						}
					}
					/*cout<<'\n'<<"zeta: "; 
					for(i=0;i<len_vndeg;i++) {
						for(j=0;j<m;j++) 
							cout<<zeta[i][j]<<" ";
						cout<<'\n';
					}*/		
				
					for(i=0;i<len_vndeg;i++) {	
						if(zeta[i][a]) {
							sum1=sum2=0;
							for(j=0;j<a;j++)
								sum1+=zeta[i][j]*mu_CN[j]; //CN to VN mean of jth CN
							for(j=a+1;j<m;j++)
								sum2+=zeta[i][j]*mu_CN[j];
							//mu_vn[i]=mu_0+(vndeg[i]-1)*mu_CN; //mean of VN to CN message	
							mu_vn[i][a][cnt]=mu_0+sum1+(zeta[i][a]-1)*mu_CN[a]+sum2;
							//cout<<"mu_vn: "<<mu_vn[i][a]<<'\n';
							phi_mu_vni[i][cnt]=exp(-0.4527*pow(mu_vn[i][a][cnt],0.86)+0.0218); //used in CN to VN mean calculation					
							//sum+=lambda[i]*phi_mu_vni;
						}
					}
					cnt++;
					//cout<<"sum, a: "<<sum<<", "<<a<<endl;
				}
				sum=0;
				for(i=0;i<len_vndeg;i++) 
					for(j=0;j<num_zeta;j++)
						sum+=lambda[i]*phi_mu_vni[i][j];
				sum/=num_zeta;
				//cout<<"sum: "<<sum<<endl;
				
				if(sum) { //if the CN is connected to this VN
					for(j=0;j<len_cndeg;j++) {
						//mu_cn[j]=phi_inv(1-pow(1-sum,cndeg[j]-1));
						tmp=1-pow(1-sum,cndeg[j]-1);
						if(tmp<0.00001)
							tmp=0.00001;
						//cout<<"tmp: "<<tmp<<'\n';
						mu_cn[j][a]=pow((0.0218-log(tmp))/0.4527,1/0.86); //here log means ln	
						//cout<<"mu_CN, a : "<<mu_cn[j][a]<<", "<<a<<'\n';	
					}
					
					mu_CN[a]=0;
					for(j=0;j<len_cndeg;j++)
						mu_CN[a]+=rho[j]*mu_cn[j][a];
					//cout<<"mu_CN, a : "<<mu_CN[a]<<", "<<a<<endl;			
				
					Perr=0;
					for(i=0;i<len_vndeg;i++) {
						//x=sqrt((mu_0+vndeg[i]*mu_CN[a])/2); 
						Qx=0;
						for(j=0;j<num_zeta;j++) {
							if(!mu_vn[i][a][j])
								mu_vn[i][a][j]=pow(10,-9);
						
							//sum3/=cnt;
							x=sqrt(mu_vn[i][a][j]/2); 
							//x=sqrt(sum3/2)
							Qx+=((1-exp(-1.98*x))*exp(-0.5*x*x))/(1.135*sqrt(2*pi)*x); //according to wiki
						}
						//cout<<"x, Qx: "<<x<<", "<<Qx<<endl;
						Perr+= lambda[i]*Qx;
					}
					//cout<<"Perr: "<<Perr<<endl;
					//rwrd[iter][a]=1-Perr;
					rwrd[iter][a]=1-Perr/num_zeta;
				}
				else
					rwrd[iter][a]=0;
			}
			//cout<<'\n'<<"rwrd: "; for(i=0;i<iter_max;i++) cout<<rwrd[i][a]<<" ";
		}
					
		//Q-value calculation
		for(i=0;i<iter_max;i++)
			sum_rwd[i]=0; //refresh
		for(j=0;j<iter_max;j++)
			for(a=0;a<m-1;a++) //assuming CN indices 0,..,m-2 are in \mathcal{A}'
				sum_rwd[j]+=rwrd[j][a];
		//cout<<'\n'<<"sum_rwd: "; for(i=0;i<m-1;i++) cout<<sum_rwd[i]<<" ";
														
		for(a=0;a<m;a++) {
			for(int lmax=0;lmax<iter_max;lmax++) {
				sum2=0; 
				for(i=1;i<lmax;i++) {
					sum_j1_to_i=0;
					for(j=1;j<=i;j++) 
						sum_j1_to_i+=(epsilon/m)*sum_rwd[j]+(1-epsilon+epsilon/m)*rwrd[j][m-1];	
					sum_j1_to_i*=pow(beta,i);
					sum2+=sum_j1_to_i;
				}	
				Q[lmax][a]=rwrd[0][a]+(1-beta)*sum2;
			}
		}
		
		for(i=0;i<iter_max;i++)
			sum_rwd[i]=0; //refresh
		for(i=0;i<iter_max;i++)
			for(a=0;a<m-1;a++)
				sum_rwd[i]+=Q[i][a];
				
		for(i=0;i<iter_max;i++) 
			V[i]=(epsilon/m)*sum_rwd[i]+(1-epsilon+epsilon/m)*Q[i][m-1]; //value function
		
		cout<<'\n'<<"V: "; for(i=0;i<iter_max;i++) cout<<V[i]<<" ";
		outf1.open("V"+func(fn)+"_"+func(i2)+".txt"); //, std::ios_base::app
		for(i=0;i<iter_max;i++) 
			//outf1<<sum_rwd[i]<< " "; 
			outf1<<V[i]<< " "; 
			//outf1<<Q[i][65]<< " "; 
		outf1.close();
		
		/*outf2.open("rwrd"+func(i2)+".txt");
		for(i=0;i<iter_max;i++) 
			outf2<<rwrd[i][54]<< " "; 
		outf2.close();*/
		
	}
	cout<<'\n';
	cout<<'\n';
}

