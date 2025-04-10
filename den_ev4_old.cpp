
//density evolution based on Gaussian approximation for sequential decoding seperately on GCN and SPCN

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
	clock_t tStart = clock();	
	int i,a,k,i2,j,l,len_vndeg,len_cndeg,iter,iter_max=120,flg,base,zeta_indx,opt=0,K,num,m=134,***zeta,num_zeta=8*pow(10,3),fn,snr_mix=0;
	float R,b,*mu_cn_vn,mu_0,*mu_0_vec,***mu_vn_cn,****mu_vn_cn_mix;
	float **phi_inv,sum,sum1,sum2=0,**phi_mu_vn_cn,**phi_mu_vn_cn_sum_k,tmp,Perr,inc,sigma,sigma_old,
	**rwrd,sum_j1_to_i,epsilon=0.6,beta=0.9,beta_i,*sum_rwd,*V,**Q, nu,
	x,Qx,Qx_vec_sum,*Qx_vec,*x_vec;
	ifstream inf;
	ofstream outf1,outf2;
	int vndeg1[]={1,2,3}, cndeg1[]={4}; //target rate [7,4] Hamming
	float lambda1[]={0.429,0.429,0.142}, rho1[]={1}; 
	
	int vndeg2[]={2}, cndeg2[]={7}; float lambda2[]={1}, rho2[]={1}; fn=8; R=0.714; nu=1; //(2,7) base LDPC code
	//int vndeg2[]={2,0,0}, cndeg2[]={15}; float lambda2[]={1,0,0}, rho2[]={1}; fn=11; R=0.867; nu=1; //(2,15) base LDPC code
	
	//fn=1; R=0.714; nu=0;
	//fn=2; R=0.501; nu=0.373;
	//fn=5; R=0.288; nu=0.746;
	//fn=7; R=0.151; nu=0.97;
	
	cout<<'\n'<<"snr_mix: "; cin>>snr_mix;
	cout<<"snr_mix: "<<snr_mix<<'\n';
	
	//float EbNo[]={0.98,1,1.1,1.2,1.3,1.4,1.5,2,3,4,5,6,7};
	float EbNo[]={0.98,0.99,1,1.25,1.5,2,3,4,5,6,6.25,6.5};
	//float EbNo[]={6,6.25,6.5};
	//float EbNo[]={0.98,1,1.05,1.1,1.5,2,3,4,5,6,6.5}; 
	K=sizeof(EbNo)/sizeof(EbNo[0]);
	cout<<"K: "<<K<<'\n';
	
	len_vndeg=1; //sizeof(vndeg1)/sizeof(vndeg1[0]);
	len_cndeg=1; //sizeof(cndeg1)/sizeof(cndeg1[0]);
	
	/*a=b=0;
	for(i=0;i<len_cndeg;i++) 
		a+=rho[i]*cndeg[i];
	for(i=0;i<len_vndeg;i++) 
		b+=lambda[i]*vndeg[i];
	R=1-b/a;*/
	
	//mu_vn_cn=new float*[len_vndeg]; 
	//for(i=0;i<len_vndeg;i++) 
		//mu_vn_cn[i]=new float[m];
		
	mu_vn_cn = new float**[len_vndeg];
	mu_vn_cn_mix = new float***[len_vndeg];
	for(i = 0; i < len_vndeg; i++) {
    	mu_vn_cn[i] = new float*[m];
    	mu_vn_cn_mix[i] = new float**[m];
   		for(j = 0; j< m; j++) {
    	    mu_vn_cn[i][j] = new float[num_zeta];
    	    mu_vn_cn_mix[i][j] = new float*[num_zeta];
    	    for(k=0; k<num_zeta; k++) { 
    	        mu_vn_cn[i][j][k] = 0;
    	        mu_vn_cn_mix[i][j][k] = new float[K];
				for(l=0;l<K;l++)
					mu_vn_cn_mix[i][j][k][l]=0; 
    	    }
    	}
	}
	Qx_vec=new float[K];
	x_vec=new float[K];
		
	phi_inv=new float*[len_cndeg];
	for(i=0;i<len_cndeg;i++) 
		phi_inv[i]=new float[m];
		
	phi_mu_vn_cn=new float*[len_vndeg];
	phi_mu_vn_cn_sum_k=new float*[len_vndeg];
	for(i=0;i<len_vndeg;i++) {
		phi_mu_vn_cn[i]=new float[num_zeta];
		phi_mu_vn_cn_sum_k[i]=new float[num_zeta];
	}
		
	Q=new float*[iter_max]; 
	for(i=0;i<iter_max;i++) 
		Q[i]=new float[m];
		
	rwrd=new float*[iter_max]; 
	for(i=0;i<iter_max;i++) 
		rwrd[i]=new float[m];
		
	mu_cn_vn=new float[m];
	sum_rwd=new float[iter_max];
	V=new float[iter_max];
	
	/*zeta=new int*[len_vndeg]; 
	for(i=0;i<len_vndeg;i++) 
		zeta[i]=new int[m];*/
	
	zeta = new int**[len_vndeg];
	for(i=0; i<len_vndeg; i++) {
    	zeta[i] = new int*[m];
   		for(j=0; j<m; j++) {
    	    zeta[i][j] = new int[num_zeta];
    	    for(k=0; k<num_zeta; k++) { 
    	        zeta[i][j][k] = 0;
    	    }
    	}
	}
	
	inf.open("zeta_mat8.txt");
	for(k=0;k<num_zeta;k++) 
	for(i=0;i<len_vndeg;i++) 
		for(j=0;j<m;j++)  
			inf>>zeta[i][j][k]; 
	inf.close();
		
	mu_0_vec= new float[K];
	
	if(snr_mix)
		num=1;
	else
		num=K;
		
	//for(i2=num-3;i2<num;i2++) { //number of different SNRs
	for(i2=0;i2<num;i2++) {
		flg=1;
		if(!snr_mix) {
			sigma=sqrt(1/(2*R*pow(10,0.1*EbNo[i2])));
			cout<<'\n'<<"EbNo: "<<EbNo[i2];
			mu_0=2/(sigma*sigma);
		}
		else {
			for(i=0;i<K;i++) {
				sigma=sqrt(1/(2*R*pow(10,0.1*EbNo[i])));
				mu_0_vec[i]=2/(sigma*sigma);
			}
			cout<<"mu_0_vec: "<<'\n';
			for(i=0;i<K;i++) 
				cout<<mu_0_vec[i]<<" ";
			cout<<'\n';
		}
		
		//refresh/////////////////////////////////////////
		for(i=0;i<len_vndeg;i++)
			for(a=0;a<m;a++) 
				for(zeta_indx=0; zeta_indx<num_zeta; zeta_indx++) 
					mu_vn_cn[i][a][zeta_indx]=0;
					
		for(i=0;i<len_vndeg;i++)
			for(zeta_indx=0; zeta_indx<num_zeta; zeta_indx++) 
				phi_mu_vn_cn[i][zeta_indx]=0;	
				
		for(j=0;j<len_cndeg;j++) 
			for(a=0;a<m;a++) 
				phi_inv[j][a]=0;
		/////////////////////////////////////////////
		
		
		for(iter=0; iter<iter_max; iter++) { //no. of decoder iters.
			for(a=0;a<m;a++) { //a is index of scheduled CN
				zeta_indx=0;				
				base=1; //1 means base code, else overall code
				//base=(rand()%100)<= 100*nu;
				//cout<<"base: "<<base<<endl;
				
				//num_zeta=0;	
				//while(zeta_indx<num_zeta_max) { //randomly generate num_zeta vectors for each VN degree 
				for(zeta_indx=0; zeta_indx<num_zeta; zeta_indx++) {
					/*for(i=0;i<len_vndeg;i++) 
						for(j=0;j<m;j++) zeta[i][j]=0; //refresh

					for(i=0;i<len_vndeg;i++) {
						sum2=0;
						while((!base && sum2<vndeg1[i]) || (base && sum2<vndeg2[i])) {
							j=rand()%m;
							//cout<<"j: "<<j<<'\n';
							if(!zeta[i][j]) {
								flg=0;
								while(!flg) {
									if(!base)
										k=rand()%vndeg1[i];
									else
										k=rand()%vndeg2[i];
									if(vndeg1[i]==1)
										k=1;
									//cout<<"k: "<<k<<endl;
									if((!base && k+sum2<=vndeg1[i]) || (base && k+sum2<=vndeg2[i])) {
										flg=1;
										zeta[i][j]=k;
										sum2+=zeta[i][j];
									}
									else
										flg=0;
								}
							}
						}
					}*/
					/*if(base) {
						cout<<'\n'<<"zeta: "; 
						for(i=0;i<len_vndeg;i++) {
							for(j=0;j<m;j++) 
								cout<<zeta[i][j][zeta_indx]<<" ";
							cout<<'\n';
						}
					}*/
					
					for(i=0;i<len_vndeg;i++) {	
						if(zeta[i][a][zeta_indx]) {
							sum1=sum2=0;
							for(j=0;j<a;j++)
								sum1+=zeta[i][j][zeta_indx]*mu_cn_vn[j]; //CN to VN mean of jth CN
							for(j=a+1;j<m;j++)
								sum2+=zeta[i][j][zeta_indx]*mu_cn_vn[j];
							//mu_vn_cn[i]=mu_0+(vndeg[i]-1)*mu_cn_vn; //mean of VN to CN message	
							if(!snr_mix) {
								mu_vn_cn[i][a][zeta_indx]=mu_0+sum1+(zeta[i][a][zeta_indx]-1)*mu_cn_vn[a]+sum2;
								phi_mu_vn_cn[i][zeta_indx]=exp(-0.4527*pow(mu_vn_cn[i][a][zeta_indx],0.86)+0.0218); 
							}
							else { 
								phi_mu_vn_cn_sum_k[i][zeta_indx]=0;                                
								for(l=0;l<K;l++) {
									mu_vn_cn_mix[i][a][zeta_indx][l]=mu_0_vec[l]+sum1+(zeta[i][a][zeta_indx]-1)*mu_cn_vn[a]+sum2;
									phi_mu_vn_cn_sum_k[i][zeta_indx]+=mu_vn_cn_mix[i][a][zeta_indx][l];
								}
							}
							//cout<<"mu_vn_cn: "<<mu_vn_cn[i][a]<<'\n';
						}
					}
					zeta_indx++;
					//cout<<"sum, a: "<<sum<<", "<<a<<endl;
				}
				sum=0;
				//cout<<"num_zeta: "<<num_zeta<<endl;
				
				for(i=0;i<len_vndeg;i++) 
					for(j=0;j<num_zeta;j++)
						if(!base) {
							if(!snr_mix)
								sum+=lambda1[i]*phi_mu_vn_cn[i][j];
							else
								sum+=lambda1[i]*phi_mu_vn_cn_sum_k[i][j];
						}
						else {
							if(!snr_mix)
								sum+=lambda2[i]*phi_mu_vn_cn[i][j];
							else
								sum+=lambda2[i]*phi_mu_vn_cn_sum_k[i][j];
						}
				if(!snr_mix)
					sum/=num_zeta;
				else
					sum/=(num_zeta*K);
				//cout<<"sum: "<<sum<<endl;
				
				if(sum) { //if the CN is connected to this VN
					for(j=0;j<len_cndeg;j++) {
						//phi_inv[j]=phi_inv(1-pow(1-sum,cndeg[j]-1));
						if(!base) 
							tmp=1-pow(1-sum,cndeg1[j]-1);
						else	
							tmp=1-pow(1-sum,cndeg2[j]-1);
						if(tmp<0.00001)
							tmp=0.00001;
						//cout<<"tmp: "<<tmp<<'\n';
						phi_inv[j][a]=pow((0.0218-log(tmp))/0.4527,1/0.86); //here log means ln	
						//cout<<"mu_cn_vn, a : "<<phi_inv[j][a]<<", "<<a<<'\n';	
					}
					
					mu_cn_vn[a]=0;
					for(j=0;j<len_cndeg;j++)
						if(!base) mu_cn_vn[a]+=rho1[j]*phi_inv[j][a];
						else mu_cn_vn[a]+=rho2[j]*phi_inv[j][a];
					//cout<<"mu_cn_vn, a : "<<mu_cn_vn[a]<<", "<<a<<endl;			
				
					Perr=0;
					for(i=0;i<len_vndeg;i++) {
						//x=sqrt((mu_0+vndeg[i]*mu_cn_vn[a])/2); 
						Qx=0; //Qx1=Qx2=Qx3=Qx4=Qx5=Qx6=Qx7
						for(l=0;l<K;l++) 
							Qx_vec[l]=0;
									
						for(j=0;j<num_zeta;j++) {
							if(!snr_mix) {	
								if(!mu_vn_cn[i][a][j])
									mu_vn_cn[i][a][j]=pow(10,-9);
								x=sqrt(mu_vn_cn[i][a][j]/2); 
								Qx+=((1-exp(-1.98*x))*exp(-0.5*x*x))/(1.135*sqrt(2*pi)*x); //according to wiki
							}
							else {
								for(l=0;l<K;l++) 
									if(!mu_vn_cn_mix[i][a][j][l])
										mu_vn_cn_mix[i][a][j][l]=pow(10,-9);
								
								for(l=0;l<K;l++) {
									x_vec[l]=sqrt(mu_vn_cn_mix[i][a][j][l]/2);
									Qx_vec[l]+=((1-exp(-1.98*x_vec[l]))*exp(-0.5*x_vec[l]*x_vec[l]))/(1.135*sqrt(2*pi)*x_vec[l]); 
								}
							}
						}
						/*if(!snr_mix) 
							cout<<"x, Qx: "<<x<<", "<<Qx<<endl;
						else {
							cout<<"x1, Qx1: "<<x1<<", "<<Qx1<<endl;
							cout<<"x2, Qx2: "<<x2<<", "<<Qx2<<endl;
							cout<<"x3, Qx3: "<<x3<<", "<<Qx3<<endl;
							cout<<"x4, Qx4: "<<x4<<", "<<Qx4<<endl;
							cout<<"x5, Qx5: "<<x5<<", "<<Qx5<<endl;
							cout<<"x6, Qx6: "<<x6<<", "<<Qx6<<endl;
						}*/
						if(!snr_mix) {
							if(!base) Perr+= lambda1[i]*Qx;
							else Perr+= lambda2[i]*Qx;
						}
						else {
							Qx_vec_sum=0;
							for(l=0;l<K;l++) 
								Qx_vec_sum+=Qx_vec[l];
							
							if(!base)  
								Perr+= lambda1[i]*Qx_vec_sum; 
							else
								Perr+= lambda2[i]*Qx_vec_sum; 
						}
					}
					//cout<<"Perr: "<<Perr<<endl;
					//rwrd[iter][a]=1-Perr;
					if(!snr_mix)
						rwrd[iter][a]=1-Perr/num_zeta;
					else
						rwrd[iter][a]=1-Perr/(num_zeta*K);
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
		//cout<<'\n'<<"avg. sum_rwd: "; for(i=0;i<iter_max;i++) cout<<sum_rwd[i]/m<<" ";
		/*if(!snr_mix) 
			outf2.open("rew"+func(fn)+"_"+func(i2)+"_b.txt"); //, std::ios_base::app
		else 
			outf2.open("rew"+func(fn)+"_mix_b.txt");
		for(i=0;i<iter_max;i++) 
			outf2<<sum_rwd[i]/m<<" "; 
		outf2.close();*/
														
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
		if(!snr_mix) outf1.open("V"+func(fn)+"_"+func(i2)+".txt"); //, std::ios_base::app
		else outf1.open("V"+func(fn)+"_mix.txt");
		for(i=0;i<iter_max;i++) 
			outf1<<V[i]<<" "; 
		outf1.close();
		
		/*outf2.open("rwrd"+func(i2)+".txt");
		for(i=0;i<iter_max;i++) 
			outf2<<rwrd[i][54]<< " "; 
		outf2.close();*/
		
	}
	cout<<'\n';
    printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n';
}

