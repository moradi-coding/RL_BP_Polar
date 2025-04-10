
//find the zeta vectors needed for den_ev4.cpp

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
	int i,a,k,i2,j,l,len_vndeg,len_cndeg,iter,iter_max=120,flg,base,zeta_indx,num_zeta,opt=0,K,num,m=134,***zeta,num_zeta_max=8*pow(10,3),fn,snr_mix=0;
	float R,b,*mu_cn_vn,mu_0,*mu_0_vec,***mu_vn_cn,****mu_vn_cn_mix;
	float **phi_inv,sum,sum1,sum2=0,**phi_mu_vn_cn,**phi_mu_vn_cn_sum_k,tmp,Perr,inc,sigma,sigma_old,
	**rwrd,sum_j1_to_i,epsilon=0.6,beta=0.9,beta_i,*sum_rwd,*V,**Q, nu,
	x,Qx,Qx_vec_sum,*Qx_vec,*x_vec;

	ofstream outf1,outf2;
	int vndeg1[]={1,2,3}, cndeg1[]={4}; //target rate [7,4] Hamming
	float lambda1[]={0.429,0.429,0.142}, rho1[]={1}; 
	
	int vndeg2[]={2}, cndeg2[]={7}; float lambda2[]={1}, rho2[]={1}; fn=8; R=0.714; nu=1; //(2,7) base LDPC code
	//int vndeg2[]={2,0,0}, cndeg2[]={15}; float lambda2[]={1,0,0}, rho2[]={1}; fn=11; R=0.867; nu=1; //(2,15) base LDPC code

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
	for(i = 0; i < len_vndeg; i++) {
    	mu_vn_cn[i] = new float*[m];
   		for(j = 0; j< m; j++) {
    	    mu_vn_cn[i][j] = new float[num_zeta_max];
    	    for(k=0; k<num_zeta_max; k++) { 
    	        mu_vn_cn[i][j][k] = 0;
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
		phi_mu_vn_cn[i]=new float[num_zeta_max];
		phi_mu_vn_cn_sum_k[i]=new float[num_zeta_max];
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
    	    zeta[i][j] = new int[num_zeta_max];
    	    for(k=0; k<num_zeta_max; k++) { 
    	        zeta[i][j][k] = 0;
    	    }
    	}
	}

	base=1;
	num_zeta=0;	
	while(num_zeta<num_zeta_max) { //randomly generate num_zeta_max vectors for each VN degree 
		//for(i=0;i<len_vndeg;i++) 
			//for(j=0;j<m;j++) zeta[i][j][num_zeta]=0; //refresh

		for(i=0;i<len_vndeg;i++) {
			sum2=0;
			while((!base && sum2<vndeg1[i]) || (base && sum2<vndeg2[i])) {
				j=rand()%m;
				cout<<"j: "<<j<<endl;
				if(!zeta[i][j][num_zeta]) {
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
							zeta[i][j][num_zeta]=k;
							sum2+=zeta[i][j][num_zeta];
						}
						else
							flg=0;
					}
				}
			}
		}
		if(base) {
			cout<<'\n'<<"zeta: "; 
			for(i=0;i<len_vndeg;i++) {
				for(j=0;j<m;j++) 
					cout<<zeta[i][j][num_zeta]<<" ";
				cout<<'\n';
			}
		}
		num_zeta++;
		
	}
					
	outf2.open("zeta_mat"+func(i2)+".txt");
	for(k=0;k<num_zeta_max;k++) 
	for(i=0;i<len_vndeg;i++) 
		for(j=0;j<m;j++)  
			outf2<<zeta[i][j][k]<<" "; 
	outf2.close();
		
	cout<<'\n';
    printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n';
}

