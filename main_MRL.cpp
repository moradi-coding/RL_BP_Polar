
//for RL and MQL

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

using namespace std;

double *H,*H2,*x,*x_hat,*y,P,*LR,*pLR,*E_c_v,*E_v_c,*E_c_v_mu,*E_v_c_mu,*res_c_v,*res_c_v_srtd,err,err2,err3,biterr,biterr2,biterr3,varn,ncv,ncv2,ncv3,ncv4,ncv5,ncv6,ncv7,ncv8,
*ncv_vec,*ncv_vec2,*ncv_vec3,*Q,*Q1,*Q2,*Q3,*Q4,*Q5,*Q6,*Q_temp,*Q_temp2,**Q2_inv,qstrt,gap,*Q_cnt,*zeta,*Gtemp,*rep,*diff,**R_sv,*param_vec,*val_vec,
*ncx2_vec,**target,F=1; 
long m,n,n2,*vns,*cns,*vns_D,*cns_D,*excl_cw,BC,S,*indx,*indx2,*states,*syn,*syn_cls,*syn_old,**cls_ind_mat,*spc_indx;
long col_wt,row_wt,num,dcw,dcw2,dcw3,hshft,vshft,mem,W,call,flg,num_cls,L=99,cnt_st=0,**Oc,dc,dv,*cn_indx,bch_sz_max,DeepRL=0,DeepRL2=0,M,cls_sz=1,boxplus=0,cw_cnt=0,samps; 
long **cns_subcode,**vns_subcode,**cns_cluster,**vns_cluster0,**vns_cluster,m_sub,num_sub,meth,**state_sv,num_vns_cls,unif,learn=1,num_gcn,clus_iter=1,*row_wt_vec,*col_wt_vec,
bcjr,*vn_sv,fn,rnd_max;
double trnd=0,L1=0,lth=0,*E_c_v_cnt,R,*E_c_v_cnt2;

double *L_out,*L_in,llr_max=25000, llr_max2=20, pnc=0, alpha, eps;

#define CW 1
#define beta 0.9

ifstream inf,inf2;

string func(long n) {
	stringstream result;
	result << n;
	return result.str();
}

#include "read_fl.cpp"
#include "RL7_meta.cpp" //latest 	

int main() {	
	srand(time(0));	
	clock_t tStart = clock();	
	long modl,snr_idx,glob,flg2,rdc; 
	long num_dat,cnt,cw,load=0,X,Y,i,j,snr_mix;
	//samps=X+Y; //X is samples for snr_mix, Y is samples per SNR
	long gup=100; //no. of global updates for M-RELDEC

	cout<<'\n'<<"fn: "; cin>>fn;
	cout<<'\n'<<"snr_mix: "; cin>>snr_mix; 
	cout<<'\n'<<"bcjr: "; cin>>bcjr; 
	/*if(fn==15) {
		alpha=0.5;
		beta=0.5;
		eps=0.6;
		rnd_max=100;
	}
	else {*/
		alpha=0.1;
		//beta=0.9;
		eps=0.6;
		rnd_max=100;
	//}
	//else {cout<<'\n'<<"matnum?: "; cin>>matnum;}*/
	//cout<<'\n'<<"load Q table: "; cin>>load; //1 loads 0 doesn't
	//cout<<'\n'<<"press 0 for RELDEC, 1 for M-RELDEC, 2 for AM-RELDEC: "; cin>>snr_mix; //0 does RELDEC, 1 does M-RELDEC, 2 does AM-RELDEC
	//snr_mix=1; 
	//bcjr=1; //1 does bcjr based RL, 0 does BP based
	rdc=0;
	//if(snr_mix==2) {cout<<'\n'<<"rdc value: "; cin>>rdc;} //1 does AM-RELDEC-7, 0 does AM-RELDEC-75
		
	
	//if(!snr_mix) {X=0; Y=15000;} //for RL	
	if(snr_mix) { //for RELDEC
		//X=15000; Y=0;
		//X=500000; //180000; 
		X=100000; 
		Y=0;
	}
	else {
		//X=1000; Y=14000;
		X=0; Y=500000;
		//X=30000; Y=150000;
	} //for M-RELDEC
	/*else {
		//X=7500; Y=7500;
		X=30000; Y=150000;
	} *///for AM-RELDEC
	modl=0;
	
	
	//if(!DeepRL) {cout<<'\n'<<"model based?: "; cin>>modl; } //1 yes, 0 no
	BC=1;

	//construct a random generator engine:
	std::random_device rd;
    std::mt19937 e2(rd()); //random seed
    //std::mt19937 e2(3); //fixed seed

	ifstream inf2,inf3,inf4,inf5,inf6,inf7,inf8,inf9,inf10,inf11;  
	
	read_fl(fn);
	E_c_v_cnt2= new double[num_cls];
	
	L_in=new double[n];
	L_out=new double[n];
	vn_sv=new long[row_wt*cls_sz];
	
	if(BC) {L=3;W=3;}

	S=pow(2,num_vns_cls); //total number of states per cluster
	cout<<'\n'<<"S: "<<S<<", bcjr: "<<bcjr;

	ofstream outf,outf2,outf3,outf4,outf5,outf6,outf7,outf8,outf9,outf10; 
	string filename;

	for(j=0;j<m;j++) for(i=0;i<n;i++) inf >> H[j*n+i];
	//cout<<'\n'<<"H: "<<'\n'; for(j=0;j<m;j++){for(i=0;i<n;i++) {cout<<H[j*n+i]<<" "; if(i>0 && (i+1)%(row_wt)==0) cout<<"  ";} cout<<'\n'; if(j>0 && (j+1)%(row_wt)==0) cout<<'\n';}
	
	//removing repeated VNs
	/*for(i=0;i<num_cls;i++) {
		cnt=0;
		for(j=0;j<row_wt*cls_sz;j++) {
			if(!cnt) 
				flg=1;
			else
				flg=0;
			for(k=0;k<cnt;k++)
				if(vns_cluster0[i][j]==vns_cluster[i][k]) {
					flg=0;
					break;
				}
				else
					flg=1;
				
			if(flg && vns_cluster0[i][j]>-1) {
				vns_cluster[i][cnt]=vns_cluster0[i][j]; 
				cnt++;
			}
		}
	}*/

	n2=n;
	hshft=n;
	vshft=m;

	//cout<<'\n'<<"hshft: "<<hshft<<" vshft: "<<vshft<<" n2: "<<n2<<endl; 
	
	//2D arrays
	//x= new double[CW*n2]; //each row has a different CW stream
	x_hat= new double[CW*n2];
	y= new double[CW*n2];  
	LR= new double[CW*n];  //Lvalues in the window
	pLR= new double[CW*n];  
	syn= new long[m];
	syn_old= new long[m];
	syn_cls= new long[cls_sz];
	
	//if(!GLDPC) batch_sz=new long[num_cls];
	//else batch_sz=new long[num_sub];
	excl_cw= new long[CW];
	indx= new long[n];
	indx2= new long[M];
	rep= new double[M];
	diff= new double[M];
	//param_vec= new double[74600];
	//ncx2_vec= new double[74600];
	//val_vec= new double[100];
	

	vns=new long[m*row_wt]; //all VNs of CN j in a row
	
	//3D arrays
	E_v_c= new double[m*row_wt*CW];
	E_c_v= new double[n*col_wt*CW];
	E_v_c_mu= new double[m*row_wt*CW];
	E_c_v_mu= new double[n*col_wt*CW];
	res_c_v= new double[m*row_wt*CW];
	res_c_v_srtd= new double[m*row_wt*CW];
	
	cns=new long[n*col_wt]; //all CNs of a VN in a row

	//if(!DeepRL && !GLDPC) {
		//if(fn==1 || fn==15) {
			Q=new double[S*num_cls]; //S rows, num_cls columns
			Q1=new double[S*num_cls];
			Q2=new double[S*num_cls];
			Q3=new double[S*num_cls];
			Q4=new double[S*num_cls];
			Q5=new double[S*num_cls];
			Q6=new double[S*num_cls];
			Q_cnt=new double[S*num_cls];
			Q_temp=new double[num_cls];
		//}		
		/*else {
			Q=new double[S*num_cls*cls_sz]; //S*num_cls rows, cls_sz columns
			Q_cnt=new double[S*num_cls*cls_sz];
			Q_temp=new double[cls_sz];
		}*/
	//}
	/*else if(GLDPC) {
		Q=new double[S*num_sub*num_sub]; 
		Q_cnt=new double[S*num_sub*num_sub];
		Q_temp=new double[num_sub];
	}*/
	
	//Q2=new double*[M]; for(i=0;i<M;i++) Q2[i]=new double[M]; 
	//tran_prob=new double*[M]; for(i=0;i<M;i++) tran_prob[i]=new double[M]; 
	Q2_inv=new double*[M]; for(i=0;i<M;i++) Q2_inv[i]=new double[2*M]; 
	
	cls_ind_mat=new long*[num_cls]; for(i=0;i<num_cls;i++) cls_ind_mat[i]=new long[cls_sz];
	states= new long[M];
	cn_indx=new long[m];
	spc_indx=new long[num_cls-num_gcn];


	if(load) { // && !DeepRL
		filename="Qtables/Q1_1.txt"; inf10.open(filename.c_str()); 
		//filename="Qtables/Q"+func(fn)+"_1.txt"; inf10.open(filename.c_str()); 
		for(i=0;i<S;i++) 
			for(j=0;j<num_cls;j++) 
				inf10 >> Q[i*num_cls+j];
		//filename="Q_8a.txt"; inf10.open(filename.c_str()); for(i=0;i<S*num_cls;i++) for(j=0;j<cls_sz;j++) inf10 >> Q[i*cls_sz+j];
	} //load Q table from prev. training
		
	//if(fn==15 && M>2) {rep[0]=-6.7; rep[1]=1.11; rep[2]=3.8; rep[3]=6.82;} //for pLLR
	//cout<<'\n'<<"rep: "; for(j=0;j<M;j++) cout<<rep[j]<<" "; cout<<'\n';

	//for(i=0;i<num_cls;i++) for(j=0;j<cls_sz;j++) inf9 >> cls_ind_mat[i][j];
	//cout<<'\n'<<"Q: "<<'\n'; for(i=0;i<S*num_cls;i++) {for(j=0;j<num_cls;j++) cout<<Q[i*num_cls+j]<<" "; cout<<'\n';}

	double EbNo_BC[]={1,2,3,4,5,6}; 
	//EbNo_BC[0]=1; EbNo_BC[1]=2; EbNo_BC[2]=3; EbNo_BC[3]=4; EbNo_BC[4]=4.5; EbNo_BC[5]=5;
	num_dat=sizeof(EbNo_BC)/sizeof(EbNo_BC[0]);
	//Oc=new long*[num_dat]; for(i=0;i<num_dat;i++) Oc[i]=new long[rnd_max];
	cout<<'\n'<<"num_dat: "<<num_dat<<endl;

	if(BC) {
		//initializing
		for(i=0;i<m;i++) for(j=0;j<row_wt;j++) vns[i*row_wt+j]=-1;
		for(i=0;i<n;i++) for(j=0;j<col_wt;j++) cns[i*col_wt+j]=-1;

		for(j=0;j<m;j++) {
			cnt=0; 
			for(i=0;i<n;i++) if(H[j*n+i]) {vns[j*row_wt+cnt]=i; cnt++;} 
		}
		//cout<<'\n'<<"vns "<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<vns[i*row_wt+j]<<" "; cout<<'\n';}
		for(i=0;i<n;i++) {
			cnt=0; 
			for(j=0;j<m;j++) if(H[j*n+i]) {cns[i*col_wt+cnt]=j; cnt++;} 
		}
		cout<<'\n'<<"cns "<<'\n'; for(i=n-1;i<n;i++) {for(j=0;j<col_wt;j++) cout<<cns[i*col_wt+j]<<" "; cout<<'\n';}
	}

	//for(i=0;i<m;i++) 
		//for(j=0;j<M;j++) 
			//rw_vec[i][j]=rw_cnt[i][j]=0; //refresh


//for RL/M-RELDEC
cnt=0;
if(snr_mix) {
	flg2=0;
	for(glob=1;glob>=0;glob--) { 
		//initializing individual Q-functions
		if(!glob) {
			for(i=0;i<S;i++) 
				for(j=0;j<num_cls;j++) {
					Q1[i*num_cls+j]=Q[i*num_cls+j];
					Q2[i*num_cls+j]=Q[i*num_cls+j];
					Q3[i*num_cls+j]=Q[i*num_cls+j];
					Q4[i*num_cls+j]=Q[i*num_cls+j];
					Q5[i*num_cls+j]=Q[i*num_cls+j];
					Q6[i*num_cls+j]=Q[i*num_cls+j];
				}
		}	
		for(snr_idx=0;snr_idx<num_dat;snr_idx++) {
			cout<<'\n'<<"snr_idx: "<<snr_idx<<endl;
			cw_cnt=0;
			if(!glob) samps=Y; 
			else samps=X/num_dat;
			while(cw_cnt<samps) {	
				varn=1/(2*R*pow(10,0.1*EbNo_BC[snr_idx])); //if EbNo not in dB, then varn=1/(2*R*EbNo[i2]);
				//cout<<'\n'<<"snr_idx, varn, sig: "<<snr_idx<<", "<<varn<<", "<<sqrt(varn)<<endl;
				std::normal_distribution<double> dist(0,sqrt(varn));
				for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) y[cw*n2+i]=1+dist(e2); //adding gaussian noise to all-zero CW		
				cw_cnt++;
				//cout<<'\n'<<"cw_cnt: "<<cw_cnt<<endl;
	
				//inf11.open("y.txt"); for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) inf11 >> y[cw*n2+i];
				//cw=0; cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";		
				RL(modl,snr_idx,glob); 		
			}
		}	
		if(snr_mix)
			break;					
	}
}

//for AM-RELDEC
else {
	while(cnt<gup) {
		flg2=0;
		//initializing individual Q-functions
		/*for(i=0;i<S;i++) 
			for(j=0;j<num_cls;j++) {
				Q1[i*num_cls+j]=Q[i*num_cls+j];
				Q2[i*num_cls+j]=Q[i*num_cls+j];
				Q3[i*num_cls+j]=Q[i*num_cls+j];
				Q4[i*num_cls+j]=Q[i*num_cls+j];
				Q5[i*num_cls+j]=Q[i*num_cls+j];
				Q6[i*num_cls+j]=Q[i*num_cls+j];
			}*/
	
		for(glob=0;glob<2;glob++) { 
			//for(snr_idx=0;snr_idx<1/*num_dat*/;snr_idx++) {
			for(snr_idx=num_dat-1;snr_idx<num_dat;snr_idx++) {
				cout<<'\n'<<"snr_idx: "<<snr_idx<<endl;
				cw_cnt=0;
				if(!glob && cnt<gup-1) samps=Y/gup;
				else if(!glob && cnt==gup-1 && rdc) samps/=10; //reduce training in the last step
				else if(glob) samps=X/(num_dat*gup);
				while(cw_cnt<samps) {	
					varn=1/(2*R*pow(10,0.1*EbNo_BC[snr_idx])); //if EbNo not in dB, then varn=1/(2*R*EbNo[i2]);
					//cout<<'\n'<<"snr_idx, varn, sig: "<<snr_idx<<", "<<varn<<", "<<sqrt(varn)<<endl;
					
					//cout<<'\n'<<"std_dev: "<<sqrt(varn)<<endl;
					std::normal_distribution<double> dist(0,sqrt(varn));
					for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) y[cw*n2+i]=1+dist(e2); //adding gaussian noise to all-zero CW		
					cw_cnt++;
					//cout<<'\n'<<"cw_cnt: "<<cw_cnt<<endl;
		
					//inf11.open("y.txt"); for(cw=0;cw<CW;cw++) for(i=0;i<n2;i++) inf11 >> y[cw*n2+i];
					//cw=0; cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<<" ";		
					
					if(glob && !flg2) {						
						//updating the global Q-function
						for(i=0;i<S;i++) 
							for(j=0;j<num_cls;j++) 
								Q[i*num_cls+j]=(Q1[i*num_cls+j]+Q2[i*num_cls+j]+Q3[i*num_cls+j]+Q4[i*num_cls+j]+Q5[i*num_cls+j]+Q6[i*num_cls+j])/6;	
						flg2=1;				
					}
					RL(modl,snr_idx,glob); 		
				}
			}						
		}
		cnt++;	
		cout<<'\n'<<"cnt: "<<cnt<<endl;
	}
}

	if(!Y) {
		if(!bcjr) filename="Q"+func(fn)+"_"+func(X/1000)+"k.txt"; 
		else filename="Q_bcjr"+func(fn)+"_"+func(X/1000)+"k.txt"; 
		outf.open(filename.c_str()); for(i=0;i<S;i++) for(j=0;j<num_cls;j++) outf<<Q[i*num_cls+j]<<" ";
	} //generate only snr_mix Q-table
	else {
		filename="Q"+func(fn)+"1_"+func(Y)+"k.txt"; outf.open(filename.c_str()); for(i=0;i<S;i++) for(j=0;j<num_cls;j++) outf<<Q1[i*num_cls+j]<<" ";
		filename="Q"+func(fn)+"2_"+func(Y)+"k.txt"; outf2.open(filename.c_str()); for(i=0;i<S;i++) for(j=0;j<num_cls;j++) outf2<<Q2[i*num_cls+j]<<" ";
		filename="Q"+func(fn)+"3_"+func(Y)+"k.txt"; outf3.open(filename.c_str()); for(i=0;i<S;i++) for(j=0;j<num_cls;j++) outf3<<Q3[i*num_cls+j]<<" ";
		filename="Q"+func(fn)+"4_"+func(Y)+"k.txt"; outf4.open(filename.c_str()); for(i=0;i<S;i++) for(j=0;j<num_cls;j++) outf4<<Q4[i*num_cls+j]<<" ";
		filename="Q"+func(fn)+"5_"+func(Y)+"k.txt"; outf5.open(filename.c_str()); for(i=0;i<S;i++) for(j=0;j<num_cls;j++) outf5<<Q5[i*num_cls+j]<<" ";
		filename="Q"+func(fn)+"6_"+func(Y)+"k.txt"; outf6.open(filename.c_str()); for(i=0;i<S;i++) for(j=0;j<num_cls;j++) outf6<<Q6[i*num_cls+j]<<" ";
	}
//}

	
	cout<<'\n'<<"trnd: "<<trnd<<endl;
	if(!snr_mix) cout<<'\n'<<"SNR specific"<<endl;
	else cout<<'\n'<<"SNR mix"<<endl;
	//else if(snr_mix==2) cout<<'\n'<<"AM-RELDEC"<<endl;
	
	delete[] x_hat;
	delete[] y;  
	delete[] LR;  
	delete[] pLR;  
	delete[] syn;
	delete[] syn_cls;
	
	delete[] excl_cw;
	delete[] indx;
	delete[] indx2;
	delete[] rep;
	delete[] diff;
	
	delete[] E_v_c;
	delete[] E_c_v;
	delete[] E_v_c_mu;
	delete[] E_c_v_mu;
	delete[] res_c_v;
	delete[] res_c_v_srtd;
	
	delete[] cns; //all CNs of a VN in a row
	delete[] Q; 
	delete[] Q1;
	delete[] Q2;
	delete[] Q3;
	delete[] Q4;
	delete[] Q5;
	delete[] Q_cnt;
	delete[] Q_temp;
	//delete[] EbNo_BC;
	
	cout<<'\n'; 
	printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n'<<"fn "<<fn;

}


