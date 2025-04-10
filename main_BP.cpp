
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
 
double *LR,*E_v_c,*E_c_v,*pLR,*H,*x,*y,err,err2,err3,err4,err5,err6,err7,err8,
biterr,biterr2,biterr3,biterr4,biterr5,biterr6,biterr7,biterr8,*Q,*Q2,*Q3,*Q4,*Q5,*Q6,*Q_cls,*Q_temp,*Q_temp2;

double *ncv_vec,*ncv_vec2,*ncv_vec3,*ncv_vec4,*ncv_vec5,*ncv_vec6,*ncv_vec7,*ncv_vec8,*nvc_vec,*nvc_vec2,*nvc_vec3,*nvc_vec4,
*vn_dist,*cn_dist;
long m,n,n2,*vns,*cns,*excl_cw,BC,S,*indx,*vn_indx,*cn_indx,*indx_cls_CN,*indx_cls,*gcn_indx,*gcn_indx_ordr,*spc_indx,*spc_indx_ordr,
*cls_s_crnt,*cls_s_prev,*pick_cls,*syn,*syn_cls,*mu_cnt,*mu_cnt_L,*perm,pnc;
double varn,ncv,ncv2,ncv3,ncv4,ncv5,ncv6,ncv7,ncv8,nvc,nvc2,nvc3,nvc4,
*E_c_v_cnt,*E_c_v_cnt2,*E_v_c_cnt,*num_iters1,*num_iters2,*num_iters3,*num_iters4,*num_iters5,*num_iters6,*num_iters7,
*num_iters8,L_thld=0,R;
long col_wt,row_wt,dcw,dcw2,dcw3,dcw4,dcw5,dcw6,dcw7,dcw8,hshft,vshft,mem,W,call,meth,num_cls,cls_idx,cls_sz=1,M,M_deep,cnt2=0,j_old,dv,dc,DeepRL=0,boxplus=0,
L=99,num2,iter,num_gcn,**a,*vn_deg;
long **cns_subcode,**vns_subcode,**cns_cluster,**vns_cluster,m_sub,num_sub,concat=0,fn,*x_hat,nan_cnt,
*j_sv,ev,num=0,num_vns_cls,unif,learn=0,clus_iter=1,*row_wt_vec,*col_wt_vec,*vn_sv,not_err=0,Ifl; 

double *L_out, *L_in, llr_max=25000, llr_max2=20; //thresholds for BCJR and BP
long prep_dataset=0, dataset_sz=pow(10,6), cnt_data=0, drl=0; //prepares dataset of pLLR and sch_ord of size max_dec_runs-1 each using MP9.cpp, drl=1 decodes using NN 

#define CW 1 
#define avg 10

ifstream inf,inf2;

string func(long n) {
	stringstream result;
	result << n;
	return result.str();
}

std::vector<double> inp;

#include "read_fl.cpp"
#include "BP.cpp"


int main() {	
	srand(time(0));	
	clock_t tStart = clock();	
	long i,i2,j,k,flg,num_dat,cnt,cnt2,cw,l1,l2,l3,l,mixture=1,m1=1,m3=1,m4=1, m2=0,m5=0,m6=0,m7=0,m8=0, //simulate BP and BCJR separately
	run,run2,run3,run4,run5,run6,run7,run8,dcw_max=pow(10,6);
	long max_dec_runs=10*pow(10,7);
	double msg_cnt, msg_cnt2, cnt_tot,*y_vec_in;
	
	double *blk_err,*blk_err2,*blk_err3,*blk_err4,*blk_err5,*blk_err6,*blk_err7,*blk_err8,
	       *bit_err,*bit_err2,*bit_err3,*bit_err4,*bit_err5,*bit_err6,*bit_err7,*bit_err8;  
	if(DeepRL) cout<<'\n'<<"Copy json files from sub-folder DeepRL/codename/ to the DeepRL2 folder";
	cout<<'\n'<<"file no: "; cin>>fn; ///////////////////////////////// 
	//fn=11; 
	pnc=0; //for [7,4] Hamming code BCJR test- use method 5 for this
	cout<<"llr_max: "<<llr_max<<", llr_max2: "<<llr_max2<<endl;
	ofstream outf1,outf2,outf3,outf4,outf5,outf6,outf7,outf8,outf9;
	if(fn==13 || fn==14 || fn==15 || fn==16 || fn==17 || fn==18 || fn==19 || fn==20 || fn==21 || fn==22 || fn==23 || fn==24 || fn==25 || fn==26
	          || fn==27 || fn==28 || fn==29)
		Ifl=200;
	else
		Ifl=50;
	
	//construct a random generator engine:
	std::random_device rd;
   	std::mt19937 e2(rd()); //random seed
   	//std::mt19937 e2(3); //fixed seed

	ifstream inf3,inf4,inf5,inf6,inf7,inf8,inf9,inf10,inf11,inf12,inf13,inf14,inf15,inf16,inf17,inf18,inf20,inf21;
	read_fl(fn);
	y_vec_in=new double[6*n]; 
	
	ofstream outf,out_file2; string filename;
	S=pow(2,num_vns_cls);

	Q=new double[S*num_cls]; 
	Q2=new double[S*num_cls];
	Q3=new double[S*num_cls]; 
	Q4=new double[S*num_cls]; 
	Q5=new double[S*num_cls]; 
	
	Q_temp=new double[num_cls];
	Q_cls=new double[num_cls];
	
	L_in=new double[n];
	L_out=new double[n];
	vn_sv=new long[row_wt*cls_sz];
	
	cout<<'\n'<<"fn: "<<fn<<" num_gcn: "<<num_gcn<<" pnc: "<<pnc<<" Ifl: "<<Ifl<<" mixture: "<<mixture<<" max_dec_runs: "<<max_dec_runs<<endl;

	//here n=WJp*p
	n2=n;
	hshft=n;
	vshft=m;
	
	num2=hshft;

	//cout<<'\n'<<"hshft: "<<hshft<<" vshft: "<<vshft<<" n2: "<<n2<<endl; 
	
	//2D arrays
	x_hat= new long[CW*n2];
	y= new double[CW*n2];  
	LR= new double[CW*n]; 
	pLR= new double[CW*n]; 
	syn=new long[m];
	syn_cls=new long[cls_sz]; 
	excl_cw=new long[CW];
	j_sv= new long[m];
	
	vn_indx=new long[m];
	cn_indx=new long[m];
	spc_indx=new long[num_cls-num_gcn];
	spc_indx_ordr=new long[num_cls-num_gcn];
	gcn_indx=new long[num_gcn];
	gcn_indx_ordr=new long[num_gcn];
	indx=new long[cls_sz];
	indx_cls=new long[num_cls];
	indx_cls_CN=new long[num_cls];
	cls_s_crnt=new long[num_cls];
	cls_s_prev=new long[num_cls];
	pick_cls=new long[num_cls];
	
	vns=new long[m*row_wt]; //all VNs of CN j in a row
	//3D arrays
	E_v_c= new double[m*row_wt*CW]; //each row is a CN, and the entries of a row are its VNs
	E_v_c_cnt= new double[m*row_wt*CW];
	
	E_c_v= new double[n*col_wt*CW]; //each row is a VN, and the entries of a row are its CNs
	E_c_v_cnt= new double[n*col_wt*CW];
	E_c_v_cnt2= new double[num_cls];
	
	cns=new long[n*col_wt]; //all CNs of a VN in a row
	
	double EbNo_BC[8];
	if(fn==-1 || fn==-4 || fn==-5) { //BCJR testing of [7,4] and [15,11] Hamming codes
		EbNo_BC[0]=1; EbNo_BC[1]=3; EbNo_BC[2]=5; EbNo_BC[3]=6; EbNo_BC[4]=7; EbNo_BC[5]=8; EbNo_BC[6]=9; num_dat=7;
	}
	else if(fn==-2) { 
		EbNo_BC[0]=1; EbNo_BC[1]=3; EbNo_BC[2]=5; EbNo_BC[3]=7; EbNo_BC[4]=9; EbNo_BC[5]=11; EbNo_BC[6]=12; num_dat=7;
		//EbNo_BC[0]=1; EbNo_BC[1]=2; EbNo_BC[2]=2.5; EbNo_BC[3]=3; EbNo_BC[4]=3.5; EbNo_BC[5]=4; EbNo_BC[6]=4.5; EbNo_BC[7]=5; num_dat=8;
	}
	else if(fn==-3) { 
		EbNo_BC[0]=0; EbNo_BC[1]=0.5; EbNo_BC[2]=1; EbNo_BC[3]=1.5; EbNo_BC[4]=1.75; num_dat=5;
		//EbNo_BC[0]=1; EbNo_BC[1]=2; EbNo_BC[2]=2.5; EbNo_BC[3]=3; EbNo_BC[4]=3.5; EbNo_BC[5]=4; EbNo_BC[6]=4.5; EbNo_BC[7]=5; num_dat=8;
	}
	else if(fn==8) {EbNo_BC[0]=1; EbNo_BC[1]=2; EbNo_BC[2]=2.5; EbNo_BC[3]=3; EbNo_BC[4]=3.5; EbNo_BC[5]=4; EbNo_BC[6]=4.5; EbNo_BC[7]=4.75; num_dat=8;}
	else if(fn==10) {EbNo_BC[0]=1; EbNo_BC[1]=2; EbNo_BC[2]=2.5; EbNo_BC[3]=2.75; EbNo_BC[4]=3; EbNo_BC[5]=3.1; num_dat=6;}
	else if(fn==11) {EbNo_BC[0]=1; EbNo_BC[1]=1.25; EbNo_BC[2]=1.5; EbNo_BC[3]=1.75; EbNo_BC[4]=2; EbNo_BC[5]=2.1; EbNo_BC[6]=2.2; num_dat=7;}
	else if(fn==12) {EbNo_BC[0]=1; EbNo_BC[1]=3; EbNo_BC[2]=5; EbNo_BC[3]=6; EbNo_BC[4]=7; num_dat=5;} //EbNo_BC[5]=8; 
	else if(fn==13 || fn==14 || fn==15 || fn==16 || fn==17 || fn==18 || fn==19 || fn==20 || fn==21 || fn==22 || fn==23 || fn==24 || fn==25 || fn==26
	               || fn==27 || fn==28 || fn==29) {
		EbNo_BC[0]=1; EbNo_BC[1]=2; EbNo_BC[2]=3; EbNo_BC[3]=4; EbNo_BC[4]=5; EbNo_BC[5]=6; EbNo_BC[6]=7; num_dat=7;
	}
	cout<<'\n'<<"EbNo_BC: "; for(i=0;i<num_dat;i++) cout<<EbNo_BC[i]<< " "; cout<<'\n';
	
	num_iters1=new double[num_dat];
	num_iters2=new double[num_dat];
	num_iters3=new double[num_dat];
	num_iters4=new double[num_dat];
	num_iters5=new double[num_dat];
	num_iters6=new double[num_dat];
	num_iters7=new double[num_dat];
	num_iters8=new double[num_dat];

	blk_err= new double[num_dat]; 
	blk_err2= new double[num_dat];
	blk_err3= new double[num_dat]; 
	blk_err4= new double[num_dat]; 
	blk_err5= new double[num_dat];
	blk_err6= new double[num_dat];
	blk_err7= new double[num_dat];
	blk_err8= new double[num_dat];

	bit_err= new double[num_dat]; 
	bit_err2= new double[num_dat]; 	
	bit_err3= new double[num_dat]; 
	bit_err4= new double[num_dat]; 
	bit_err5= new double[num_dat];
	bit_err6= new double[num_dat];
	bit_err7= new double[num_dat];
	bit_err8= new double[num_dat];

	ncv_vec=new double[num_dat];
	ncv_vec2=new double[num_dat];
	ncv_vec3=new double[num_dat];
	ncv_vec4=new double[num_dat];
	ncv_vec5=new double[num_dat];
	ncv_vec6=new double[num_dat];
	ncv_vec7=new double[num_dat];
	ncv_vec8=new double[num_dat];
	
	nvc_vec=new double[num_dat];
	nvc_vec2=new double[num_dat];
	nvc_vec3=new double[num_dat];
	nvc_vec4=new double[num_dat];
	
	vn_deg=new long[n];
	vn_dist=new double[n];
	cn_dist=new double[m];
	
	for(i=0;i<m;i++) for(j=0;j<row_wt;j++) vns[i*row_wt+j]=-1;
	for(i=0;i<n;i++) for(j=0;j<col_wt;j++) cns[i*col_wt+j]=-1;

	cnt_tot=0;
	for(j=0;j<m;j++) {
		cnt=0; 
		for(i=0;i<n;i++) 
			if(H[j*n+i]) {
				vns[j*row_wt+cnt]=i; 
				cnt++;
			} 
		cnt_tot+=cnt;
		cn_dist[cnt-1]++;
	}
	//cout<<'\n'<<"vns "<<'\n'; for(i=0;i<m;i++) {for(j=0;j<row_wt;j++) cout<<vns[i*row_wt+j]<< " "; cout<<'\n';}
	cout<<'\n'<<"avg. CN deg: "<<cnt_tot/m; 
	
	cnt_tot=0;
	for(i=0;i<n;i++) {
		cnt=0; 
		for(j=0;j<m;j++) 
			if(H[j*n+i]) {cns[i*col_wt+cnt]=j; cnt++;} 
		cnt_tot+=cnt;
		//cout<<'\n'<<"VN deg: "<<cnt; 
		vn_deg[i]=cnt;
		vn_dist[cnt-1]++;
	}
	//cout<<'\n'<<"cns "<<'\n'; for(i=n-1;i<n;i++) {for(j=0;j<col_wt;j++) cout<<cns[i*col_wt+j]<< " "; cout<<'\n';}
	cout<<'\n'<<"avg. VN deg: "<<cnt_tot/n; 
	
	cout<<'\n'<<"vn_dist: "<<'\n';
	for(i=0;i<n;i++) 
		if(vn_dist[i])
			cout<<"prob: "<<vn_dist[i]/n<<", deg: "<<i+1<<'\n';
			
	cout<<'\n'<<"cn_dist: "<<'\n';
	for(i=0;i<m;i++) 
		if(cn_dist[i])
			cout<<"prob: "<<cn_dist[i]/m<<", deg: "<<i+1<<'\n';
	
	cnt=cnt2=0;
	for(i=0;i<num_cls;i++) {
		if(cls_sz==1) {
			spc_indx[cnt]=i;
			cnt++;
		}
		else {
			for(j=0;j<cls_sz;j++) 
				if(cns_cluster[i][j]<0) {
					flg=0;
					break;
				}
				else 
					flg=1;
			if(flg) {
				gcn_indx[cnt2]=i;
				cnt2++;
			}
			else {
				spc_indx[cnt]=i;
				cnt++;
			}
		}
	}
	
	//cout<<'\n'<<"num_cls-num_gcn: "<<num_cls-num_gcn<<'\n'; 
	//cout<<'\n'<<"spc_indx: "; for(j=0;j<num_cls-num_gcn;j++) cout<<spc_indx[j]<<" "; 
	//cout<<'\n'<<"gcn_indx: "; for(j=0;j<num_gcn;j++) cout<<gcn_indx[j]<<" "; 
	
	//finding the 'a' vector for 
	/*a=new long*[n]; 
	for(i=0;i<n;i++) 
		a[i]=new long[num_cls]; 
	for(i=0;i<n;i++) {
		for(j=0;j<col_wt;j++) {
			for(k=0;k<num_cls;k++) { 
				for(j2=0;j2<cls_sz;j2++) 
					if(cns_cluster[k][j2]>-1 && cns_cluster[k][j2]==cns[i*col_wt+j]) {
						clsidx=k;
						flg=1;
						break;
					}
					else
						flg=0;
			
			//for(k=0;k<num_cls-num_gcn;k++) 
				//if(clsidx==spc_indx[k])
				if(flg)
					break;
			}
			if(flg) 
				a[i][clsidx]++;
		}
	}
	cout<<'\n'<<"a: "<<'\n'; 
	for(i=0;i<n;i++) {
		for(j=0;j<num_cls;j++) 
			cout<<a[i][j]<<" "; 
		cout<<'\n';
	}
	int a_notdv=0,a_dv=0;
	for(i=0;i<n;i++) {
		if(vn_deg[i]==3) {
			for(j=0;j<num_cls;j++) {
				if(a[i][j]==3) {
					a_dv++; 
					flg=0;
					break;
				}
				else 
					flg=1;
			}
			if(flg)
				a_notdv++;
		}
	}
	cout<<'\n'<<"a_notdv,a_dv: "<<a_notdv<<", "<<a_dv;*/
	//cout<<'\n'<<"vn_deg: "; for(i=0;i<n;i++) cout<<vn_deg[i]<<" "; 
	
	if(mixture) { //not SNR specific
		if(fn==1) 
			inf15.open("Qtables/Q1_1.txt");
		else if(fn==2) 
			inf15.open("Qtables/Q2_1.txt");
		else if(fn==3) 
			inf15.open("Qtables/Q3_1.txt");
		else if(fn==4) 
			inf15.open("Qtables/Q4_1.txt");
		else if(fn==5) 
			inf15.open("Qtables/Q5_1.txt");
		else if(fn==6) {
			inf15.open("Qtables/Q6_300k.txt");
			inf9.open("Qtables/Q_bcjr6_450k.txt");
		}
		else if(fn==7) 
			inf15.open("Qtables/Q7_1.txt"); //trained with clus_iter=1
		else if(fn==8) {
			//inf15.open("Qtables/Q8_300k.txt"); //for BP
			inf15.open("Qtables/Q8_500k.txt"); cout<<'\n'<<"Qtable: 1-6dB mix";
			//inf9.open("Qtables/Q_bcjr8_450k_1dB.txt"); cout<<'\n'<<"Qtable: 1dB";
			inf9.open("Qtables/Q_bcjr8_500k_1_6dB.txt"); cout<<'\n'<<"Qtable: 1-6dB mix";
			//inf9.open("Qtables/Q_bcjr8_500k_6dB.txt"); cout<<'\n'<<"Qtable: 6dB";
		}
		else if(fn==9) {
			inf15.open("Qtables/Q9_300k.txt");
			inf9.open("Qtables/Q_bcjr9_300k.txt");
		}
		else if(fn==10) {
			inf15.open("Qtables/Q10_500k.txt");
			inf9.open("Qtables/Q_bcjr10_500k.txt");
		}
		else if(fn==11) {
			inf15.open("Qtables/Q11_500k.txt");
			inf9.open("Qtables/Q_bcjr11_450k.txt");
		}
		else if(fn==12) {
			inf15.open("Qtables/Q12_5k.txt");
			//inf9.open("Qtables/Q_bcjr11_450k.txt");
		}
		else if(fn==13) 
			inf15.open("Qtables/Q13_100k.txt");
		else if(fn==14) 
			inf15.open("Qtables/Q14_100k.txt");
		else if(fn==15) 
			inf15.open("Qtables/Q15_100k.txt");
		else if(fn==16) 
			inf15.open("Qtables/Q16_100k.txt");
		else if(fn==17) 
			inf15.open("Qtables/Q17_100k.txt");
		else if(fn==18) 
			inf15.open("Qtables/Q18_100k.txt");
		else if(fn==19) 
			inf15.open("Qtables/Q19_100k.txt");
		else if(fn==20) 
			inf15.open("Qtables/Q20_100k.txt");
		else if(fn==21) 
			inf15.open("Qtables/Q21_100k.txt");
		else if(fn==22) 
			inf15.open("Qtables/Q22_100k.txt");
		else if(fn==23) 
			inf15.open("Qtables/Q23_100k.txt");
		else if(fn==24) 
			inf15.open("Qtables/Q24_100k.txt");
		else if(fn==25) 
			inf15.open("Qtables/Q25_100k.txt");
		else if(fn==26) 
			inf15.open("Qtables/Q26_100k.txt");
		else if(fn==27) 
			inf15.open("Qtables/Q27_100k.txt");
		else if(fn==28) 
			inf15.open("Qtables/Q28_100k.txt");
		else if(fn==29) 
			inf15.open("Qtables/Q29_100k.txt");
			
		for(i=0;i<S;i++) 
			for(j=0;j<num_cls;j++) {
				inf15 >> Q[i*num_cls+j]; //for RELDEC mixture- meth4
				inf9 >> Q2[i*num_cls+j]; //for BCJR based learning
			}
	}
	
	else {
		if(fn==5) {
			inf3.open("Qtab_transfer/Q5_1.txt");
			inf4.open("Qtab_transfer/Q5_2.txt");
			inf5.open("Qtab_transfer/Q5_3.txt");
			inf6.open("Qtab_transfer/Q5_4.txt");
			inf7.open("Qtab_transfer/Q5_5.txt");
			inf8.open("Qtab_transfer/Q5_6.txt");
		}
		else if(fn==6) {
			inf3.open("Qtab_transfer/Q6_1.txt");
			inf4.open("Qtab_transfer/Q6_2.txt");
			inf5.open("Qtab_transfer/Q6_3.txt");
			inf6.open("Qtab_transfer/Q6_4.txt");
			inf7.open("Qtab_transfer/Q6_5.txt");
			inf8.open("Qtab_transfer/Q6_6.txt");
		}
		else if(fn==8) {
			inf3.open("Qtab_transfer/Q8_1.txt");
			inf4.open("Qtab_transfer/Q8_2.txt");
			inf5.open("Qtab_transfer/Q8_3.txt");
			inf6.open("Qtab_transfer/Q8_4.txt");
			inf7.open("Qtab_transfer/Q8_5.txt");
			inf8.open("Qtab_transfer/Q8_6.txt");
		}
	
	}	
	
	double avg_Q[num_cls];
	for(j=0;j<num_cls;j++)
		avg_Q[j]=0;
	
	if(mixture) {	
		for(j=0;j<num_cls;j++)
			for(i=0;i<S;i++) 
				avg_Q[j]+=Q2[i*num_cls+j];
				
		for(j=0;j<num_cls;j++)
			avg_Q[j]/=S;
			
		//sorting average	
		/*for(i=0;i<num_cls;i++) indx[i]=i; //refresh
		for(i=0;i<num_cls;i++) {
			for(i2=i+1;i2<num_cls;i2++) {
				if(avg_Q[i]>avg_Q[i2]) {   
					tmp=avg_Q[i];
					avg_Q[i]=avg_Q[i2];
					avg_Q[i2]=tmp;
					
					//tmp2=indx[i];
					//indx[i]=indx[i2];
					//indx[i2]=tmp2;
				}
							
			}
		}*/
			
		//cout<<'\n'<<"avg_Q spc_indx: "; for(j=0;j<num_cls-num_gcn;j++) cout<<avg_Q[spc_indx[j]]<<" "; 
		//cout<<'\n'<<'\n'<<"avg_Q gcn_indx: "; for(j=0;j<num_gcn;j++) cout<<avg_Q[gcn_indx[j]]<<" "; 
		//cout<<'\n'<<"indx: "; for(j=0;j<num_cls;j++) cout<<indx[j]<<" "; 
	}
	
	/*inf11.open("y_mat.txt"); 
	for(j=0;j<6;j++) 
		for(i=0;i<n;i++) 
			inf11 >> y_vec_in[j*n+i];
			
	inf11.open("y_special.txt"); 
	for(i=0;i<n;i++) 
		inf11 >> y_vec_in[0*n+i];*/
	
	/*cout<<"y_vec_in: "<<'\n';
	for(j=0;j<13455;j++) {	
		for(i=0;i<n;i++) 
			cout<<y_vec_in[j*n+i]<<" "; 
		cout<<'\n';
	}*/
	
	for(i2=0;i2<num_dat;i2++) {
	//for(i2=0;i2<5;i2++) {
	//for(i2=num_dat-2;i2<num_dat;i2++) {
		run=run2=run3=run4=run5=run6=run7=run8=1;		
		if(fn==11 && pnc) {
			R=0.5;
			varn=1/(2*R*pow(10,0.1*(EbNo_BC[i2]+log10(R)))); //doing Es/No 
		}
		else varn=1/(2*R*pow(10,0.1*EbNo_BC[i2]));
		cout<<'\n'<<"i2, varn, sig: "<<i2<<", "<<varn<<", "<<sqrt(varn)<<endl;
		//min. no. of error events
		if(fn==-1 || fn==-4) {
			//if(i2<3) ev=200; else if(i2==3) ev=150; else if(i2>=4) ev=100;
			//if(i2<3) ev=100; else if(i2==3) ev=50; else if(i2>=4) ev=25;
			ev=5;
		}
		else if(fn==-3) {
			if(i2<3) ev=25; else if(i2<5) ev=25; else if(i2<5) ev=10; else ev=3;
		}
		else if(fn==8) {
			//if(i2<3) ev=10; else if(i2<5) ev=10; else if(i2<6) ev=5; else ev=1;
			if(m1 || m3 || m4) {if(i2<3) ev=100; else if(i2<5) ev=75; else ev=50;}
			//if(m2 || m5 || m6) {if(i2<3) ev=10; else if(i2<5) ev=5; else if(i2<6) ev=1; else ev=1;}
			if(m2 || m5 || m6) {if(i2<3) ev=50; else if(i2<5) ev=25; else if(i2<6) ev=10; else ev=3;}
		}
		else if(fn==10){
			if(m1 || m3 || m4) {if(i2<3) ev=50; else if(i2<5) ev=25; else ev=10;}
			if(m2 || m5) {if(i2<2) ev=50; else if(i2<3) ev=25; else ev=10;}
			else if(m6) {if(i2<2) ev=5; else if(i2<3) ev=3; else ev=1;}
			//if(i2<3) ev=100; else if(i2<5) ev=50; else ev=20;
		}
		else if(fn==11) {
			if(m1 || m3 || m4) {if(i2<3) ev=25; else if(i2<5) ev=25; else ev=25;}
			//if(m2 || m5 || m6) {if(i2<3) ev=25; else if(i2<6) ev=3; else ev=1;}
			if(m2 || m5 || m6) {if(i2<3) ev=10; else if(i2<6) ev=3; else ev=1;}
		}
		else if(fn==12) {
			if(i2<3) ev=25; else if(i2<5) ev=10; else ev=3;
		}
		else if(fn==13 || fn==14 || fn==16 || fn==17 || fn==18 || fn==19 || fn==20 || fn==21 || fn==22 || fn==23 || fn==24 || fn==25
		               || fn==28 || fn==29) {
			if(i2<3) ev=150; else if(i2<5) ev=100; else if(i2<5) ev=25; else ev=10;
		}
		else if(fn==26) {
			if(i2<3) ev=75; else if(i2<5) ev=50; else if(i2<5) ev=12; else ev=5;
		}
		else if(fn==15) {
			if(i2<3) ev=10; else if(i2<5) ev=10; else if(i2<5) ev=5; else ev=3;
		}
		else if(fn==27) {
			if(i2<3) ev=10; else if(i2<4) ev=5; else if(i2<5) ev=3; else ev=1;
		}
		//else varn=1/(2*R*pow(10,0.1*EbNo_kk[i2])); //if EbNo not in dB, then varn=1/(2*R*EbNo_kk[i2]);
							
		if(!mixture) {					
			for(i=0;i<S;i++) 
				for(j=0;j<num_cls;j++) 
					if(!i2) inf3 >> Q[i*num_cls+j]; 
					else if(i2==1) inf4 >> Q[i*num_cls+j]; 
					else if(i2==2) inf5 >> Q[i*num_cls+j]; 
					else if(i2==3) inf6 >> Q[i*num_cls+j]; 
					else if(i2==4) inf7 >> Q[i*num_cls+j]; 
					else if(i2==5) inf8 >> Q[i*num_cls+j]; 
					
			for(j=0;j<num_cls;j++)
				for(i=0;i<S;i++) 
					avg_Q[j]+=Q[i*num_cls+j];
				
			for(j=0;j<num_cls;j++)
				avg_Q[j]/=S;
			
			//cout<<'\n'<<"avg_Q spc_indx: "; for(j=0;j<num_cls-num_gcn;j++) cout<<avg_Q[spc_indx[j]]<<" "; 
			//cout<<'\n'<<'\n'<<"avg_Q gcn_indx: "; for(j=0;j<num_gcn;j++) cout<<avg_Q[gcn_indx[j]]<<" ";			
		}
		
		std::normal_distribution<double> dist(0,sqrt(varn));
		ncv=ncv2=ncv3=ncv4=ncv5=ncv6=ncv7=ncv8=0;
		nvc=nvc2=nvc3=nvc4=0;
		dcw=dcw2=dcw3=dcw4=dcw5=dcw6=dcw7=dcw8=nan_cnt=0;
		err=err2=err3=err4=err5=err6=err7=err8=biterr=biterr2=biterr3=biterr4=biterr5=biterr6=biterr7=biterr8=l1=l2=l3=0;
		//for(cnt=0;cnt<6;cnt++) {
		while((m1 && err<ev && run) || (m2 && err2<ev && run2) || (m3 && err3<ev && run3) || (m4 && err4<ev && run4) || 
		      (m5 && err5<ev && run5) || (m6 && err6<ev && run6) || (m7 && err7<ev && run7)) { 
			if(fn==0) { //for test
				cw=0; y[cw*n+0]=1.1; y[cw*n+1]=1.1; y[cw*n+2]=0.5; y[cw*n+3]=1.1; y[cw*n+4]=0.5; //for SPA
				//cw=0; y[cw*n+0]=-1.1; y[cw*n+1]=3.1; y[cw*n+2]=-2.5; y[cw*n+3]=1.1; y[cw*n+4]=0.5; //for min-sum 128 172 183 211 338 362 427 473 513 526 564 854
			}
			else {
				for(cw=0;cw<CW;cw++) 
					for(i=0;i<n2;i++) y[cw*n2+i]=1+dist(e2); //adding gaussian noise to all-zero CW
			}
			//inf11.open("y.txt"); for(cw=0;cw<CW;cw++) for(i=0;i<n;i++) inf11 >> y[cw*n+i];
			//for(i=0;i<n;i++) 
				//y[i]=y_vec_in[cnt*n+i];
			//cw=0; cout<<'\n'<<"y: "; for(i=0;i<n2;i++) cout<<y[cw*n2+i]<< " "; cout<<'\n';
			//cw=0; ofstream outf; outf.open("y.txt"); for(i=0;i<n2;i++) outf<<y[cw*n2+i]<< " "; outf<<std::endl; outf.close();

			for(meth=1;meth<=8;meth++) {	
				if((meth==1 && m1 && run && err<ev) || (meth==2 && m2 && run2 && err2<ev) || (meth==3 && m3 && run3 && err3<ev) || (meth==4 && m4 && run4 && err4<ev) || 
				   (meth==5 && m5 && run5 && err5<ev) || (meth==6 && m6 && run6 && err6<ev) || (meth==7 && m7 && run7 && err7<ev)) {				
					cw=0; 
					for(j=0;j<n;j++) 
						for(k=0;k<col_wt;k++) 
							E_c_v_cnt[CW*(col_wt*j+k)+cw]=0; //refresh	
					for(j=0;j<num_cls;j++) E_c_v_cnt2[j]=0;			

					for(j=0;j<m;j++) 
						for(k=0;k<row_wt;k++) 						
							E_v_c_cnt[CW*(row_wt*j+k)+cw]=0;										
														
					l=BP(i2,num_dat/*,model0*/); 

					cw=msg_cnt=msg_cnt2=0; 
					//for(j=0;j<n;j++) for(k=0;k<col_wt;k++) msg_cnt+=E_c_v_cnt[CW*(col_wt*j+k)+cw];
					for(i=0;i<num_cls;i++) msg_cnt+=E_c_v_cnt2[i];
					//cout<<"E_c_v_cnt2: "; for(i=0;i<num_cls;i++) cout<<E_c_v_cnt2[i]<< " "; cout<<endl;
					
					for(j=0;j<m;j++) for(k=0;k<row_wt;k++) msg_cnt2+=E_v_c_cnt[CW*(row_wt*j+k)+cw];
					
					if(meth==1) ncv+=msg_cnt;
					else if(meth==2) ncv2+=msg_cnt;
					else if(meth==3) ncv3+=msg_cnt;
					else if(meth==4) ncv4+=msg_cnt;
					else if(meth==5) ncv5+=msg_cnt;
					else if(meth==6) ncv6+=msg_cnt;
					else if(meth==7) ncv7+=msg_cnt;
					else if(meth==8) ncv8+=msg_cnt;
					
					if(meth==1) nvc+=msg_cnt2;
					else if(meth==2) nvc2+=msg_cnt2;
					else if(meth==3) nvc3+=msg_cnt2;
					else if(meth==4) nvc4+=msg_cnt2;
					
					if(meth==1) num_iters1[i2]+=l;
					else if(meth==2) num_iters2[i2]+=l;
					else if(meth==3) num_iters3[i2]+=l;
					else if(meth==4) num_iters4[i2]+=l;
					else if(meth==5) num_iters5[i2]+=l;
					else if(meth==6) num_iters6[i2]+=l;
					else if(meth==7) num_iters7[i2]+=l;
					else if(meth==8) num_iters8[i2]+=l;
					//cout<<'\n'<<"err: "<<err<<" err2: "<<err2<<" err3: "<<err3<<" err4: "<<err4<<" meth: "<<meth<<endl;
					if(meth==1 && !(dcw%dcw_max)) 
						cout<<'\n'<<"err: "<<err<<endl;
					else if(meth==2 && !(dcw2%dcw_max)) 
						cout<<'\n'<<"err2: "<<err2<<endl;
					else if(meth==3 && !(dcw3%dcw_max)) 
						cout<<'\n'<<"err3: "<<err3<<endl;
					else if(meth==4 && !(dcw4%dcw_max)) 
						cout<<'\n'<<"err4: "<<err4<<endl;
					else if(meth==5 && !(dcw5%dcw_max)) 
						cout<<'\n'<<"err5: "<<err5<<endl;
					else if(meth==6 && !(dcw6%dcw_max)) 
						cout<<'\n'<<"err6: "<<err6<<endl;
					else if(meth==7 && !(dcw7%dcw_max)) 
						cout<<'\n'<<"err7: "<<err7<<endl;
					else if(meth==8 && !(dcw8%dcw_max)) 
						cout<<'\n'<<"err8: "<<err8<<endl;
				}
			}
			
			if(m1 && dcw>max_dec_runs)
				run=0;
			 else
			  	run=1; 
			 if(m2 && dcw2>max_dec_runs)
				run2=0;
			 else
			  	run2=1; 
			 if(m3 && dcw3>max_dec_runs)
				run3=0;
			 else
			  	run3=1; 
			 if(m4 && dcw4>max_dec_runs)
				run4=0;
			 else
			  	run4=1; 
			 if(m5 && dcw5>max_dec_runs)
				run5=0;
			 else
			  	run5=1; 
			 if(m6 && dcw6>max_dec_runs)
				run6=0;
			 else
			  	run6=1; 
			 if(m7 && dcw7>max_dec_runs)
				run7=0;
			 else
			  	run7=1; 
			 if(m8 && dcw8>max_dec_runs)
				run8=0;
			 else
			  	run8=1; 	  	
		} //while loop
		
		blk_err[i2]=err/dcw;
		bit_err[i2]=biterr/(num2*dcw);
		ncv_vec[i2]=ncv/(1*dcw); //no. of CN to VN msg. updates
		nvc_vec[i2]=nvc/(1*dcw); 
		//if(i2==num_dat-1) {for(i=0;i<num2;i++) loc_err1[i]/=dcw; /**/} 
		//cout<<'\n'<<"dcw, nan_cnt, err: "<<dcw<<", "<<nan_cnt<<", "<<err<<endl;
		
		blk_err2[i2]=err2/dcw2;		
		bit_err2[i2]=biterr2/(num2*dcw2);
		ncv_vec2[i2]=ncv2/(1*dcw2); 
		nvc_vec2[i2]=nvc2/(1*dcw2); 
		//if(i2==num_dat-1) {for(i=0;i<num2;i++) loc_err2[i]/=dcw2;}

		blk_err3[i2]=err3/dcw3;		
		bit_err3[i2]=biterr3/(num2*dcw3);
		ncv_vec3[i2]=ncv3/(1*dcw3); 
		nvc_vec3[i2]=nvc3/(1*dcw3); 

		blk_err4[i2]=err4/dcw4;		
		bit_err4[i2]=biterr4/(num2*dcw4);
		ncv_vec4[i2]=ncv4/(1*dcw4);
		nvc_vec4[i2]=nvc4/(1*dcw4); 

		blk_err5[i2]=err5/dcw5;		
		bit_err5[i2]=biterr5/(num2*dcw5);
		//cout<<"bit_err5: "<<bit_err5[i2]<<'\n'; 
		
		ncv_vec5[i2]=ncv5/(1*dcw5); 

		blk_err6[i2]=err6/dcw6;		
		bit_err6[i2]=biterr6/(num2*dcw6);
		ncv_vec6[i2]=ncv6/(1*dcw6); 	

		blk_err7[i2]=err7/dcw7;		
		bit_err7[i2]=biterr7/(num2*dcw7);
		ncv_vec7[i2]=ncv7/(1*dcw7); 
		//if(i2==num_dat-1) {for(i=0;i<num2;i++) loc_err3[i]/=dcw7;}
		//cout<<'\n'<<"ncv7: "<<ncv7<<" dcw7: "<<dcw7<<endl;
		
		blk_err8[i2]=err8/dcw8;		
		bit_err8[i2]=biterr8/(num2*dcw8);
		ncv_vec8[i2]=ncv8/(1*dcw8); 
		//if(i2==num_dat-1) {for(i=0;i<num2;i++) loc_err4[i]/=dcw8;}
		//cout<<'\n'<<"ncv8: "<<ncv8<<" dcw8: "<<dcw8<<endl;	

		num_iters1[i2]/=dcw; 
		num_iters2[i2]/=dcw2; 
		num_iters3[i2]/=dcw3; 
		num_iters4[i2]/=dcw4; 
		num_iters5[i2]/=dcw5;
		num_iters6[i2]/=dcw6; 
		num_iters7[i2]/=dcw7; 
		num_iters8[i2]/=dcw8;

		if(m1) {cout<<'\n'<<"err: "<<err<<" dcw: "<<dcw<<" nan_cnt: "<<nan_cnt<<", "<<'\n';}
		else if(m2) {cout<<'\n'<<"err2: "<<err2<<" dcw2: "<<dcw2<<" nan_cnt: "<<nan_cnt<<'\n';}
		else if(m3) {cout<<'\n'<<"err3: "<<err3<<" dcw2: "<<dcw3<<" nan_cnt: "<<nan_cnt<<'\n';}
		else if(m4) {cout<<'\n'<<"err4: "<<err4<<" dcw2: "<<dcw4<<" nan_cnt: "<<nan_cnt<<'\n';}
		else if(m5) {cout<<'\n'<<"err5: "<<err5<<" dcw2: "<<dcw5<<" nan_cnt: "<<nan_cnt<<'\n';}
		else if(m6) {cout<<'\n'<<"err6: "<<err6<<" dcw2: "<<dcw6<<" nan_cnt: "<<nan_cnt<<'\n';}
		else if(m7) {cout<<'\n'<<"err7: "<<err7<<" dcw2: "<<dcw7<<" nan_cnt: "<<nan_cnt<<'\n';}
		else if(m8) {cout<<'\n'<<"err8: "<<err8<<" dcw2: "<<dcw8<<" nan_cnt: "<<nan_cnt<<'\n';}
	
		if(m1) {cout<<"blk_err: "; for(i=0;i<num_dat;i++) cout<<blk_err[i]<< " "; cout<<endl;}
		if(m2) {cout<<"blk_err2: "; for(i=0;i<num_dat;i++) cout<<blk_err2[i]<< " "; cout<<endl;}
		if(m3) {cout<<"blk_err3: "; for(i=0;i<num_dat;i++) cout<<blk_err3[i]<< " "; cout<<endl;}
		if(m4) {cout<<"blk_err4: "; for(i=0;i<num_dat;i++) cout<<blk_err4[i]<< " "; cout<<endl;}
		if(m5) {cout<<"blk_err5: "; for(i=0;i<num_dat;i++) cout<<blk_err5[i]<< " "; cout<<endl;}
		if(m6) {cout<<"blk_err6: "; for(i=0;i<num_dat;i++) cout<<blk_err6[i]<< " "; cout<<endl;}
		if(m7) {cout<<"blk_err7: "; for(i=0;i<num_dat;i++) cout<<blk_err7[i]<< " "; cout<<endl;}
		if(m8) {cout<<"blk_err8: "; for(i=0;i<num_dat;i++) cout<<blk_err8[i]<< " "; cout<<endl;}	
	
		if(m1) {cout<<"bit_err: "; for(i=0;i<num_dat;i++) cout<<bit_err[i]<< " "; cout<<endl;}
		if(m2) {cout<<"bit_err2: "; for(i=0;i<num_dat;i++) cout<<bit_err2[i]<< " "; cout<<endl;}
		if(m3) {cout<<"bit_err3: "; for(i=0;i<num_dat;i++) cout<<bit_err3[i]<< " "; cout<<endl;}
		if(m4) {cout<<"bit_err4: "; for(i=0;i<num_dat;i++) cout<<bit_err4[i]<< " "; cout<<endl;}
		if(m5) {cout<<"bit_err5: "; for(i=0;i<num_dat;i++) cout<<bit_err5[i]<< " "; cout<<endl;}
		if(m6) {cout<<"bit_err6: "; for(i=0;i<num_dat;i++) cout<<bit_err6[i]<< " "; cout<<endl;}
		if(m7) {cout<<"bit_err7: "; for(i=0;i<num_dat;i++) cout<<bit_err7[i]<< " "; cout<<endl;}
		if(m8) {cout<<"bit_err8: "; for(i=0;i<num_dat;i++) cout<<bit_err8[i]<< " "; cout<<endl;}
		
		if(m1) {cout<<"ncv_vec: "; for(i=0;i<num_dat;i++) cout<<ncv_vec[i]<< " "; cout<<endl;}
		if(m2) {cout<<"ncv_vec2: "; for(i=0;i<num_dat;i++) cout<<ncv_vec2[i]<< " "; cout<<endl;}
		if(m3) {cout<<"ncv_vec3: "; for(i=0;i<num_dat;i++) cout<<ncv_vec3[i]<< " "; cout<<endl;}
		if(m4) {cout<<"ncv_vec4: "; for(i=0;i<num_dat;i++) cout<<ncv_vec4[i]<< " "; cout<<endl;}
		if(m5) {cout<<"ncv_vec5: "; for(i=0;i<num_dat;i++) cout<<ncv_vec5[i]<< " "; cout<<endl;}
		if(m6) {cout<<"ncv_vec6: "; for(i=0;i<num_dat;i++) cout<<ncv_vec6[i]<< " "; cout<<endl;}
		if(m7) {cout<<"ncv_vec7: "; for(i=0;i<num_dat;i++) cout<<ncv_vec7[i]<< " "; cout<<endl;}
		if(m8) {cout<<"ncv_vec8: "; for(i=0;i<num_dat;i++) cout<<ncv_vec8[i]<< " "; cout<<endl;}
		
		if(m1) {cout<<"nvc_vec: "; for(i=0;i<num_dat;i++) cout<<nvc_vec[i]<< " "; cout<<endl;}
		if(m2) {cout<<"nvc_vec2: "; for(i=0;i<num_dat;i++) cout<<nvc_vec2[i]<< " "; cout<<endl;}
		if(m3) {cout<<"nvc_vec3: "; for(i=0;i<num_dat;i++) cout<<nvc_vec3[i]<< " "; cout<<endl;}
		if(m4) {cout<<"nvc_vec4: "; for(i=0;i<num_dat;i++) cout<<nvc_vec4[i]<< " "; cout<<endl;}
		
		if(m1) {cout<<"num_iters1: "; for(i=0;i<num_dat;i++) cout<<num_iters1[i]<< " "; cout<<endl;}
		if(m2) {cout<<"num_iters2: "; for(i=0;i<num_dat;i++) cout<<num_iters2[i]<< " "; cout<<endl;}
		if(m3) {cout<<"num_iters3: "; for(i=0;i<num_dat;i++) cout<<num_iters3[i]<< " "; cout<<endl;}
		if(m4) {cout<<"num_iters4: "; for(i=0;i<num_dat;i++) cout<<num_iters4[i]<< " "; cout<<endl;}
		if(m5) {cout<<"num_iters5: "; for(i=0;i<num_dat;i++) cout<<num_iters5[i]<< " "; cout<<endl;}
		if(m6) {cout<<"num_iters6: "; for(i=0;i<num_dat;i++) cout<<num_iters6[i]<< " "; cout<<endl;}
		if(m7) {cout<<"num_iters7: "; for(i=0;i<num_dat;i++) cout<<num_iters7[i]<< " "; cout<<endl;}
		if(m8) {cout<<"num_iters8: "; for(i=0;i<num_dat;i++) cout<<num_iters8[i]<< " "; cout<<endl;}
	}
		
							
	for(meth=1;meth<=8;meth++) {
		for(i=0;i<num_dat;i++) {	 
			if(meth==1 && m1) {outf1.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf1<<blk_err[i]<< " "; outf1.close();} 
			else if(meth==2 && m2) {outf2.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf2<<blk_err2[i]<< " "; outf2.close();}
			else if(meth==3 && m3) {outf3.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf3<<blk_err3[i]<< " "; outf3.close();} 
			else if(meth==4 && m4) {outf4.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf4<<blk_err4[i]<< " "; outf4.close();} 
			else if(meth==5 && m5) {outf5.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf5<<blk_err5[i]<< " "; outf5.close();} 
			else if(meth==6 && m6) {outf6.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf6<<blk_err6[i]<< " "; outf6.close();} 
			else if(meth==7 && m7) {outf7.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf7<<blk_err7[i]<< " "; outf7.close();} 
			else if(meth==8 && m8) {outf8.open("FER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf8<<blk_err8[i]<< " "; outf8.close();}
			
			if(meth==1 && m1) {outf1.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf1<<bit_err[i]<< " "; outf1.close();} 
			else if(meth==2 && m2) {outf2.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf2<<bit_err2[i]<< " "; outf2.close();}
			else if(meth==3 && m3) {outf3.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf3<<bit_err3[i]<< " "; outf3.close();} 
			else if(meth==4 && m4) {outf4.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf4<<bit_err4[i]<< " "; outf4.close();} 
			else if(meth==5 && m5) {outf5.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf5<<bit_err5[i]<< " "; outf5.close();} 
			else if(meth==6 && m6) {outf6.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf6<<bit_err6[i]<< " "; outf6.close();} 
			else if(meth==7 && m7) {outf7.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf7<<bit_err7[i]<< " "; outf7.close();} 
			else if(meth==8 && m8) {outf8.open("BER_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf8<<bit_err8[i]<< " "; outf8.close();}
						
			if(meth==1 && m1) {outf1.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf1<<ncv_vec[i]<< " "; outf1.close();} 
			else if(meth==2 && m2) {outf2.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf2<<ncv_vec2[i]<< " "; outf2.close();}
			else if(meth==3 && m3) {outf3.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf3<<ncv_vec3[i]<< " "; outf3.close();} 
			else if(meth==4 && m4) {outf4.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf4<<ncv_vec4[i]<< " "; outf4.close();} 
			else if(meth==5 && m5) {outf5.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf5<<ncv_vec5[i]<< " "; outf5.close();} 
			else if(meth==6 && m6) {outf6.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf6<<ncv_vec6[i]<< " "; outf6.close();} 
			else if(meth==7 && m7) {outf7.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf7<<ncv_vec7[i]<< " "; outf7.close();}
			else if(meth==8 && m8) {outf8.open("ncv_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf8<<ncv_vec8[i]<< " "; outf8.close();} 
			
			if(meth==1 && m1) {outf1.open("nvc_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf1<<nvc_vec[i]<< " "; outf1.close();} 
			else if(meth==2 && m2) {outf2.open("nvc_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf2<<nvc_vec2[i]<< " "; outf2.close();}
			else if(meth==3 && m3) {outf3.open("nvc_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf3<<nvc_vec3[i]<< " "; outf3.close();} 
			else if(meth==4 && m4) {outf4.open("nvc_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf4<<nvc_vec4[i]<< " "; outf4.close();}
						
			if(meth==1 && m1) {outf1.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf1<<num_iters1[i]<< " "; outf1.close();} 
			else if(meth==2 && m2) {outf2.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf2<<num_iters2[i]<< " "; outf2.close();}
			else if(meth==3 && m3) {outf3.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf3<<num_iters3[i]<< " "; outf3.close();} 
			else if(meth==4 && m4) {outf4.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf4<<num_iters4[i]<< " "; outf4.close();} 
			else if(meth==5 && m5) {outf5.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf5<<num_iters5[i]<< " "; outf5.close();} 
			else if(meth==6 && m6) {outf6.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf6<<num_iters6[i]<< " "; outf6.close();} 
			else if(meth==7 && m7) {outf7.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf7<<num_iters7[i]<< " "; outf7.close();}
			else if(meth==8 && m8) {outf8.open("niter_"+func(fn)+"_"+func(meth)+".txt", std::ios_base::app); outf8<<num_iters8[i]<< " "; outf8.close();} 
		}
	}
	
	outf9.open("snr"+func(fn)+".txt"); 
	for(i=0;i<num_dat;i++) 
		outf9<<EbNo_BC[i]<<" "; 
	outf9.close();
	
	cout<<'\n'<<"fn, R: "<<fn<<", "<<R<<endl;
			
	/*delete[] x_hat;
	delete[] y;  
	delete[] LR;  
	delete[] LR_srtd;
	delete[] pLR;  
	delete[] syn;
	delete[] syn_cls;
	
	delete[] excl_cw;
	delete[] indx;
	delete[] indx2;
	delete[] rep;
	
	delete[] E_v_c;
	delete[] E_c_v;
	delete[] E_v_c_mu;
	delete[] E_c_v_mu;
	delete[] res_c_v;
	delete[] res_c_v_srtd;
	delete[] res_v_c;
	delete[] res_v_c_srtd;
	
	delete[] cns; 
	delete[] Q;
	delete[] Q_cnt;
	delete[] Q_temp;
	
	delete[] blk_err; 
	delete[] blk_err2;
	delete[] blk_err3; 
	delete[] blk_err4; 
	delete[] blk_err5;
	delete[] blk_err6;
	delete[] blk_err7;
	delete[] blk_err8;

	delete[] bit_err; 
	delete[] bit_err2; 	
	delete[] bit_err3; 
	delete[] bit_err4; 
	delete[] bit_err5;
	delete[] bit_err6;
	delete[] bit_err7;
	delete[] bit_err8;

	delete[] ncv_vec;
	delete[] ncv_vec2;
	delete[] ncv_vec3;
	delete[] ncv_vec4;
	delete[] ncv_vec5;
	delete[] ncv_vec6;
	delete[] ncv_vec7;
	delete[] ncv_vec8;
	
	delete[] vns;
	delete[] cns;
	delete[] excl_cw;
	delete[] indx3;
	delete[] vn_indx;
	delete[] cn_indx;
	delete[] indx_cls_CN;
	delete[] indx_cls;
	delete[] gcn_indx;
	delete[] gcn_indx_ordr;
	delete[] spc_indx;
	delete[] spc_indx_ordr;
	delete[] cls_s_crnt;
	delete[] cls_s_prev;
	delete[] pick_cls;
	delete[] syn;
	delete[] syn_cls;
	
	delete[] mu_cnt;
	delete[] mu_cnt_L;
	delete[] perm;
	delete[] E_c_v_cnt;
	delete[] num_iters1;
	delete[] num_iters2;
	delete[] num_iters3;
	delete[] num_iters4;
	delete[] num_iters5;
	delete[] num_iters6;
	delete[] num_iters7;
	delete[] num_iters8;
	delete[] inp_vec;*/
	
	//delete[] j_sv;
	
	//**cls_ind_mat_contg,**cls_ind_mat_ran,**cls_ind_mat_opt,**cns_subcode,**vns_subcode,**cns_cluster,**vns_cluster
	cout<<'\n';
    printf("executed in: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	cout<<'\n';

}





