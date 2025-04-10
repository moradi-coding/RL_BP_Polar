//tries to generate a random regular LDPC codes without 4-cycles with col. wt. gama and row. wt p

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <sstream>
using namespace std;

int **H_proto,**H,**Hsc,**Hsc_sv,**Hsc_proto,**H0,**H1,**I,**I2,**I3,**B,**sigma,**sigma2,**sigma3,**tau2,**Hsc_2,**Hsc_3,**lsft_mat,*E; 
int gama,p,q,m2,n2,J=1,m,n,L=4,opt=0,randm; //H has total m rows and n cols

//params for find_cycl4
long par_num=0,chld_num=0,tcyc=pow(10,6),revs=pow(10,6),abs_64=0; //max. no. of reverses in a rnd can be 12
long *snk_r,*snk_c,**fbd_r,**fbd_c,*par_r,*par_c,**chld_r,**chld_c,**cyc_sv_r,**cyc_sv_c,**cyc_sv_rg,**cyc_sv_cg,*cyc_sel_r,*cyc_sel_c,**cyc_param,**cyc_cg,**cyc_r_nel,
**cyc_c_nel; 
long strt_r,strt_c,stp_r,stp_c,row_strt,col_strt,row_stp,col_stp,*edge_r,*edge_c; //H has total m rows and n cols
long while_tot=0,rev_tot=0,tot_r,tot_c,sum=0,ccnt=0,ccnt_old=0,num_edgs=pow(10,6),ge=1,mem=1,s,cyc_len=4,**cycol; 

//#include "H_matb.cpp"
//#include "main_cyc4.cpp"
#include "find_cycl4.cpp"

string func(int n) {
    stringstream result;
    result << n;
    return result.str();
}

int main() {	
	srand(time(0));	
	double t2,t1=time(NULL),cnt3,cnt3_old;	
	long i,i5,i2,i3,i4,j,j2,k,l,count=0,A,*outp,*k_sv; 
	int fn=0,totf=1,param,ccnt2=0,flg=0,flg2=0,flg3=0,rr,ass,ncyc=0,BC=0,blks,ers_r=-1,rf,glob=0,meth=1,col_wt,row_wt,btrck=0,
	rnd,rnd_max=200,iter_max=1,iter,i_btrck,cnt,cnt2;   
	int cyc_chck=1; //1 checks for 4-cycles, 0 does not
	
	ifstream inf4,inf5,inf6,inf7,inf10,inf11,inf12; 
	//cout<<'\n'<<" f: "; cin>>f; 
	m2=50; n2=100; fn=3; gama=3; p=6;
	//m2=21; n2=286; fn=3; gama=5; p=69; //try to keep p as close as possible to theretical value for a given gama and rate
	//m2=134*2; n2=469*2; fn=1; gama=2; p=7;
	//m2=4; n2=14; fn=1; gama=2; p=7;
	//m2=4; n2=30; fn=1; gama=2; p=15;
	//m2=15; n2=30; fn=1; gama=3; p=6;
	
	//if(f==1 && !s) {cout<<'\n'<<" rg erase: "; cin>>ers_r;} //if -1, then no row groups to erase
	//inf11.open("ge.txt"); inf11>>ge; //general edge spreading if 1, not otherwise

	m=m2; n=n2;
	tot_r=m; tot_c=n;
	
	E=new int[gama]; 
	H_proto=new int*[m2]; for(i5=0;i5<m2;i5++) H_proto[i5]=new int[n2]; //m rows n cols
	H=new int*[m]; for(i5=0;i5<m;i5++) H[i5]=new int[n]; //m rows n cols
	H0=new int*[m]; for(i5=0;i5<m;i5++) H0[i5]=new int[n]; H1=new int*[m]; for(i5=0;i5<m;i5++) H1[i5]=new int[n]; 
	Hsc_proto=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc_proto[i5]=new int[tot_c]; 
	Hsc=new int*[m]; for(i5=0;i5<m;i5++) Hsc[i5]=new int[n]; 
	Hsc_sv=new int*[m]; for(i5=0;i5<m;i5++) Hsc_sv[i5]=new int[n]; 
	I=new int*[q]; for(i5=0;i5<q;i5++) I[i5]=new int[q];
	I2=new int*[L]; for(i5=0;i5<L;i5++) I2[i5]=new int[L];
	tau2=new int*[L]; for(i5=0;i5<L;i5++) tau2[i5]=new int[L];
	sigma2=new int*[L]; for(i4=0;i4<L;i4++) sigma2[i4]=new int[L];
	//B=new int*[gama]; for(i5=0;i5<gama;i5++) B[i5]=new int[p];
	sigma=new int*[q]; for(i5=0;i5<q;i5++) sigma[i5]=new int[q];
	sigma3=new int*[J]; for(i4=0;i4<J;i4++) sigma3[i4]=new int[J];
	lsft_mat=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) lsft_mat[i5]=new int[tot_c];

	edge_r=new long[num_edgs]; edge_c=new long[num_edgs]; outp=new long[12]; 
	par_r=new long[cyc_len]; par_c=new long[cyc_len]; 
	snk_r= new long[cyc_len]; snk_c= new long[cyc_len]; 
	chld_r=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) chld_r[i4]=new long[cyc_len]; 
	chld_c=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) chld_c[i4]=new long[cyc_len]; 
	cycol=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cycol[i4]=new long[cyc_len/2];
	cyc_sv_r=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_r[i4]=new long[cyc_len]; 
	cyc_sv_c=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_c[i4]=new long[cyc_len];
	cyc_sv_rg=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_rg[i4]=new long[cyc_len]; 
	cyc_sv_cg=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_cg[i4]=new long[cyc_len];
	cyc_cg=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_cg[i4]=new long[cyc_len];
	//col_wt=new long[tot_c];
	k_sv=new long[m];
	fbd_r=new long*[revs]; for(i5=0;i5<revs;i5++) fbd_r[i5]=new long[cyc_len]; 
	fbd_c=new long*[revs]; for(i5=0;i5<revs;i5++) fbd_c[i5]=new long[cyc_len]; 
	for(i2=0;i2<revs;i2++) for(i3=0;i3<cyc_len;i3++) fbd_r[i2][i3]=fbd_c[i2][i3]=-1; //initializing

	ifstream inf1,inf2,inf3;
	ofstream outf,outf2,outf3,outf4,outf5,outf6,outf7,outf8,outf9,outf10,outf11,outf12,outf13,outf14,outf15,outf16; string filename; 

	//circulant shift factors for sub-code 1. For random B matrix, circulant shifts are random

	cout<<'\n'<<"m "<<m<<" n "<<n;	
	row_strt=col_strt=0; row_stp=tot_r; col_stp=tot_c; //for cycle counting
 
iter=0;
cnt3=cnt3_old=1000;
//while(cnt3>m/2-1) {
while(iter<iter_max) {
	for(i=0;i<m;i++)	
		for(j=0;j<n;j++) 
			Hsc[i][j]=0;
	
	for(i=0;i<n;i++) {
		cout<<'\n'<<"i: "<<i; 					
		//cout<<'\n'<<"Hsc bef. btrck: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) {cout<<Hsc[i3][j]<<" "; } cout<<'\n';}
		col_wt=ass=rnd=0;
		//cyc_chck=1;	
		while(col_wt<gama) {
			flg=flg3=0;
			ass=ass%m;
			if(!ass || rnd>rnd_max) {
				col_wt=0;
				for(j=0;j<m;j++) {
					Hsc[j][i]=0; 
					k_sv[j]=-1; //refresh
				}
				if(rnd>rnd_max) {
					if(!btrck)
						i_btrck=i; //saving the col. when btrck starts
						
					btrck++;	
					for(i2=i;i2>i-btrck;i2--) 	
						for(j=0;j<m;j++) 
							Hsc[j][i]=0;
					if(i-btrck>0) 
						i-=btrck; //do backtracking by removing cols.
					else {
						btrck=0;
						i--;
						cnt++;
					}
					//cout<<'\n'<<"Hsc: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) {cout<<Hsc[i3][j]<<" "; } cout<<'\n';}	
					//cout<<'\n'<<"btrck, i_btrck, i: "<<btrck<<", "<<i_btrck<<", "<<i<<endl;
					//flg3=1; //backtracking is done
					if(cnt>1) 
						cyc_chck=0; //stop cycle check as it is not possible to remove cycles anymore
					break;
				}
			}
			
			while(ass>0 && ass<m && !flg) {
				k=rand()%m;
				//cout<<'\n'<<"k: "<<k;
				for(j=0;j<ass;j++) {
					if(k!=k_sv[j]) 
						flg=1; 
					else {
						flg=0; 
						break;
					}  
				}
				//cout<<'\n'<<"k_sv: "; for(j=0;j<ass;j++) cout<<k_sv[j]<<" ";
				//cout<<'\n'<<"ass: "<<ass;
			}
			if(!ass) {
				k=rand()%m;
				col_wt++;
			}
			
			if(ass<m) { 
				Hsc[k][i]=1;
				k_sv[ass]=k; 
				ass++; 
			}
			
			if(ass==m)
				rnd++;	
			
			//cout<<'\n'<<"ass: "<<ass; 
			//cout<<'\n'<<"k_sv: "; for(j=0;j<ass;j++) cout<<k_sv[j]<<" ";

			row_wt=flg=0;
			for(j=0;j<=i;j++) {		
				if(Hsc[k][j]) 
					row_wt++;
				if(row_wt>p) {
					flg=1;
					break;
				}
				else 
					flg=0;
			}
			
			if(flg) {//if row_wt>p due to a 1 assignment in col. i
				Hsc[k][i]=0;
				if(ass==1)
					col_wt--;
				//cout<<'\n'<<"col_wt after edge rem.: "<<col_wt<<endl;
			}
				
			else if(!flg && ass>1) { 
				//check if the new edge generates any 4-cycle
				//if(ass>1) { 
					par_num=chld_num=0; 
					for(i3=0;i3<cyc_len;i3++) 
						snk_r[i3]=snk_c[i3]=par_r[i3]=par_c[i3]=-1; //refresh
					strt_r=k; 
					strt_c=i;					
								
					ccnt=0;			
					if(cyc_chck) 
						find_cycl();
					
					if(cyc_chck && ccnt) { //if there is a 4-cycle
						//cout<<'\n'<<"start: "<<strt_r<<","<<strt_c;	
						//cout<<'\n'<<"no. of cycles: "<<ccnt<<endl; 	
						if(ass<m) 
							Hsc[k][i]=0;
						else {
							for(j=0;j<m;j++) {	
								Hsc[j][i]=0; //remove the column because no more 1 can be assigned to col. i
								k_sv[j]=-1; //refresh
							}
							ass=col_wt=0;	
						}	
					}
					else
						col_wt++; 
			}
			//cout<<'\n'<<"col_wt: "<<col_wt<<endl;
			
			if(col_wt==gama && i>i_btrck) 
				btrck=0;
			
			/*if(btrck>100) {
				cout<<"cannot assign!!!";
				flg2=1;
				break;
			}*/					
		}
	}	
	
	/*for(i=0;i<m;i++) 
		for(j=0;j<n;j++) 
			Hsc_sv[i][j]=Hsc[i][j];

	cyc_len=6;
	main_cyc4(cyc_len);
	cout<<'\n'<<"ccnt: "<<ccnt;*/
	
		cnt3=0;
		for(j=0;j<m;j++) {
			cnt2=0; 
			for(i=0;i<n;i++) 
				if(Hsc[j][i]) 
					cnt2++;
			if(cnt2==p)
				cnt3++; //no. of CNs with max degree
		}
		//cout<<'\n'<<"CNs with max deg: "<<cnt3<<endl;
		
	//output matrix 
	/*if(!iter || ccnt<ccnt_old) { 
		if(cnt3<cnt3_old) {
			//ccnt_old=ccnt;
			//cnt3_old=cnt3;
			outf4.open("mat_"+func(fn)+".txt");
			for(i=0;i<m;i++) {
				for(j=0;j<n;j++) 
					outf4<<Hsc[i][j]<<" ";  //Hsc_sv[i][j]
				if(i<m)
					outf4<<'\n';
			}
			outf4.close(); 
		}
	}*/
	//cout<<'\n'<<"Hsc: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) {cout<<Hsc[i3][j]<<" "; } cout<<'\n';}
	iter++;	
}

	//check validity of matrix
	for(j=0;j<n;j++) {
		ass=0; 
		for(i=0;i<m;i++) if(Hsc[i][j]) ass++; //if(Hsc_sv[i][j])
		if(ass!=gama) {flg=1; cout<<'\n'<<"col no.: "<<j; break;}
		else flg=0;
	}
	if(!flg) cout<<'\n'<<"valid col. wt."; 
	else cout<<'\n'<<"invalid col. wt: "<<ass;

	for(i=0;i<m;i++) {
		ass=0; 
		for(j=0;j<n;j++) if(Hsc[i][j]) ass++; //if(Hsc_sv[i][j])
		if(ass!=p) {flg=1; cout<<'\n'<<"row no.: "<<i; break;}
		else flg=0;
	}
	if(!flg) cout<<'\n'<<"valid row. wt."; 
	else cout<<'\n'<<"invalid row. wt.: "<<ass;
	//}
	
	outf4.open("mat_"+func(fn)+".txt");
			for(i=0;i<m;i++) {
				for(j=0;j<n;j++) 
					outf4<<Hsc[i][j]<<" ";  //Hsc_sv[i][j]
				if(i<m)
					outf4<<'\n';
			}
			outf4.close();
	
	if(!BC) {cout<<'\n'<<"tot_r "<<m<<" tot_c "<<n;}
	else {cout<<'\n'<<"tot_r "<<gama*p*J<<" tot_c "<<n;} 


	//cout<<'\n'<<"mval: "<<m<<" nval: "<<n;

	t2= time(NULL); cout<<'\n'<<"Executed in "<<t2-t1<<"s"<<'\n';
	cout<<'\n';
}
