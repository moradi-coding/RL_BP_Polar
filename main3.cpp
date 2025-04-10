//creat GLDPC matrix from a base matrix

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

int **Hbase,**hamming,**Hlift,**Hsc_sv,**Hsc_proto,**H0,**H1,**I,**I2,**I3,**B,**sigma,**sigma2,**sigma3,**tau2,**Hsc_2,**Hsc_3,**lsft_mat,*E,**cns_subcode,**vns_subcode; 
int m,n,opt,randm,m_sub,n_sub,row_wt;  //H has total m rows and n cols

//params for find_cycl4
long par_num=0,chld_num=0,tcyc=pow(10,5),revs=pow(10,5),abs_64=0; //max. no. of reverses in a rnd can be 12
int *snk_r,*snk_c,**fbd_r,**fbd_c,*par_r,*par_c,**chld_r,**chld_c,*col_wt,**cyc_sv_r,**cyc_sv_c,**cyc_sv_rg,**cyc_sv_cg,*cyc_sel_r,*cyc_sel_c,**cyc_param,**cyc_cg,**cyc_r_nel,**cyc_c_nel; 
int strt_r,strt_c,stp_r,stp_c,row_strt,col_strt,row_stp,col_stp,*edge_r,*edge_c; //H has total m rows and n cols
long while_tot=0,rev_tot=0,tot_r,tot_c,sum=0,ccnt=0,ccnt_old=0,num_edgs=pow(10,5),ge=1,mem,s,cyc_len=6; 

#include "H_GLDPC.cpp"

string func(int n) 
{
    stringstream result;
    result << n;
    return result.str();
}

int main() {	
	
	srand(time(0));	
	double t2,t1= time(NULL);	
	long i,i5,i2,i3,i4,j,j2,k,l,count=0,A,*outp; 
	int f=1,totf=1,param,ccnt2=0,flg=0,rr,cnt,ncyc=0,Allison=1,unc=0,blks,ers_r=-1,rf,glob=0;   
	
	ifstream inf4,inf5,inf6,inf7,inf10,inf11,inf12; 

	/*cout<<'\n'<<" f: "; cin>>f; 
	if(f==1){
		cout<<'\n'<<" p: "; 
		cin>>p; cout<<'\n'<<" gama: "; 
		cin>>gama; cout<<'\n'<<" sub-code: "; 
		cin>>s; cout<<'\n'<<" J: "; 
		cin>>J; 
		L=1;
	} 
	else {inf6.open("p.txt"); inf6>>p; inf10.open("gama.txt"); inf10>>gama; inf12.open("s.txt"); inf12>>s;}*/

	//dimension of base matrix , m is also the no. of sub-codes, each with m_sub CNs
	m=40; n=300;	
	//m=68; n=510;
	//m=132; n=990;
	
	//dimensions of sub-code parity-check matrix
	m_sub=4; n_sub=15;
	row_wt=15;
	
	//dimension of matrix lifted with sub-code matrix column
	tot_r=m*m_sub; tot_c=n;
	
	cout<<'\n'<<"tot_r "<<tot_r; 
	cout<<'\n'<<"tot_c "<<tot_c; 

	Hbase=new int*[m]; for(i5=0;i5<m;i5++) Hbase[i5]=new int[n]; //m rows n cols
	hamming=new int*[m_sub]; for(i5=0;i5<m_sub;i5++) hamming[i5]=new int[n_sub]; 
	Hlift=new int*[tot_r]; for(i5=0;i5<tot_r;i5++) Hlift[i5]=new int[tot_c]; 
	cns_subcode=new int*[m]; for(i5=0;i5<m;i5++) cns_subcode[i5]=new int[m_sub];
	vns_subcode=new int*[m]; for(i5=0;i5<m;i5++) vns_subcode[i5]=new int[row_wt];
	
	//for(i3=0;i3<m;i3++) for(j=0;j<row_wt;j++) vns_subcode[i3][j]=-1;
	
	inf4.open("matrices/hamming2.txt"); for(j=0;j<m_sub;j++) for(i=0;i<n_sub;i++) inf4 >> hamming[j][i];
	inf5.open("mat.txt"); for(j=0;j<m;j++) for(i=0;i<n;i++) inf5 >> Hbase[j][i];

	ifstream inf1,inf2,inf3;
	ofstream outf,outf2,outf3,outf4,outf5,outf6,outf7,outf8,outf9,outf10,outf11,outf12,outf13,outf14,outf15,outf16; string filename; 

	//cout<<'\n'<<"hamming: "<<'\n'; for(i3=0;i3<m_sub;i3++) {for(j=0;j<n_sub;j++) cout<<hamming[i3][j]<<" "; cout<<'\n'; }
	//cout<<'\n'<<"Hbase: "<<'\n'; for(i3=0;i3<m;i3++) {for(j=0;j<n;j++) cout<<Hbase[i3][j]<<" "; cout<<'\n'; }

	H_GLDPC();
	
	//cout<<'\n'<<"Hlift: "<<'\n'; for(i3=0;i3<tot_r;i3++) {for(j=0;j<tot_c;j++) cout<<Hlift[i3][j]<<" "; cout<<'\n';}
	cout<<'\n'<<"cns_subcode: "<<'\n'; for(i3=0;i3<m;i3++) {for(j=0;j<m_sub;j++) cout<<cns_subcode[i3][j]<<" "; cout<<'\n';}
	cout<<'\n'<<"vns_subcode: "<<'\n'; for(i3=0;i3<m;i3++) {for(j=0;j<row_wt;j++) cout<<vns_subcode[i3][j]<<" "; cout<<'\n';}
	
	outf.open("matrices/mat_gldpc.txt"); for(i=0;i<tot_r;i++) for(j=0;j<tot_c;j++) outf<<Hlift[i][j]<<" "; outf<<std::endl; outf.close();
	outf2.open("matrices/gldpc_cns_subcode.txt"); for(i=0;i<m;i++) for(j=0;j<m_sub;j++) outf2<<cns_subcode[i][j]<<" "; outf2<<std::endl; outf2.close();
	outf3.open("matrices/gldpc_vns_subcode.txt"); for(i=0;i<m;i++) for(j=0;j<row_wt;j++) outf3<<vns_subcode[i][j]<<" "; outf3<<std::endl; outf3.close();
	
	t2= time(NULL); cout<<'\n'<<"Executed in "<<t2-t1<<"s"<<'\n';
	cout<<'\n';
}
