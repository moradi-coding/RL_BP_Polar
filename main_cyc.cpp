//just counts cycles from the input H matrix

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

long par_num=0,chld_num=0,tcyc=5*pow(10,6),revs=5*pow(10,6),abs_64=0; //max. no. of reverses in a rnd can be 12
long *snk_r,*snk_c,*E,**E_vec,**fbd_r,**fbd_c,**H,**Hsc,**Hsc2,**H0,**H1,**Hsc_2,**Hsc_3,**I,**I2,**sigma,**sigma2,**sigma3,**tau2,**B,
*par_r,*par_c,**chld_r,**chld_c,*col_wt,**cycol,
**cyc_sv_r,**cyc_sv_c,**cyc_sv_rg,**cyc_sv_cg,*cyc_sel_r,*cyc_sel_c,**cyc_param,**cyc_cg,**cyc_r_nel,**cyc_c_nel,**regn; 
long p,q,gama,J=1,m,n,L,cyc_len,strt_r,strt_c,stp_r,stp_c,row_strt,col_strt,row_stp,col_stp; //H has total m rows and n cols
long while_tot=0,rev_tot=0,tot_r,tot_c,sum=0,ccnt=0; 


#include "find_cycl4.cpp"

string func(long n) 
{
    stringstream result;
    result << n;
    return result.str();
}

int main() {	
	
	srand(time(0));	
	double t2,t1= time(NULL);	
	long i,i5,i2,i3,i4,j,j2,k,l,count=0,A,b2,*outp,fn; 
	long f,fl=1,totf=1,*edge_r,*edge_c,num_edgs=5*pow(10,6),param,ccnt2=0,flg=0,rr,cnt,ncyc=0,c,BC,*sset_r,*sset_c,mu1,mu2,mu3,*clusters,*cn_cnt,*indx,tmp;   
	
	ifstream inf,inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf10; 
	cout<<'\n'<<"file no: "; cin>>fn; 

	cyc_len=4; 

	if(fn==11) {m=292; n=2190; inf.open("matrices/mat_base_liva.txt"); col_stp=n;} 
	else if(fn==13) {m=76; n=108; inf.open("matrices/mat_76_108_polar.txt"); col_stp=n;} 
	else if(fn==14) {m=104; n=136; inf.open("matrices/mat_104_136_polar.txt"); col_stp=n;} 
	else if(fn==16) {m=85; n=117; inf.open("matrices/mat_85_117_polar.txt"); col_stp=n;} 
	else if(fn==17) {m=187; n=251; inf.open("matrices/mat_187_251_polar.txt"); col_stp=n;} 
	else if(fn==18) {m=77; n=109; inf.open("matrices/mat_77_109_polar.txt"); col_stp=n;} 
	else if(fn==19) {m=84; n=116; inf.open("matrices/mat_84_116_polar.txt"); col_stp=n;} 
	else if(fn==20) {m=191; n=255; inf.open("matrices/mat_191_255_polar.txt"); col_stp=n;} 
	else if(fn==21) {m=191; n=255; inf.open("matrices/mat_191_255_polar2.txt"); col_stp=n;} 
	else if(fn==23) {m=79; n=111; inf.open("matrices/mat_79_111_polar.txt"); col_stp=n;} 
	else if(fn==25) {m=83; n=115; inf.open("matrices/mat_83_115_polar.txt"); col_stp=n;} 

	cout<<'\n'<<"m: "<<m<<" n: "<<n<<" fn: "<<fn<<endl;
	tot_r=m; tot_c=n; 

	edge_r=new long[num_edgs]; edge_c=new long[num_edgs]; outp=new long[12]; 
	par_r=new long[cyc_len]; par_c=new long[cyc_len]; clusters=new long[m]; 
	snk_r= new long[cyc_len]; snk_c= new long[cyc_len]; 
	chld_r=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) chld_r[i4]=new long[cyc_len]; 
	chld_c=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) chld_c[i4]=new long[cyc_len]; 
	cyc_sv_r=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_r[i4]=new long[cyc_len]; 
	cyc_sv_c=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_c[i4]=new long[cyc_len];
	cycol=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cycol[i4]=new long[cyc_len/2];
	cyc_sv_rg=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_rg[i4]=new long[cyc_len]; 
	cyc_sv_cg=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_sv_cg[i4]=new long[cyc_len];
	cyc_cg=new long*[tcyc]; for(i4=0;i4<tcyc;i4++) cyc_cg[i4]=new long[cyc_len];
	col_wt=new long[tot_c];
	cn_cnt=new long[m];
	indx=new long[m];
	fbd_r=new long*[revs]; for(i5=0;i5<revs;i5++) fbd_r[i5]=new long[cyc_len]; 
	fbd_c=new long*[revs]; for(i5=0;i5<revs;i5++) fbd_c[i5]=new long[cyc_len]; 
	for(i2=0;i2<revs;i2++) for(i3=0;i3<cyc_len;i3++) fbd_r[i2][i3]=fbd_c[i2][i3]=-1; //initializing

	regn=new long*[gama]; for(i=0;i<gama;i++) regn[i]=new long[p];

	H=new long*[m]; for(i5=0;i5<m;i5++) H[i5]=new long[n]; //m rows n cols
	Hsc=new long*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc[i5]=new long[tot_c];
	//Hsc2=new long*[tot_r]; for(i5=0;i5<tot_r;i5++) Hsc2[i5]=new long[tot_c];

	ofstream outf,outf2,outf3,outf4,outf5,outf6; string filename; 
	
	for(i=0;i<m;i++) for(j=0;j<n;j++) inf >> Hsc[i][j];

	row_strt=col_strt=0; row_stp=tot_r; //col_stp=tot_c; 
	//cout<<'\n'<<"Hsc: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) cout<<Hsc[i3][j]<<""; cout<<'\n';} 
	/*cout<<'\n'<<"H1: "<<'\n'; 
	for(i=0;i<m;i++) {
		for(j=0;j<n;j++) {
			cout<<Hsc[i][j]<<""; 
			if(j>0 && (j+1)%(6)==0) cout<<" ";
		} 
		cout<<'\n'; if(i>0 && (i+1)%(6)==0) cout<<'\n';
	}*/

	//cout<<'\n'<<"row_stp "<<row_stp<<" col_stp "<<col_stp<<" cyc_len "<<cyc_len<<std::endl;

	sum=count=ccnt=0; for(i2=0;i2<tot_c;i2++) col_wt[i2]=0; //refresh
	for(i=row_strt;i<row_stp;i++) 
		for(i2=col_strt;i2<col_stp;i2++) 
			if(Hsc[i][i2]) {
				edge_r[count]=i; 
				edge_c[count]=i2; 
				count++;
			} 

	//**************imp for standard dec. testing, col_wt=3 always**********************/
	for(i2=col_strt;i2<col_stp;i2++) 
		for(i=row_strt;i<row_stp;i++) 
			if(Hsc[i][i2]) col_wt[i2]=3/*++*/; 
	//cout<<'\n'<<"count: "<<count<<std::endl;
	if(totf>1) {
		A=(fl-1)*(count/(totf-1)); 
		if(fl<totf) 
			b2=A+count/(totf-1); 
		else b2=count;
	}
	else {
		A=0; 
		b2=count;
	}
	for(i=0;i<A;i++) Hsc[edge_r[i]][edge_c[i]]=0; //decimating edges used by other processors (if any)
	cout<<'\n'<<"A: "<<A<<" b2-1: "<<b2-1<<std::endl;

	outf.open("cycle_vns.txt", std::ios_base::app);

	for(i=A;i<b2;i++) { 
		//while_tot=0
		par_num=chld_num=0; for(i3=0;i3<cyc_len;i3++) snk_r[i3]=snk_c[i3]=par_r[i3]=par_c[i3]=-1; //refresh
		strt_r=edge_r[i]; strt_c=edge_c[i];
		//cout<<'\n'<<"start: "<<strt_r<<","<<strt_c;
		//cout<<"i: "<<i<<endl;	
		find_cycl();
		Hsc[strt_r][strt_c]=0; //decimating unnecessary edges
		//cout<<'\n'<<"while_tot: "<<while_tot<<" rev_tot: "<<rev_tot;
		sum+=par_num+chld_num;
		//cout<<"ccnt: "<<ccnt<<endl;

		//for(j=0;j<ccnt;j++) 
			/*for(i3=0;i3<cyc_len;i3++) 
				if(!(i3%2)) 
					outf<<cyc_sv_c[ccnt][i3]<<" "; */

		//if(count==30) break;
	}

	outf<<std::endl; 
	outf.close();

	cout<<'\n'<<"no. of "<<cyc_len<<" cycles: "<<ccnt<<" fn: "<<fn<<endl;
	row_strt=row_stp=col_strt=col_stp=0;

	//cout<<'\n'; for(j=0;j<ccnt;j++) {for(i3=0;i3<cyc_len;i3++) cout<<cyc_sv_r[j][i3]<<","<<cyc_sv_c[j][i3]<<" "; cout<<'\n';} //all cycles

	//save set of VNs per cycle

	/*cout<<'\n'; 
	for(j=0;j<ccnt;j++) {
		for(i3=0;i3<cyc_len;i3++) 
			if(!(i3%2)) 
				cout<<cyc_sv_c[j][i3]<<" "; 
			cout<<'\n';
	}*/

	/*outf.open("cycle_vns.txt"); //, std::ios_base::app
	for(j=0;j<ccnt;j++) {
		for(i3=0;i3<cyc_len;i3++) 
			if(!(i3%2)) 
				outf<<cyc_sv_c[j][i3]<<" "; 
	}
	outf<<std::endl; 
	outf.close();*/

	//cout<<'\n'<<"indx: "; for(i=0;i<m;i++) cout<<indx[i]<<" "; 

	//filename="cg"+func(fl)+".txt"; 
	//outf4.open("cls_ind_mat_opt_AB_3.txt"); for(i=0;i<m;i++) outf4<<clusters[i]<<" "; outf4<<std::endl; outf4.close();
	//then apply read2 to get trg and tcg etc

	t2= time(NULL); cout<<'\n'<<"Executed in "<<t2-t1<<"s"<<'\n';
	cout<<'\n';
}
