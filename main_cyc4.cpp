//called only by main4.cpp

#include "find_cycl4.cpp"

int main_cyc4(int cyc_len) {

	long i,i5,i2,i3,i4,j,j2,k,l,count=0,A,b2; 
	long f,fl=1,totf=1,param,ccnt2=0,flg=0,rr,cnt,ncyc=0,c,BC;   
	
	//ifstream inf1,inf2,inf3,inf4,inf5,inf6,inf7,inf10; 
	//inf1.open("mval.txt"); inf1>>m; inf2.open("nval.txt"); inf2>>n; inf3.open("fval.txt"); inf3>>f; inf6.open("p.txt"); inf6>>p; inf10.open("gama.txt"); inf10>>gama; 
	for(i2=0;i2<revs;i2++) for(i3=0;i3<cyc_len;i3++) fbd_r[i2][i3]=fbd_c[i2][i3]=-1; //initializing

	row_strt=col_strt=0; row_stp=tot_r*J; col_stp=tot_c*J; 
	//cout<<'\n'<<"Hsc: "<<'\n'; for(i=row_strt;i<row_stp;i++) {for(j=col_strt;j<col_stp;j++) {cout<<Hsc[i][j]<<""; if(j>0 && (j+1)%p==0) cout<<" ";} cout<<'\n'; if(i>0 && (i+1)%p==0) cout<<'\n';}
	//cout<<'\n'<<"Hsc: "<<'\n'; for(i3=0;i3<m;i3++){for(j=0;j<n;j++) {cout<<Hsc[i3][j]<<""; if(j>0 && (j+1)%(p*J)==0) cout<<" ";} cout<<'\n'; if(i3>0 && (i3+1)%(p*J)==0) cout<<'\n';}

	//cout<<'\n'<<"row_stp "<<row_stp<<" col_stp "<<col_stp<<" cyc_len "<<cyc_len<<std::endl;
	sum=count=ccnt=0; 
	//for(i2=0;i2<tot_c;i2++) col_wt[i2]=0; //refresh
	for(i=row_strt;i<row_stp;i++) for(i2=col_strt;i2<col_stp;i2++) if(Hsc[i][i2]) {edge_r[count]=i; edge_c[count]=i2; count++;} 
	//**************imp for standard dec. testing, col_wt=3 always**********************/
	//for(i2=col_strt;i2<col_stp;i2++) for(i=row_strt;i<row_stp;i++) if(Hsc[i][i2]) col_wt[i2]=3; 
	//cout<<'\n'<<"count: "<<count<<std::endl;
	if(totf>1) {A=(fl-1)*(count/(totf-1)); if(fl<totf) b2=A+count/(totf-1); else b2=count;}
	else {A=0; b2=count;}
	for(i=0;i<A;i++) Hsc[edge_r[i]][edge_c[i]]=0; //decimating unnecessary edges
	//cout<<'\n'<<"A: "<<A<<" b2-1: "<<b2-1<<std::endl;

	for(i=A;i<b2;i++) { 
		//while_tot=0
		par_num=chld_num=0; for(i3=0;i3<cyc_len;i3++) snk_r[i3]=snk_c[i3]=par_r[i3]=par_c[i3]=-1; //refresh
		strt_r=edge_r[i]; strt_c=edge_c[i];
		//cout<<'\n'<<"start: "<<strt_r<<","<<strt_c;	
		find_cycl();
		Hsc[strt_r][strt_c]=0; //decimating unnecessary edges
		//cout<<'\n'<<"cycle length: "<<cyc_len<<" parent: "<<par_num<<" children: "<<chld_num;
		//cout<<'\n'<<"while_tot: "<<while_tot<<" rev_tot: "<<rev_tot;
		sum+=par_num+chld_num;
		//cout<<'\n'<<i<<std::endl;
		//if(count==30) break;
	}
	//cout<<'\n'<<"cycles: "<<ccnt;
	row_strt=row_stp=col_strt=col_stp=0;
	return ccnt;
	//cout<<'\n'<<"cycles: "<<'\n'; for(j=0;j<ccnt;j++) {for(i3=0;i3<cyc_len;i3++) cout<<cyc_sv_rg[j][i3]<<","<<cyc_sv_cg[j][i3]<<" "; cout<<'\n';} 

}
