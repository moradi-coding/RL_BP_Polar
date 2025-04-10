
//BCJR calclations using the [7,4] Hamming trellis in Vivian paper, calculations based on Ryan-Lin book
//flooding GLDPC decoder for all GCNs
void bcjr_log(int a) { 
		
	long i,j;
	double tmp1,tmp2,tmp3,tmp4,tmp1b,tmp2b,tmp3b,tmp4b,max1,max2,max3,max1b,max2b,max3b,c,y_sum,b,d,num;

	double gama1_0, gama1_1, gama2_0, gama2_1, gama3_0, gama3_1, gama4_0, gama4_1, gama5_0, gama5_1, gama6_0, gama6_1, gama7_0, gama7_1, gama8_0, gama8_1, gama9_0, gama9_1, gama10_0, gama10_1, gama11_0, gama11_1, gama12_0, gama12_1, gama13_0, gama13_1, gama14_0, gama14_1, gama15_0, gama15_1; 

	double alph, alph1_0, alph1_1, alph2_00, alph2_01, alph2_10, alph2_11,
	alph3_00, alph3_01, alph3_10, alph3_11, alph4_000, alph4_001, alph4_010, alph4_011, alph4_100, alph4_101, alph4_110, alph4_111, 
	alph5_00, alph5_01, alph5_10, alph5_11, 
	alph6_0, alph6_1, alph7_0;
	
	double beta1_0, beta1_1, beta2_00, beta2_01, beta2_10, beta2_11, 
	beta3_00, beta3_01, beta3_10, beta3_11, beta4_000, beta4_001, beta4_010, beta4_011, beta4_100, beta4_101, beta4_110, beta4_111, 
	beta5_00, beta5_01, beta5_10, beta5_11, 
	beta6_0, beta6_1, beta7_0;
	
	d=0.5; 
	
	c=-1/(2*varn);
	//for(a=0;a<num_cls;a++) {
		//initialization of branch metrics
		for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a', j is stage of trellis
			i=vns_cluster[a][j]; //i is a VN of cluster (CN) 'a'
			//cout<<"j, b, LLR: "<<j<<", "<<b<<", "<<L_in[i]<<endl;
			b=L_in[i]-LR[i]; //gives the extrinsic LLR
			
			if(!j) {
				gama1_0=th(c*pow(abs(y[i]-1),2)+d*1*b); //0 mapped to 1	
				gama1_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); //1 mapped to -1	
			}	
			else if(j==1) {
				gama2_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama2_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==2) {
				gama3_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama3_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==3) {
				gama4_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama4_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==4) {
				gama5_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama5_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==5) {
				gama6_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama6_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}	
			else if(j==6) {
				gama7_0=th(c*pow(abs(y[i]-1),2)+d*1*b); 
				gama7_1=th(c*pow(abs(y[i]+1),2)+d*-1*b); 
			}		
		}
		//cout<<"L_in: "; for(i=0;i<n;i++) cout<<L_in[i]<< " "; cout<<'\n'<<'\n';
		
		//forward metrics
		alph=log(1);
		
		alph1_0=th(alph+gama1_0); //threshold the result
		alph1_1=th(alph+gama1_1); 
		
		alph2_00=th(alph1_0+gama2_0);
		alph2_01=th(alph1_1+gama2_1);
		alph2_10=th(alph1_1+gama2_0);
		alph2_11=th(alph1_0+gama2_1);
		
		alph3_00=th2(alph2_00+gama3_0, alph2_01+gama3_1);
		alph3_01=th2(alph2_00+gama3_1, alph2_01+gama3_0);
		alph3_10=th2(alph2_10+gama3_0, alph2_11+gama3_1);
		alph3_11=th2(alph2_10+gama3_1, alph2_11+gama3_0);
		
		alph4_000=th(alph3_00+gama4_0);
		alph4_001=th(alph3_01+gama4_0);
		alph4_010=th(alph3_10+gama4_0);
		alph4_011=th(alph3_11+gama4_0);
		alph4_100=th(alph3_00+gama4_1);
		alph4_101=th(alph3_01+gama4_1);
		alph4_110=th(alph3_10+gama4_1);
		alph4_111=th(alph3_11+gama4_1);
		
		alph5_00=th2(alph4_000+gama5_0, alph4_110+gama5_1);
		alph5_01=th2(alph4_010+gama5_0, alph4_100+gama5_1);
		alph5_10=th2(alph4_011+gama5_1, alph4_101+gama5_0);
		alph5_11=th2(alph4_001+gama5_1, alph4_111+gama5_0);
		
		alph6_0=th2(alph5_00+gama6_0, alph5_11+gama6_1);
		alph6_1=th2(alph5_01+gama6_1, alph5_10+gama6_0);
		
		alph7_0=th2(alph6_0+gama7_0, alph6_1+gama7_1);
		
		
		//backward metrics
		beta7_0=log(1);
		
		beta6_0=th(beta7_0+gama7_0);
		beta6_1=th(beta7_0+gama7_1);
		
		beta5_00=th(beta6_0+gama6_0);
		beta5_01=th(beta6_1+gama6_1);
		beta5_10=th(beta6_1+gama6_0);
		beta5_11=th(beta6_0+gama6_1);
		
		beta4_000=th(beta5_00+gama5_0);
		beta4_001=th(beta5_11+gama5_1);
		beta4_010=th(beta5_01+gama5_0);
		beta4_011=th(beta5_10+gama5_1);
		beta4_100=th(beta5_01+gama5_1);
		beta4_101=th(beta5_10+gama5_0);
		beta4_110=th(beta5_00+gama5_1);
		beta4_111=th(beta5_11+gama5_0);
		
		beta3_00=th2(beta4_000+gama4_0, beta4_100+gama4_1);
		beta3_01=th2(beta4_001+gama4_0, beta4_101+gama4_1);
		beta3_10=th2(beta4_010+gama4_0, beta4_110+gama4_1);
		beta3_11=th2(beta4_011+gama4_0, beta4_111+gama4_1);
		
		beta2_00=th2(beta3_00+gama3_0, beta3_01+gama3_1);
		beta2_01=th2(beta3_00+gama3_1, beta3_01+gama3_0);
		beta2_10=th2(beta3_10+gama3_0, beta3_11+gama3_1);
		beta2_11=th2(beta3_10+gama3_1, beta3_11+gama3_0);
		
		beta1_0=th2(beta2_00+gama2_0, beta2_11+gama2_1);
		beta1_1=th2(beta2_01+gama2_1, beta2_10+gama2_0);
		

	for(j=0;j<num_vns_cls;j++) { //num_vns_cls is the no. of VNs connected to cluster 'a', j is stage of trellis
		i=vns_cluster[a][j]; //i is a VN of cluster (CN) 'a'
		if(!j) {
			//L_out[i]= log(exp(alph_0+gama1_0+s000_1_beta)) - log(exp(alph_0+gama1_1+s010_1_beta));
			tmp1=th(alph+gama1_0+beta1_0);
			tmp1b=th(alph+gama1_1+beta1_1);
			
			//L_out[i]= log(exp(tmp1)) - log(exp(tmp1b)); 
			//cout<<"j: "<<j<<'\n';
			//cout<<exp(tmp1)<<endl;
		}
		
		else if(j==1) {
			//L_out[i]= log(exp(s000_1_alpha+gama1_00+s000_2_beta)+exp(s010_1_alpha+gama1_10+s010_2_beta)) - log(exp(s000_1_alpha+gama1_01+s011_2_beta)+exp(s010_1_alpha+gama1_11+s001_2_beta));
			tmp1=th(alph1_0+gama2_0+beta2_00);
			tmp2=th(alph1_1+gama2_0+beta2_10);
			
			tmp1b=th(alph1_1+gama2_1+beta2_01);
			tmp2b=th(alph1_0+gama2_1+beta2_11);
			//L_out[i]= max1+log(1+exp(-1*abs(tmp1-tmp2))) - (max2+log(1+exp(-1*abs(tmp1b-tmp2b))));
		}
	
		else if(j==2) {
			//L_out[i]= log(exp(s000_2_alpha+gama2_00+s000_3_beta)+exp(s001_2_alpha+gama2_10+s001_3_beta)+exp(s010_2_alpha+gama2_00b+s010_3_beta)+exp(s011_2_alpha+gama2_10b+s011_3_beta)) - log(exp(s000_2_alpha+gama2_01+s001_3_beta)+exp(s001_2_alpha+gama2_11+s000_3_beta)+exp(s010_2_alpha+gama2_01b+s011_3_beta)+exp(s011_2_alpha+gama2_11b+s010_3_beta));
			
			tmp1=th(alph2_00+gama3_0+beta3_00);
			tmp2=th(alph2_01+gama3_0+beta3_01);
			tmp3=th(alph2_10+gama3_0+beta3_10);
			tmp4=th(alph2_11+gama3_0+beta3_11);
			
			tmp1b=th(alph2_00+gama3_1+beta3_01);
			tmp2b=th(alph2_01+gama3_1+beta3_00);
			tmp3b=th(alph2_10+gama3_1+beta3_11);
			tmp4b=th(alph2_11+gama3_1+beta3_10);
			
			//L_out[i]= log(exp(tmp1)+exp(tmp2)+exp(tmp3)+exp(tmp4)) - log(exp(tmp1b)+exp(tmp2b)+exp(tmp3b)+exp(tmp4b));
			//max(a,b,c,d)=max((a,b,c),d)=max(max2(max1(a,b),c)),d)=max(max2,d)=max(max2,tmp4)=max3+log(1+exp(-1*abs(max2-tmp4)))
		}
	
		else if(j==3) {
			tmp1=th(alph3_00+gama4_0+beta4_000);
			tmp2=th(alph3_01+gama4_0+beta4_001);
			tmp3=th(alph3_10+gama4_0+beta4_010);
			tmp4=th(alph3_11+gama4_0+beta4_011);
			
			tmp1b=th(alph3_00+gama4_1+beta4_100);
			tmp2b=th(alph3_01+gama4_1+beta4_101);
			tmp3b=th(alph3_10+gama4_1+beta4_110);
			tmp4b=th(alph3_11+gama4_1+beta4_111);
		}
	
		else if(j==4) {
			tmp1=th(alph4_000+gama5_0+beta5_00);
			tmp2=th(alph4_010+gama5_0+beta5_01);
			tmp3=th(alph4_101+gama5_0+beta5_10);
			tmp4=th(alph4_111+gama5_0+beta5_11);
			
			tmp1b=th(alph4_001+gama5_1+beta5_11);
			tmp2b=th(alph4_011+gama5_1+beta5_10);
			tmp3b=th(alph4_100+gama5_1+beta5_01);
			tmp4b=th(alph4_110+gama5_1+beta5_00);
		}
	
		else if(j==5) {
			tmp1=th(alph5_00+gama6_0+beta6_0);
			tmp2=th(alph5_10+gama6_0+beta6_1);
			
			tmp1b=th(alph5_01+gama6_1+beta6_1);
			tmp2b=th(alph5_11+gama6_1+beta6_0);
		}

		else if(j==6) {
			tmp1=th(alph6_0+gama7_0+beta7_0);
			tmp1b=th(alph6_1+gama7_1+beta7_0);
		}
		
		if(!j || j==6) 
			L_out[i]= th(tmp1+log(1+exp(-1*abs(tmp1))) - (tmp1b+log(1+exp(-1*abs(tmp1b)))));
		else if(j==1 || j==5) 
			L_out[i]= th(maxx(tmp1,tmp2)+log(1+exp(-1*abs(tmp1-tmp2))) - (maxx(tmp1b,tmp2b)+log(1+exp(-1*abs(tmp1b-tmp2b)))));
		else {
			max2=maxx(maxx(tmp1,tmp2),tmp3);	
			max3=maxx(max2,tmp4);	
			max2b=maxx(maxx(tmp1b,tmp2b),tmp3b);	
			max3b=maxx(max2b,tmp4b);
			L_out[i]=th(max3+log(1+exp(-1*abs(max2-tmp4))) - (max3b+log(1+exp(-1*abs(max2b-tmp4b)))));
		}
	
	}
	//cout<<"a: "<<a<<endl;
	//cout<<'\n'<<"L_out: "<<endl;
	for(i=0;i<n;i++)
		if(isnan(L_out[i])) 
			cout<<L_out[i]<<" ";
	
	
}

